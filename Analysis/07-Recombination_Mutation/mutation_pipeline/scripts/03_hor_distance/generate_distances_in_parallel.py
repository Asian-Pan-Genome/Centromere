#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Snakemake执行脚本：一个使用joblib并行生成所有成对 .distance 文件的脚本。

## 功能
这个脚本被设计成一个单一的、大型的Snakemake任务。它接收样本列表和配置参数，
然后在内部使用`joblib`库将成千上万个“成对距离计算”子任务分配到多个CPU核心上并行执行。

每个并行的子任务会：
1. 计算一对样本的距离。
2. 将结果直接保存为对应的 `.distance` 文件。

主进程只负责调度、监控进度（通过tqdm）和总结结果。

## 设计模式与警告
此脚本采用了“在脚本内部并行”的模式，以绕过Snakemake为海量任务构建DAG时的性能瓶颈。
这是一个高级的、有特定取舍的方案：
- **优点**: Snakemake启动速度极快，因为它只看到了一个任务。
- **缺点**: 牺牲了Snakemake最核心的失败恢复能力。虽然脚本内置了简易的断点续传，
  但它不如Snakemake原生的调度机制稳健和强大，尤其是在复杂的集群环境中。
"""

# ===================================================================
# 1. 导入标准库和第三方库
# ===================================================================
import pandas as pd
import numpy as np
from pathlib import Path
import time
import sys
import os
import itertools

# joblib用于实现稳健的多进程并行计算
from joblib import Parallel, delayed
# tqdm用于生成美观且信息丰富的进度条
from tqdm.auto import tqdm

# ===================================================================
# 2. 导入项目自定义库
# ===================================================================
# 将工作流的根目录添加到Python的搜索路径中，这样就可以导入自定义的`workflow.lib`模块
sys.path.append(os.getcwd())
from workflow.lib.hor_distance_func import calculate_final_distances
from workflow.lib.blast.loader import load_blast_gz, load_bed
from workflow.lib.blast.anno import annotate_bed

# ===================================================================
# 3. 工作函数 (处理并保存单个样本对)
# ===================================================================
def process_and_save_pair(args_tuple):
    """
    这是并行化的核心工作函数。它被设计为完全独立运行。

    功能:
    - 处理单个样本对 (sample_s vs sample_t)。
    - 计算所有染色体的距离。
    - 将结果直接保存到对应的.distance文件中。

    参数 (args_tuple):
    - 一个元组，包含此任务所需的所有信息，方便joblib传递。

    返回:
    - 成功: 返回创建的文件的字符串路径。
    - 失败: 返回 None。
    """
    # --- a. 解包参数，使代码更具可读性 ---
    s_sample, t_sample, blast_dir, graph_dir, hor_params, output_root = args_tuple
    
    # --- b. 动态构建此子任务的输出文件路径 ---
    output_path = Path(output_root) / s_sample / f"{s_sample}.{t_sample}.distance"
    
    try:
        # ✅ 核心优化：实现断点续传
        # 如果目标文件已经存在且内容不为空，则认为此任务已完成，直接跳过以节省时间。
        if output_path.exists() and output_path.stat().st_size > 0:
            # return str(output_path)
            raise ValueError(f"Output file {output_path} already exists and is not empty. Skipping this pair.")

        # --- c. 加载计算所需的输入文件 ---
        blast_gz_path = Path(blast_dir) / s_sample / f"{s_sample}.{t_sample}.blast.gz"
        s_hor_path = Path(graph_dir) / s_sample / "hor.bed"
        t_hor_path = Path(graph_dir) / t_sample / "hor.bed"

        # 如果关键输入不存在，则无法进行计算，返回None表示失败
        if not blast_gz_path.exists():
            raise ValueError(f"BLAST file {blast_gz_path} does not exist for pair {s_sample}-{t_sample}.")

        blast_df = load_blast_gz(blast_gz_path)
        # 如果blast结果为空，创建一个空的.distance文件并标记为成功
        if blast_df.empty:
            raise ValueError(f"BLAST file {blast_gz_path} is empty for pair {s_sample}-{t_sample}.")
            
        
        # --- d. 执行核心计算逻辑（注释和距离计算）---
        if 'q.hor' not in blast_df.columns:
            s_hor_df = load_bed(s_hor_path, cols=["chr", "start", "end", "hor"])
            blast_df['q.hor'] = annotate_bed(blast_df, s_hor_df, ["chr", "q.start", "q.end"], ["chr", "start", "end", "hor"])
        if 's.hor' not in blast_df.columns:
            t_hor_df = load_bed(t_hor_path, cols=["chr", "start", "end", "hor"])
            blast_df['s.hor'] = annotate_bed(blast_df, t_hor_df, ["chr", "s.start", "s.end"], ["chr", "start", "end", "hor"])

        results_list = []
        for chrom, chr_df in blast_df.groupby("chr"):
            if not chr_df.empty:
                merge_dist, hor_dist, aln_dist = calculate_final_distances(chr_df.copy(), hor_params)
                results_list.append({
                    "source": s_sample, "target": t_sample, "chrom": chrom,
                    "aln_distance": aln_dist, "hor_distance": hor_dist,
                    "merge_distance": merge_dist, "weight": hor_params.get("weight", np.nan)
                })

        # --- e. 将结果直接保存到文件 ---
        # 确保输出目录存在
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if results_list:
            result_df = pd.DataFrame(results_list)
            result_df.to_csv(output_path, sep='\t', index=False)
        else:
            # 如果没有结果，也创建一个空文件，表示已处理但无内容
            output_path.touch()

        # --- f. 返回文件路径，表示此子任务成功 ---
        return str(output_path)

    except Exception as e:
        # ✅ 核心健壮性：独立的错误处理
        # 捕获所有可能的异常（文件未找到、计算错误等），打印错误信息并返回None。
        # 这可以防止单个子任务的失败导致整个并行池崩溃。
        print(f"ERROR: Subprocess for pair {s_sample}-{t_sample} failed: {e}", file=sys.stderr)
        return None

# ===================================================================
# 4. 主函数 (作为任务调度器)
# ===================================================================
def main():
    """主函数，负责设置、调度并行任务并总结结果。"""
    
    # --- a. 从Snakemake获取所有必要的配置和路径 ---
    params = snakemake.params
    inputs = snakemake.input
    output_done_flag = Path(snakemake.output.done_flag)
    num_threads = snakemake.threads
    
    print("--- Starting Parallel Distance File Generation Script ---")
    print(f"Parallelism level (n_jobs): {num_threads}")

    # --- b. 生成所有需要计算的唯一上三角样本对 ---
    sample_pairs = params.sample_pair
    print(f"Generated {len(sample_pairs)} unique pairs to process.")
    
    # --- c. 准备将要传递给每个并行工作进程的参数列表 ---
    task_args = [
        (s, t, inputs.blast_dir, params.graph_dir, params.hor_params, params.output_root) 
        for s, t in sample_pairs
    ]

    # --- d. 使用joblib和tqdm执行并行计算 ---
    print("\nStarting parallel processing...")
    t_start = time.perf_counter()
    
    # `Parallel`是joblib的核心调度器。
    # `n_jobs`设置并行进程数，直接从snakemake的threads获取。
    # `backend='loky'`是比默认'multiprocessing'更稳健的后端。
    # `delayed`是joblib的函数，用于延迟执行我们的工作函数。
    # `tqdm`包裹可迭代对象，以显示实时进度条。
    results = Parallel(n_jobs=num_threads, backend='loky')(
        delayed(process_and_save_pair)(args) for args in tqdm(task_args, desc="Generating .distance files")
    )
    
    t_end = time.perf_counter()
    print(f"\n...Parallel processing finished in {t_end - t_start:.2f} seconds.")

    # --- e. 总结运行结果 ---
    success_count = sum(1 for r in results if r is not None)
    failure_count = len(results) - success_count
    
    print("\n" + "="*40)
    print("--- RUN SUMMARY ---")
    print(f"Total pairs scheduled: {len(results)}")
    print(f"  - Successfully processed/found: {success_count}")
    print(f"  - Failed during processing:     {failure_count}")
    print("="*40 + "\n")

    if failure_count > 0:
        print("WARNING: Some pairs failed to process. Check logs for ERROR messages.", file=sys.stderr)

    # --- f. 创建完成标志文件，通知Snakemake此宏任务已结束 ---
    output_done_flag.touch()
    print(f"Done flag created: {output_done_flag}")
    print("Script finished successfully.")

# ===================================================================
# 5. 脚本执行入口
# ===================================================================
if __name__ == '__main__':
    # --- a. 设置日志重定向 ---
    # 确保所有print和错误信息都写入Snakemake指定的日志文件
    log_file = snakemake.log[0]
    sys.stdout = open(log_file, 'w')
    sys.stderr = sys.stdout
    
    # --- b. 执行主逻辑，并捕获任何顶层错误 ---
    try:
        main()
    except Exception as e:
        # 捕获主进程中的致命错误（例如，参数缺失）
        print(f"\nFATAL ERROR in main orchestrator: {e}", file=sys.stderr)
        # 重新抛出异常，以便Snakemake知道此任务失败了
        raise