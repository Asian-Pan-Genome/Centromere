#!/usr/bin/env python
"""
Snakemake执行脚本：合并所有样本对的距离计算结果，并为多种距离指标
分别生成距离矩阵。

**版本 2.0 - 改进版**
- **正确性**: 确保在成功执行后创建 Snakemake 的 .done 标志文件。
- **性能**: 使用 groupby 替代嵌套循环，避免多次扫描大型DataFrame。
- **健壮性**: 更好的错误处理和日志记录，并使用命名的Snakemake输出。
- **可读性**: 将矩阵处理逻辑封装到辅助函数中。
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List

# ===================================================================
# 1. Snakemake 和日志记录设置
# ===================================================================
try:
    # ✅ 最佳实践: 使用命名的输入/输出/参数，使代码更清晰
    distance_dir: Path = Path(snakemake.params.distance_dir)
    output_dir: Path = Path(snakemake.output.dist_dir)
    done_flag: Path = Path(snakemake.output.done_flag)
    log_file: Path = Path(snakemake.log[0])
    sample_pair = snakemake.params.sample_pair
    all_samples: List[str] = snakemake.params.sample_order
except NameError:
    sys.exit("This script is meant to be run via Snakemake with named I/O.")

# --- 日志重定向 ---
# 使用'a'模式追加日志，而不是'w'覆盖，这样可以保留Snakemake的初始日志
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

# ===================================================================
# 2. 辅助函数
# ===================================================================


def process_and_save_matrix(
    data: pd.DataFrame,
    metric: str,
    samples: List[str],
    out_dir: Path,
    chrom: str
):
    """为一个染色体和一个指标创建、对称化并保存距离矩阵。"""

    # a. 使用pivot_table创建稀疏矩阵
    matrix = data.pivot_table(index='source', columns='target', values=metric)

    # b. 重新索引以确保是完整的 N x N 矩阵，并按指定顺序排列
    matrix = matrix.reindex(index=samples, columns=samples)

    # c. 对称化矩阵 (Aij = Aji)
    matrix = matrix.fillna(matrix.T)

    # d. 最终填充：对角线为0，其余NaN（如果仍存在）也填充为0或一个合适的默认值
    np.fill_diagonal(matrix.values, 0)
    matrix.fillna(0, inplace=True)

    # e. 保存矩阵
    output_path = out_dir / f"{metric}_{chrom}.csv"
    matrix.to_csv(output_path)
    print(f"    - Saved matrix to {output_path}")

# ===================================================================
# 3. 主执行逻辑
# ===================================================================


def main():
    print(f"--- Starting Distance Matrix Generation ---")
    print(f"Output directory: {output_dir}")
    print(f"Found {len(all_samples)} samples in specified order.")

    # --- 确保输出目录存在 ---
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- 读取并合并所有距离文件 ---
    distance_files = []
    for source, target in sample_pair:
        distance_files.append(distance_dir / source /
                              f"{source}.{target}.distance")

    print(f"Reading {len(distance_files)} input distance files...")
    # 若 f 不存在报告该路径不存在即可
    non_exist_files = [f for f in distance_files if not f.exists()]
    distance_files = [f for f in distance_files if f.exists()]
    if non_exist_files:
        print(f"Warning: The following distance files do not exist: {non_exist_files}")
    df_list = [pd.read_csv(f, sep='\t')
               for f in distance_files if f.stat().st_size > 0]

    if not df_list:
        print("Warning: No valid/non-empty distance files found. No matrices will be generated.")
        # ✅ 关键修复: 即使没有文件，也要创建完成标志，表示“成功处理了0个文件”
        done_flag.touch()
        print(f"Created empty done flag at {done_flag}. Exiting.")
        sys.exit(0)

    long_df = pd.concat(df_list, ignore_index=True)
    print(f"Total records loaded: {len(long_df)}")

    # --- 为每个指标的每个染色体创建矩阵 ---
    distance_metrics = ['aln_distance', 'hor_distance', 'merge_distance']

    # ✅ 性能优化: 使用 groupby('chrom')，只遍历一次DataFrame
    # 而不是为每个染色体反复过滤整个大型DataFrame。
    all_chromosomes = sorted(long_df['chrom'].unique())
    print(f"Processing {len(all_chromosomes)} chromosomes...")

    for chrom, chrom_df in long_df.groupby('chrom'):
        print(f"\n- Processing chromosome: {chrom}")
        for metric in distance_metrics:
            if metric not in chrom_df.columns:
                print(
                    f"  - Warning: Metric '{metric}' not found for chrom '{chrom}'. Skipping.")
                continue

            print(f"  - Generating matrix for metric: {metric}")
            print(f" - Miss value : {chrom_df[metric].isna().sum()}")
            process_and_save_matrix(
                chrom_df, metric, all_samples, output_dir, chrom)

    print("\n--- All distance matrices have been generated successfully. ---")

    # ✅ 关键修复: 在脚本成功完成所有工作后，创建.done标志文件
    # 这是告诉Snakemake任务已完成且所有输出都已就位的信号。
    done_flag.touch()
    print(f"Successfully created done flag: {done_flag}")


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        # 将任何未捕获的异常记录到日志文件中
        print(
            f"\nFATAL ERROR: A top-level exception occurred: {e}", file=sys.stderr)
        # 抛出异常，以便Snakemake知道任务失败了
        raise
