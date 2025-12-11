"""
此脚本由 Snakemake 的 'monomer_mutation_calculate' 规则调用。

功能:
1.  输入指定的筛选样本list，并找到对应的BLAST和FASTA文件。
2.  对每个样本对，按染色体进行处理。
3.  对于 MTH/MIS 类型的比对，使用edlib计算精确的cigar, snp, indel，并处理不等长序列。
4.  对于 INS/DEL 类型的比对，直接计算indel并生成简化cigar。
5.  将增强后的比对结果，为每个样本对的每个染色体单独保存为一个gzipped文件。
"""

import os
import sys
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import edlib
import re
from tqdm import tqdm
from typing import Optional, List

# 假设您的自定义模块位于 workflow/lib/ 目录下
# 您需要确保这些模块可以被正确导入
sys.path.append(os.getcwd())

from workflow.lib.blast.loader import load_blast_gz
from workflow.lib.fasta import MapFasta


# ===================================================================
# 1. 核心功能函数 (包含完整的不等长处理逻辑)
# ===================================================================

def parse_cigar_snp(cigar: str):
    """从CIGAR字符串中解析SNP和Indel计数。"""
    if not cigar: return 0, 0
    operations = re.findall(r'(\d+)([=XIDSM])', cigar)
    snp_count = sum(int(length) for length, op in operations if op == 'X')
    indel_count = sum(int(length) for length, op in operations if op in ('I', 'D', 'S'))
    return snp_count, indel_count

def calculate_blast_cigar(blast_df_for_chrom: pd.DataFrame, source_fa: MapFasta, target_fa: MapFasta) -> Optional[pd.DataFrame]:
    """
    为单个染色体的BLAST结果DataFrame添加cigar, snp, indel信息。
    内部区分处理 MTH/MIS 和 INS/DEL，并包含不等长序列的处理逻辑。
    """
    if blast_df_for_chrom.empty:
        return None

    # a) 将DataFrame按类型拆分
    mth_mis_df = blast_df_for_chrom[blast_df_for_chrom['type'].isin(['MTH', 'MIS'])].copy()
    ins_del_df = blast_df_for_chrom[blast_df_for_chrom['type'].isin(['INS', 'DEL'])].copy()
    
    processed_dfs = []

    # b) 处理 MTH/MIS (需要序列比对)
    if not mth_mis_df.empty:
        results = []
        for idx, row in mth_mis_df.iterrows():
            try:
                source_seq = source_fa[row['chr']][row['q.start']:row['q.end']]
                target_seq = target_fa[row['chr']][row['s.start']:row['s.end']]
                
                # --- 完整的不等长处理逻辑 ---
                if len(source_seq) == len(target_seq):
                    if len(source_seq) == 0: continue
                    alignment = edlib.align(source_seq, target_seq, task="path", mode="NW")
                    cigar = alignment['cigar']
                else:
                    equal_length = min(len(source_seq), len(target_seq))
                    if equal_length == 0: continue
                    
                    # 前N个碱基比对
                    alignment_front = edlib.align(source_seq[:equal_length], target_seq[:equal_length], task="path", mode="NW")
                    front_cigar = alignment_front['cigar']
                    front_snp, front_indel = parse_cigar_snp(front_cigar)
                    
                    # 后N个碱基比对
                    alignment_after = edlib.align(source_seq[-equal_length:], target_seq[-equal_length:], task="path", mode="NW")
                    after_cigar = alignment_after['cigar']
                    after_snp, after_indel = parse_cigar_snp(after_cigar)
                    
                    # 选择差异更小的一个
                    if (front_snp + front_indel) < (after_snp + after_indel):
                        cigar = front_cigar
                    else:
                        cigar = after_cigar
                # --- 不等长处理逻辑结束 ---

                snp, indel = parse_cigar_snp(cigar)
                row['cigar'], row['snp'], row['indel'] = cigar, snp, indel
                results.append(row)

            except (KeyError, TypeError):
                # print(f"  [WARNING] 在FASTA文件中找不到染色体 '{row['chr']}' 或坐标无效。跳过此行。")
                continue
            except Exception as e:
                # print(f"  [ERROR] 在计算CIGAR时发生错误: {e}。跳过此行。")
                continue
        if results:
            processed_dfs.append(pd.DataFrame(results))

    # c) 处理 INS/DEL (直接计算)
    if not ins_del_df.empty:
        ins_del_df['q.len'] = ins_del_df['q.end'] - ins_del_df['q.start']
        ins_del_df['s.len'] = ins_del_df['s.end'] - ins_del_df['s.start']
        ins_del_df['indel'] = (ins_del_df['q.len'] - ins_del_df['s.len']).abs()
        ins_del_df['snp'] = 0
        ins_del_df['cigar'] = ins_del_df.apply(
            lambda row: f"{row['indel']}I" if row['type'] == 'INS' else f"{row['indel']}D", axis=1)
        processed_dfs.append(ins_del_df)

    # d) 合并并返回
    if not processed_dfs:
        return None
    
    final_chrom_df = pd.concat(processed_dfs, ignore_index=True)
    final_chrom_df.sort_values(by=['q.start', 's.start'], inplace=True)
    final_chrom_df['q.len'] = final_chrom_df['q.end'] - final_chrom_df['q.start']
    final_chrom_df['s.len'] = final_chrom_df['s.end'] - final_chrom_df['s.start']
    return final_chrom_df

# ===================================================================
# 2. 并行处理的 "工人" 函数
# ===================================================================
def process_pair_group(pair, group_df: pd.DataFrame, blast_dir: str, fasta_map: dict, output_root_dir: str) -> List[str]:
    source, target = pair
    created_files = []
    try:
        source_fa_path = fasta_map[source]
        target_fa_path = fasta_map[target]
    except KeyError as e:
        print(f"  [ERROR] 在样本清单中找不到样本 {e} 的FASTA路径。跳过配对 ({source}, {target}).")
        return created_files

    try:
        source_fa = MapFasta(source_fa_path)
        target_fa = MapFasta(target_fa_path)
        
        blast_gz_path = os.path.join(blast_dir, source, f"{source}.{target}.blast.gz")
        
        if not blast_gz_path: return created_files
        
        blast_df_full = load_blast_gz(blast_gz_path)
        if blast_df_full.empty: return created_files
        
        for chrom in group_df['chrom'].unique():
            try:
                blast_df_filtered = blast_df_full[
                    (blast_df_full['chr'] == chrom)
                ].copy()
                
                if blast_df_filtered.empty: continue
                
                processed_chrom_df = calculate_blast_cigar(blast_df_filtered, source_fa, target_fa)
                
                if processed_chrom_df is not None and not processed_chrom_df.empty:
                    output_subdir = os.path.join(output_root_dir, source)
                    os.makedirs(output_subdir, exist_ok=True)
                    output_path = os.path.join(output_subdir, f"{source}.{target}.{chrom}.blast.cigar.gz")
                    processed_chrom_df.to_csv(output_path, sep='\t', index=False, compression='gzip')
                    created_files.append(output_path)
            except Exception as e:
                # print(f"  [ERROR] Error processing chrom {chrom} for pair ({source}, {target}): {e}")
                continue
    except Exception as e:
        # print(f"  [CRITICAL ERROR] An unexpected error occurred for pair ({source}, {target}): {e}")
        return created_files

    return created_files

# ===================================================================
# 3. 主逻辑 (由Snakemake执行)
# ===================================================================
def main():
    # a) 从 snakemake 对象中读取参数
    filtered_sample_file = snakemake.input.filtered_sample
    output_dir = snakemake.output.monomer_snp_dir
    sample_sheet_file = snakemake.params.sample_bed
    blast_dir = snakemake.input.blast_dir
    threads = snakemake.threads
    
    # b) 准备数据
    sample_sheet = pd.read_csv(sample_sheet_file, sep='\t', comment='#') # 假设有 'sample' 和 'fasta_path' 列
    fasta_map = sample_sheet.set_index('sample')['assembly_fa'].to_dict()
    pass_df = pd.read_csv(filtered_sample_file, sep='\t', comment='#') # 假设列名为 'source', 'target', 'chrom'
    
    # c) 按 (source, target) 分组
    grouped_pairs = pass_df.groupby(['source', 'target'])
    
    os.makedirs(output_dir, exist_ok=True)
    
    # d) 设置并行任务
    tasks = []
    for pair, group_df in grouped_pairs:
        source, target = pair
        if source in fasta_map and target in fasta_map:
             tasks.append(delayed(process_pair_group)(pair, group_df, blast_dir, fasta_map, output_dir))
    
    # e) 并行执行并显示进度条
    print(f"Starting mutation calculation for {len(tasks)} sample pairs on {threads} threads...")
    results_list_of_lists = Parallel(n_jobs=threads)(
        tqdm(tasks, total=len(tasks), desc="Calculating Mutations")
    )
    
    # f) 报告结果
    all_created_files = [item for sublist in results_list_of_lists for item in sublist]
    print(f"Processing complete. A total of {len(all_created_files)} files were created in {output_dir}")
    with open(os.path.join(output_dir, "manifest.txt"), "w") as f:
        for path in sorted(all_created_files):
            f.write(f"{path}\n")

if __name__ == '__main__':
    # --- Snakemake 执行入口 ---
    sys.stderr = sys.stdout = open(snakemake.log[0], 'w')
    try:
        main()
    except Exception as e:
        import traceback
        print(f"A critical error occurred in the main script execution: {e}")
        traceback.print_exc()
        sys.exit(1)