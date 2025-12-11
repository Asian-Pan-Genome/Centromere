"""
Snakemake执行脚本：执行样本对之间的HORSCAN比对。
它会自动查找两个输入目录中的共同染色体，并为每一对执行HORSCAN。
"""

import sys
import subprocess
from pathlib import Path

# --- 1. 从snakemake对象中获取所有信息 ---
# 当Snakemake通过`script:`指令调用此文件时，
# 它会自动创建一个名为`snakemake`的全局对象。
# 我们可以通过这个对象的属性来访问输入、输出、通配符等。

# 输入和输出路径
input_dir_s = Path(snakemake.input.dir_s)
input_dir_t = Path(snakemake.input.dir_t)
output_dir = Path(snakemake.output[0])

# 日志文件路径
log_file = Path(snakemake.log[0])

# 通配符和参数
wildcards = snakemake.wildcards
params = snakemake.params

# --- 日志重定向 ---
# 将所有print输出和错误都写入Snakemake指定的日志文件
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

print(f"--- Starting Pairwise Alignment for {wildcards.sample_s} vs {wildcards.sample_t} ---")

# --- 2. 确保输出目录存在 ---
output_dir.mkdir(parents=True, exist_ok=True)

# --- 3. 动态查找共有的染色体 ---
chroms_s = {f.stem.split('.')[0] for f in input_dir_s.glob("*.mon.bed")}
chroms_t = {f.stem.split('.')[0] for f in input_dir_t.glob("*.mon.bed")}
common_chroms = sorted(list(chroms_s.intersection(chroms_t)))

if not common_chroms:
    print("Warning: No common chromosomes found. Exiting successfully.")
    sys.exit(0)  # 正常退出

print(f"Found {len(common_chroms)} common chromosomes to align: {common_chroms}")

# --- 4. 循环执行HORSCAN命令 ---
for chrom in common_chroms:
    print(f"\n  -- Aligning chromosome: {chrom} --")
    
    input_s_file = input_dir_s / f"{chrom}.mon.bed"
    input_t_file = input_dir_t / f"{chrom}.mon.bed"
    output_prefix = output_dir / f"{chrom}"
    
    command = (
        f"{params.executable} "
        f"-s {input_s_file} "
        f"-t {input_t_file} "
        f"-o {output_prefix} "
        f"--mode {params.mode}"
    )
    
    print(f"Executing command: {command}")
    
    try:
        # 执行命令并捕获输出
        result = subprocess.run(
            command, shell=True, check=True, text=True, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        print(f"  -- Successfully aligned {chrom}. --")
        # 打印HORSCAN本身的输出到日志
        if result.stdout: print(f"Stdout:\n{result.stdout}")
        if result.stderr: print(f"Stderr:\n{result.stderr}")
            
    except subprocess.CalledProcessError as e:
        print(f"ERROR: HORSCAN failed for chromosome {chrom}.")
        print(f"Return code: {e.returncode}")
        print(f"Stdout:\n{e.stdout}")
        print(f"Stderr:\n{e.stderr}")
        # 抛出异常，让Snakemake知道此任务失败
        raise e

print(f"\n--- All alignments for {wildcards.sample_s} vs {wildcards.sample_t} complete. ---")