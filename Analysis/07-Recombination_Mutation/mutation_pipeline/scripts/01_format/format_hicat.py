"""
Snakemake 执行脚本: 格式化HiCAT数据

本脚本由Snakemake规则调用，用于处理单个样本。
它从snakemake对象获取输入/输出路径，并调用共享库中的核心功能函数。
"""
import sys
import os
sys.path.append(os.getcwd())
from workflow.lib.hor_preprocess.hicat import load_hicat_hor

# --- Snakemake API ---
# 日志重定向: 将所有print输出写入Snakemake指定的日志文件
log_file = snakemake.log[0]
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

# 从snakemake对象获取信息
sample_name = snakemake.wildcards.sample
hicat_dir = snakemake.input.hicat_dir
output_dir = snakemake.output[0] # Snakemake保证这个目录存在

print(f"--- Starting HiCAT formatting for sample: {sample_name} ---")
print(f"Input HiCAT dir: {hicat_dir}")
print(f"Output directory: {output_dir}")

# --- 执行核心功能 ---
try:
    load_hicat_hor(hicat_dir, output_dir)
    print(f"--- Successfully completed HiCAT formatting for sample: {sample_name} ---")
except Exception as e:
    print(f"ERROR: An exception occurred while processing sample {sample_name}: {e}")
    # 抛出异常以让Snakemake知道任务失败
    raise e