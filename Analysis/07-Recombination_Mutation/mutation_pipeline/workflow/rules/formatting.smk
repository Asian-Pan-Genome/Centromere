# workflow/rules/formatting.smk

rule format_hicat:
    """
    为每个样本，将其HiCAT原始输出，格式化为按染色体分割的BED-like CSV文件。
    """
    input:
        hicat_dir = lambda wildcards: SAMPLES.loc[wildcards.sample, "HiCAT_dir"],
    output:
        directory("results/{run_name}/sample_data/hicat/{sample}")
    log:
        "logs/format_hicat/{run_name}/{sample}.log"
    conda:
        "APG"
    script:
        "../../scripts/01_format/format_hicat.py"

rule format_graph:
    """
    为每个样本，将其Graph原始输出，格式化为按染色体分割的BED-like CSV文件
    """
    input:
        graph_dir = lambda wildcards: SAMPLES.loc[wildcards.sample, "graph_dir"],
    output:
        directory("results/{run_name}/sample_data/graph/{sample}")
    log:
        "logs/format_graph/{run_name}/{sample}.log"
    conda:
        "APG"
    group: 
        "APG"
    script:
        "../../scripts/01_format/format_graph.py"

rule format_asat_bed:
    """
    为每个样本，将其ASAT原始输出，格式化为按染色体分割的BED-like CSV文件
    """
    input:
        asat_bed = lambda wildcards: SAMPLES.loc[wildcards.sample, "mon_bed"],
    output:
        directory("results/{run_name}/sample_data/asat_bed/{sample}")
    log:
        "logs/format_asat_bed/{run_name}/{sample}.log"
    conda:
        "APG"
    group: 
        "APG"
    script:
        "../../scripts/01_format/format_asat_bed.py"

rule format_all:
    """
    统一调用所有格式化规则，确保每个样本的所有数据都被正确处理。
    """
    input:
        expand("results/{run_name}/sample_data/hicat/{sample}",
               run_name=config["run_name"],
               sample=SAMPLES.index),
        expand("results/{run_name}/sample_data/graph/{sample}",
               run_name=config["run_name"],
               sample=SAMPLES.index),
        expand("results/{run_name}/sample_data/asat_bed/{sample}",
               run_name=config["run_name"],
               sample=SAMPLES.index)
    output:
        touch("results/{run_name}/flags/format_all.done")