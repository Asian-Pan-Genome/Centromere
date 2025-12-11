# workflow/rules/horscan.smk


rule horscan_pairwise:
    """
    对一对比对样本（sample_s vs sample_t），自动查找其共有的染色体文件，
    并为每一个共同的染色体执行HORSCAN比对。
    """
    input:
        dir_s="results/APG_v1_analysis/sample_data/asat_bed/{sample_s}",
        dir_t="results/APG_v1_analysis/sample_data/asat_bed/{sample_t}",
    output:
        directory("results/APG_v1_analysis/horscan_pairwise/{sample_s}/{sample_t}"),
    log:
        "logs/horscan_pairwise/APG_v1_analysis/{sample_s}/{sample_t}.log",
    params:
        executable=config["tools"]["horscan_v2"],
        mode=config["tools"]["horscan_mode"],
    conda:
        "APG"
    group: 
        "horscan"  
    script:
        "../../scripts/02_horscan/pairwise_horscan.py"


rule merge_horscan_results:
    """
    合并样本对的所有HORSCAN比对结果。
    脚本内部会智能判断是否需要执行。
    """
    input:
        # ✅ 依赖关系被简化：总是、无条件地依赖上游的输出目录
        alignment_dir="results/APG_v1_analysis/horscan_pairwise/{sample_s}/{sample_t}"
    output:
        # 输出是最终的blast.gz文件
        config["blast_output_dir"] + "/{sample_s}/{sample_s}.{sample_t}.blast.gz"
    log:
        "logs/horscan_pairwise/APG_v1_analysis/{sample_s}_vs_{sample_t}.merge.log"
    conda:
        "APG"
    group: 
        "merge"
    script:
        "../../scripts/02_horscan/merge_horscan.py"

rule horscan_check:
    """
    检查所有样本对的HORSCAN比对结果是否已完成。
    """
    input:
        # 使用 expand 和 zip 来为我们生成的所有样本对创建目标文件路径
       expand(
            config["blast_output_dir"] + "/{sample_s}/{sample_s}.{sample_t}.blast.gz",
            zip,
            sample_s=SOURCE_SAMPLES,
            sample_t=TARGET_SAMPLES,
        ),
    output:
        touch("results/APG_v1_analysis/flags/horscan_check.done"),

# rule horscan_test:
#     input:
#         # 使用 expand 和 zip 来为我们生成的所有样本对创建目标文件路径
#         # "results/APG_v1_analysis/NA20129_Mat/NA20129_Mat.NA20355_hap2.blast.gz.done",
#         expand(
#             config["blast_output_dir"] + "/NA20129_Mat/NA20129_Mat.NA24385_hap1.blast.gz",
#         ),
#     output:
#         touch("results/APG_v1_analysis/flags/horscan_stage_test.txt"),