# rule pairwise_distance:
#     """
#     对每一对样本进行HOR distance的计算，最终分染色体保存
#     """
#     input:
#         blast_gz= ancient(config["blast_output_dir"] + "/{sample_s}/{sample_s}.{sample_t}.blast.gz"),
#         s_hor="results/APG_v1_analysis/sample_data/graph/{sample_s}",
#         t_hor="results/APG_v1_analysis/sample_data/graph/{sample_t}",
#     output:
#         "results/APG_v1_analysis/horscan_distance/{sample_s}/{sample_s}.{sample_t}.distance",
#     log:
#         "logs/horscan_distance/APG_v1_analysis/{sample_s}/{sample_t}.log",
#     conda:
#         "APG"
#     group: 
#         "distance"
#     script:
#         "../../scripts/03_hor_distance/pairwise_distance.py"


# ==============================================================================
# STAGE 1: 使用joblib并行生成所有.distance文件 (1个大任务)
# ==============================================================================
rule generate_all_distances:
    """
    运行一个超级脚本，在内部使用joblib并行生成所有.distance文件。
    这个规则的输出只是一个标志文件，表明此阶段完成。
    """
    input:
        blast_dir=ancient(config["blast_output_dir"]),
        format_flag=f"results/{config['run_name']}/flags/format_all.done"
    output:
        # ✅ 输出只是一个标志，而不是14万个文件
        done_flag=touch(f"results/{config['run_name']}/flags/all_distances_generated.done")
    params:
        graph_dir=f"results/{config['run_name']}/sample_data/graph/",
        sample_pair=SAMPLE_PAIR,
        hor_params=config['HORDistance'],
        output_root=f"results/{config['run_name']}/horscan_distance"
    threads:
        workflow.cores  
    log:
        f"logs/horscan_distance/{config['run_name']}/generate_all_distances.log"
    group: 
        "distance"
    conda:
        "APG"
    script:
        "../../scripts/03_hor_distance/generate_distances_in_parallel.py"


# ==============================================================================
# STAGE 2: 合并所有.distance文件并生成最终矩阵 (1个任务)
# ==============================================================================
rule merge_distance_matrix:
    """
    将所有由上一阶段生成的.distance文件合并成最终矩阵。
    """
    input:
        # ✅ 它依赖于上一阶段的完成标志
        flag=f"results/{config['run_name']}/flags/all_distances_generated.done"
    params:
        sample_order = sample_list,
        sample_pair= SAMPLE_PAIR,
        distance_dir=f"results/{config['run_name']}/horscan_distance"
    output:
        dist_dir=directory(f"results/{config['run_name']}/horscan_distance/distance"),
        done_flag=f"results/{config['run_name']}/flags/horscan_distance_matrix.done"
    log:
        f"logs/horscan_distance/{config['run_name']}/matrix.log",
    conda:
        "APG"
    script:
        "../../scripts/03_hor_distance/merge_matrix.py" 


rule plot_distance_heatmap:
    """
    为每个染色体和每个距离指标，生成一个聚类热图和对齐的HOR结构图。
    """
    input:
        # results/APG_v1_analysis/all_distance/distance/merge_distance_chr1.csv
        matrix=ancient("results/{run_name}/horscan_distance/distance/{metric}_distance_{chrom}.csv"),
    output:
        # ✅ 输出的PDF文件，文件名清晰地反映了其内容
        "plots/{run_name}/heatmap_{metric}_{chrom}.pdf"
    params:
        # ✅ 将config中的所有绘图参数作为一个整体传递给脚本
        qc_filter=config["sample_QC"],
        hor_summary=config["hor_summary"],
        sample_sheet=sample_list,
        plotting_params=config["plotting"]
    log:
        "logs/plotting/{run_name}/heatmap_{metric}_{chrom}.log"
    conda:
        "APG" 
    script:
        "../../scripts/plot_func/plot_heatmap_bar.py" # 指向我们适配后的脚本

rule plot_distance_heatmap_all:
    """
    生成所有染色体和所有指标的热图。
    """
    input:
        expand("plots/{run_name}/heatmap_{metric}_{chrom}.pdf",
               metric=['merge', 'aln', 'hor'],
               chrom=chromosomes,
               run_name=config["run_name"])
    output:
        touch(f"plots/{config['run_name']}/heatmap_all.done")
    log:
        f"logs/plotting/{config['run_name']}/heatmap_all.log"
    