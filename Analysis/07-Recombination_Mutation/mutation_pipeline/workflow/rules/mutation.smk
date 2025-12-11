# 读取blast结果，计算blast的分数，结合QC文件，筛选出得分超过0.9以上的样本对
rule filter_synteny_pair:
    input:
        blast_dir=config["blast_output_dir"],
    output:
        filtered_sample=protected(f"results/{config['run_name']}/mutation/filtered_synteny_sample.bed"),
    params:
        qc_file=config["sample_QC"],
        sample_list=sample_list
    log:
        f"logs/{config['run_name']}/mutation/filtered_synteny_sample.log",
    threads:
        workflow.cores
    conda:
        "APG"
    script:
        "../../scripts/04_mutation/filter_synteny_sample.py"


rule filter_synteny_threshold:
    input:
        filtered_sample=ancient(f"results/{config['run_name']}/mutation/filtered_synteny_sample.bed"),
    output:
        filtered_sample=f"results/{config['run_name']}/mutation/filtered_synteny_threshold.bed",
    log:
        f"logs/{config['run_name']}/mutation/filtered_synteny_threshold.log",
    params:
        threshold=config["synteny_threshold"],
    conda:
        "APG"
    script:
        "../../scripts/04_mutation/filter_synteny_threshold.py"


rule monomer_mutation_calculate:
    input:
        filtered_sample=f"results/{config['run_name']}/mutation/filtered_synteny_threshold.bed",
        blast_dir=config["blast_output_dir"],
    output:
        monomer_snp_dir=directory(
            f"results/{config['run_name']}/mutation/monomer_mutation_calculate"
        ),
    params:
        sample_bed=config["sample_sheet"],
    log:
        f"logs/{config['run_name']}/mutation/monomer_mutation_calculate.log",
    threads: workflow.cores
    conda:
        "APG"
    script:
        "../../scripts/04_mutation/monomer_mutation_calculate.py"


rule variant_monomer_mask:
    input:
        filtered_sample=f"results/{config['run_name']}/mutation/filtered_synteny_threshold.bed",
        monomer_snp_dir=f"results/{config['run_name']}/mutation/monomer_mutation_calculate",
    output:
        monomer_mask_dir=directory(
            f"results/{config['run_name']}/mutation/variant_monomer_mask"
        ),
    params:
        mask_params=config["mask_params"],
    threads: workflow.cores
    log:
        f"logs/{config['run_name']}/mutation/variant_monomer_mask.log",
    conda:
        "APG"
    script:
        "../../scripts/04_mutation/variant_monomer_mask.py"


rule analyze_mask_stats:
    input:
        filtered_sample=f"results/{config['run_name']}/mutation/filtered_synteny_threshold.bed",
        monomer_mask_dir=f"results/{config['run_name']}/mutation/variant_monomer_mask",
    output:
        f"results/{config['run_name']}/mutation/filtered_synteny_with_mask_stats.bed",
    params:
        mask_column="mask",
    log:
        f"logs/{config['run_name']}/mutation/analyze_mask_stats.log",
    threads: workflow.cores  # 脚本将使用此数量的线程进行并行处理
    conda:
        "APG"  # 与您的环境保持一致
    script:
        "../../scripts/04_mutation/variant_mask_filter.py"


rule filter_mask_stats:
    input:
        f"results/{config['run_name']}/mutation/filtered_synteny_with_mask_stats.bed",
    output:
        f"results/{config['run_name']}/mutation/filter_mask_stats.bed",
    params:
        mask_filter=config["mask_filter"],
    log:
        f"logs/{config['run_name']}/mutation/filter_mask_stats.log",
    conda:
        "APG"
    script:
        "../../scripts/04_mutation/filter_mask_stats.py"


rule plot_synteny_sample_pair:
    input:
        filtered_sample=f"results/{config['run_name']}/mutation/filtered_synteny_threshold.bed",
    output:
        plot_dir=directory(
            f"results/{config['run_name']}/mutation/synteny_sample_pair_plot"
        ),
    log:
        f"logs/{config['run_name']}/mutation/synteny_sample_pair_plot.log",
    conda:
        "APG"
    script:
        "../../scripts/04_mutation/plot_synteny_sample_all.py"


rule phylogenetic_monomer_validation:
    input:
        pairwise_mask_dir=f"results/{config['run_name']}/mutation/sv_high_variant_monomer_mask",
        qc_file=config["sample_QC"],
    output:
        phylogenetic_validation_dir=directory(
            f"results/{config['run_name']}/mutation/phylogenetic_monomer_validation"
        ),
    log:
        f"logs/{config['run_name']}/mutation/phylogenetic_monomer_validation.log",
    conda:
        "APG"
    script:
        "../../scripts/04_mutation/phylogenetic_monomer_validation.py"


rule plot_synteny_sample_pair_cent:
    input:
        filter_sample=f"results/{config['run_name']}/mutation/filter_mask_stats.bed",
        mask_dir=f"results/{config['run_name']}/mutation/variant_monomer_mask",
    params:
        sample_bed=config["sample_sheet"],
    output:
        plot_dir=directory(
            f"results/{config['run_name']}/mutation/synteny_sample_pair_cent_plot"
        ),
    log:
        f"logs/{config['run_name']}/mutation/synteny_sample_pair_cent_plot",
    threads:
        workflow.cores
    conda:
        "APG"
    script:
        "../../scripts/04_mutation/synteny_sample_pair_cent_plot.py"
