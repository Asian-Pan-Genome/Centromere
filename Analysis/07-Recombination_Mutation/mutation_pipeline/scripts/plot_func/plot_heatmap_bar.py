import os
import argparse
from typing import List, Dict, Tuple, Any

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib.collections import PatchCollection
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from tqdm import tqdm

# ==============================================================================
# 1. 全局配置 (Global Configuration)
# ==============================================================================

# 设置 Matplotlib 全局字体，以确保PDF/PS文件的字体嵌入正确，增强兼容性。
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

CONFIG = {
    "plot_params": {
        "row_height": 40,                   # 图中每个样本（行）的高度（单位：像素）
        "hor_plot_ratio": 1.5,              # HOR图相对于Heatmap的宽度比例
        "heatmap_cmap": sns.diverging_palette(10, 250, as_cmap=True), # 热图的颜色映射     
        "label_fontsize": 16,               # 坐标轴标签的字体大小
    },
    # --- 数据处理相关参数 ---
    "hor_params": {
        "clustering_method": "average",     # 层次聚类使用的算法 (e.g., 'average', 'ward', 'complete')
        "hor_layer_to_use": "top",          # 从HOR注释文件中选择要使用的层级
        "missing_hor_color": "#adadad",   # 当HOR数据缺失时使用的颜色
    }
}
# ==============================================================================
# 2. 数据加载与处理函数 (Data Loading & Processing Functions)
# ==============================================================================

def load_and_filter_matrix(
    matrix_path: str, sample_list: List[str], filter_path: str, chrom: str
) -> Tuple[np.ndarray, List[str]]:
    """
    加载距离CSV，根据QC文件进行筛选，并对矩阵进行预处理。

    Args:
        matrix_path (str): .csv格式的距离矩阵文件路径，第一行是表头，第一列是索引。
        sample_list (List[str]): 预期的、完整的样本名称列表，用于和QC表对齐。
        filter_path (str): 包含筛选信息的Excel文件路径。
        chrom (str): 当前处理的染色体名称。

    Returns:
        Tuple[np.ndarray, List[str]]: 一个元组，包含筛选和处理后的距离矩阵(NumPy数组)以及对应的样本列表。
    """
    print(f"Loading distance matrix from: {matrix_path}")
    if not os.path.exists(matrix_path):
        raise FileNotFoundError(f"Matrix file not found: {matrix_path}")
        
    # header=0表示第一行为列名, index_col=0表示第一列为行索引
    distance_matrix = pd.read_csv(matrix_path, header=0, index_col=0)

    # --- ✅ 优化点 1: 使用矢量化操作替代for循环，大幅提升性能 ---
    print(f"Loading and processing QC filter from: {filter_path}")
    filter_df = pd.read_excel(filter_path)
    chr_filter_df = filter_df[filter_df['chrom'] == chrom]
    
    # 1. 直接从QC表中筛选出所有通过QC的样本 (filterflag == 0)
    qc_passed_samples = set(chr_filter_df[chr_filter_df['filterflag'] == 0]['sample_hap'])
    
    # 2. 从传入的完整样本列表中，筛选出那些在QC通过列表中的样本
    #    这确保了样本的顺序与原始 sample_list 一致
    samples_to_keep = [s for s in sample_list if s in qc_passed_samples]
    
    # 3. 另外，确保这些样本也存在于距离矩阵的索引中，防止KeyError
    samples_in_matrix = set(distance_matrix.index)
    final_filtered_samples = [s for s in samples_to_keep if s in samples_in_matrix]

    if not final_filtered_samples:
        print(f"Warning: No samples passed QC for chromosome {chrom}. Returning empty matrix.")
        return np.array([]), []

    print(f"Filtered samples for {chrom}. Original: {len(sample_list)}, Passed QC: {len(samples_to_keep)}, Final in matrix: {len(final_filtered_samples)}")
    
    # --- ✅ 核心错误修复 ---
    # 1. 使用 .loc 按标签（行名和列名）进行索引
    # 2. 使用 .values 将筛选后的DataFrame转换为NumPy数组，以匹配函数的返回类型提示
    filtered_matrix_df = distance_matrix.loc[final_filtered_samples, final_filtered_samples]
    filtered_matrix_np = filtered_matrix_df.values
    
    # ✅ 注释修正：原注释中提到的 np.ix_ 是用于NumPy数组的，对于Pandas DataFrame，.loc是正确的选择。
    
    return filtered_matrix_np, final_filtered_samples

def load_and_process_hor_data(
    hor_path: str, sample_list: List[str], params: Dict[str, Any]
) -> List[pd.DataFrame]:
    """
    加载并处理HOR summary文件，为每个样本生成一个包含位置和颜色信息的DataFrame。
    
    Args:
        hor_path (str): HOR summary文件的路径。
        sample_list (List[str]): 经过筛选的样本列表。
        params (Dict[str, Any]): 包含处理所需参数的字典。

    Returns:
        List[pd.DataFrame]: 一个列表，每个元素是对应样本的HOR数据DataFrame。
    """
    print(f"Loading HOR data from: {hor_path}")
    hor_df = pd.read_csv(hor_path, header=None, sep='\t')
    # 为长HOR文件分配列名
    hor_df.columns = [
        "sample", "hap", "blockindex", "block_start", "block_end", "blocklen", 
        "mnnum", "chrom", "horstart", "horend", "hor_index_start", "hor_index_end", 
        "nrepeat", "hor", "layer", "reorder_hor", "sample_hap_chrom", "sample_hap", 
        "chromosome", "nmer", "hor_flag", "pan_hor_type", "HOR_class", "color", 
        "HORstv_index", "HORstv_color", "hierarchy", "bhlayer"
    ]
    
    # 筛选指定层级的数据，并创建副本以避免后续的SettingWithCopyWarning
    hor_df_filtered = hor_df.copy()
    
    # 填充颜色列中的缺失值，使用赋值操作替换inplace=True，以避免FutureWarning
    hor_df_filtered['HORstv_color'] = hor_df_filtered['HORstv_color'].fillna(params["missing_hor_color"])
    hor_df_filtered['HORstv_color'] = hor_df_filtered['HORstv_color'].replace('-', params["missing_hor_color"])

    processed_data = []
    for sample in tqdm(sample_list, desc="Processing HOR data for samples"):
        # 对CHM13样本进行特殊处理
        query = 'CHM13' if 'CHM13' in sample else sample
        df_sample = hor_df_filtered[hor_df_filtered["sample_hap"].str.contains(query, na=False)]
        
        if not df_sample.empty:
            df = df_sample[["horstart", "horend", "HORstv_color"]].copy()
            df.columns = ["start", "end", "color"]
            # 对坐标进行归一化，使每个样本的HOR图都从0开始
            offset = df["start"].min()
            df["start"] -= offset
            df["end"] -= offset
            processed_data.append(df)
        else:
            print(f"Missing HOR data for: {sample}")
            # 如果样本没有HOR数据，添加一个空的DataFrame以保持列表长度一致
            processed_data.append(pd.DataFrame(columns=["start", "end", "color"]))
            
    return processed_data

def cluster_and_reorder_data(
    matrix: np.ndarray, labels: List[str], hor_data: List[pd.DataFrame], method: str
) -> Tuple[np.ndarray, List[str], List[pd.DataFrame]]:
    """
    对数据进行层次聚类，并根据聚类结果重新排序矩阵、标签和HOR数据。
    """
    print("Performing hierarchical clustering...")
    # 使用pdist计算成对距离，然后用linkage进行聚类
    linkage_matrix = linkage(pdist(matrix), method=method)
    # 获取叶节点的顺序，并反转以符合视觉习惯（通常相似的聚在一起）
    ordered_indices = leaves_list(linkage_matrix)[::-1]

    # 根据得到的索引重新排序所有相关数据
    reordered_matrix = matrix[ordered_indices][:, ordered_indices]
    reordered_labels = [labels[i] for i in ordered_indices]
    reordered_hor_data = [hor_data[i] for i in ordered_indices]
    
    return reordered_matrix, reordered_labels, reordered_hor_data

# ==============================================================================
# 3. 绘图核心函数 (Core Plotting Functions)
# ==============================================================================

def draw_heatmap(ax: plt.Axes, matrix: np.ndarray, labels: List[str], params: Dict[str, Any]):
    """
    使用matplotlib的imshow函数高效地绘制热图。
    """
    # 截取1%到99%分位点的值作为颜色映射的范围，以减少极端值的影响
    vmin = np.percentile(matrix, 1)
    vmax = np.percentile(matrix, 99)
    
    # imshow是绘制矩阵/图像的最高效方法
    ax.imshow(matrix, cmap=params["heatmap_cmap"], vmin=vmin, vmax=vmax, interpolation='none', aspect='equal')
    
    # 设置刻度和标签
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=params["label_fontsize"])
    ax.set_yticklabels(labels, rotation=0, fontsize=params["label_fontsize"])
    
    # 将X轴刻度移动到顶部，这是热图的常见样式
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

def draw_hor_plot(ax: plt.Axes, hor_data: List[pd.DataFrame], labels: List[str], params: Dict[str, Any]):
    """
    绘制与热图对齐的HOR条带图。
    """
    row_height = params["row_height"]
    all_patches = []
    max_x = 0
    
    for i, df in enumerate(tqdm(hor_data, desc="Drawing HOR plot")):
        if df.empty: continue
            
        y_pos = i * row_height
        max_x = max(df["end"].max(), max_x) if not df.empty else max_x
        
        # 为每个HOR片段创建一个矩形补丁
        for _, row in df.iterrows():
            rect = patches.Rectangle(
                (row["start"], y_pos), 
                row["end"] - row["start"], 
                row_height * 0.8, # 使条带之间有空隙
                facecolor=row["color"],
                linewidth=0
            )
            all_patches.append(rect)

    # 一次性将所有补丁添加到坐标轴上，比逐个添加效率高
    ax.add_collection(PatchCollection(all_patches, match_original=True))
    
    # 设置坐标轴范围和刻度
    ax.set_xlim(0, max_x)
    ax.set_ylim(0, len(labels) * row_height)
    
    # 反转Y轴，使其0点在顶部，与热图的imshow行为保持一致
    ax.invert_yaxis()
    
    # 设置Y轴刻度位置，但隐藏标签，因为标签已在热图上显示
    ax.set_yticks(np.arange(len(labels)) * row_height + row_height / 2)
    ax.set_yticklabels([])
    ax.tick_params(axis='y', length=0) # 隐藏Y轴的小刻度线

    # 隐藏不需要的边框线，使图形更简洁
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)

def plot_main_figure(
    matrix: np.ndarray, labels: List[str], hor_data: List[pd.DataFrame], 
    params: Dict[str, Any], output_file: str
):
    """
    创建并组合热图和HOR图，并保存为文件。
    """
    num_labels = len(labels)
    row_height = params["row_height"]
    
    # 动态计算图形的尺寸（单位：英寸）
    heatmap_size_inch = num_labels * row_height / 100 
    hor_width_inch = heatmap_size_inch * params["hor_plot_ratio"]
    figsize = (heatmap_size_inch + hor_width_inch, heatmap_size_inch)

    # 创建子图，并使用gridspec_kw精确控制两个子图的宽度比例
    fig, (ax_heatmap, ax_hor) = plt.subplots(
        ncols=2, 
        figsize=figsize, 
        gridspec_kw={"width_ratios": [heatmap_size_inch, hor_width_inch]},
        constrained_layout=True # 自动调整布局防止重叠
    )
    
    print("Drawing heatmap...")
    draw_heatmap(ax_heatmap, matrix, labels, params)
    
    print("Drawing HOR plot...")
    draw_hor_plot(ax_hor, hor_data, labels, params)
    
    print(f"Saving figure to {output_file}...")
    plt.savefig(output_file, format="pdf", bbox_inches="tight")
    plt.close(fig)
    

# ==============================================================================
# 4. 主执行逻辑 (Main Execution Logic) - ✅ 新增部分
# ==============================================================================
def main():
    """
    脚本的主执行函数。
    它会解析Snakemake对象，加载数据，进行聚类，然后绘图并保存。
    """
    # --- a. 从Snakemake对象获取输入、输出和参数 ---
    # 这种方式比命令行参数解析更简洁，是Snakemake `script:`指令的推荐用法
    try:
        # 输入文件
        matrix_path = snakemake.input.matrix
        full_sample_list = snakemake.params.sample_sheet
        filter_path = snakemake.params.qc_filter
        hor_path = snakemake.params.hor_summary
        
        # 输出文件
        output_file = snakemake.output[0]
        
        # 通配符
        chrom = snakemake.wildcards.chrom
        
        # 日志文件
        log_file = snakemake.log[0]

    except NameError:
        sys.exit("This script is designed to be run via the Snakemake `script:` directive.")

    # --- b. 设置日志重定向 ---
    sys.stdout = open(log_file, 'w')
    sys.stderr = sys.stdout
    
    print(f"--- Starting Heatmap Generation for Chromosome: {chrom} ---")

    # --- c. 执行核心工作流 ---
    # 1. 加载完整的样本列表

    # 2. 加载并筛选距离矩阵
    filtered_matrix, filtered_samples = load_and_filter_matrix(
        matrix_path, full_sample_list, filter_path, chrom
    )
    
    # 3. 加载并处理HOR数据
    hor_data = load_and_process_hor_data(hor_path, filtered_samples, CONFIG['hor_params'])
    
    # 4. 对数据进行聚类和重排序
    reordered_matrix, reordered_labels, reordered_hor_data = cluster_and_reorder_data(
        filtered_matrix, filtered_samples, hor_data, CONFIG['hor_params']["clustering_method"]
    )
    
    # 5. 绘制主图
    plot_main_figure(
        reordered_matrix, reordered_labels, reordered_hor_data, CONFIG['plot_params'], output_file
    )
    
    print(f"--- Figure generation for {chrom} complete. ---")


# ===================================================================
# 5. 脚本执行入口 - ✅ 新增部分
# ===================================================================
if __name__ == "__main__":
    # 当此脚本由Snakemake的`script:`指令调用时，会直接执行main()函数。
    # Snakemake会在执行前自动将`snakemake`对象注入到全局命名空间中。
    main()