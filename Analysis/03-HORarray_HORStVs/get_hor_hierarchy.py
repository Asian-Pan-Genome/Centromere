import pandas as pd
import numpy as np
import sys

class Node:
    def __init__(self, value, mnnum=None, hierarchy=None, horindex=None, bhlayer=None):
        self.value = value
        self.mnnum = mnnum
        self.hierarchy = hierarchy
        self.horindex = horindex
        self.bhlayer = bhlayer
        self.children = []
        self.parent = None

    def add_child(self, child_node):
        child_node.parent = self
        self.children.append(child_node)

    def get_parent(self):
        return self.parent

    def __repr__(self, level=0):
        ret = "\t" * level + repr(self.value) + repr(self.hierarchy) + f" (horindex: {self.horindex})" + f" (bhlayer: {self.bhlayer})""\n"
        for child in self.children:
            ret += child.__repr__(level + 1)
        return ret


def load_hor_stv():
    df = pd.read_csv("../../step-11/chrom_hor_phy/all.HiCAT.hor.summary.final.HOR.class.05.strandmerge.color.xls", sep="\t", header=0)
    horclass_color = dict(zip(df['HOR_class'], df['color']))
    hor_class = dict(zip(df['HOR'], df['HOR_class']))

    df = pd.read_csv("../../step-11/chrom_hor_phy/all.HiCAT.hor.summary.final.HOR.class.05.strandmerge.horcolor.xls", sep="\t", header=0)
    df['n-mer'] = df['HOR'].apply(lambda x : len(x.split('_')))
    df['index'] = range(1, len(df) + 1)
    df['HOR_class'] = df['HOR_class'].astype(str)
    df['index'] = df['index'].astype(str)
    df['n-mer'] = df['n-mer'].astype(str)
    df['horindex'] = "C" + df['HOR_class'] + "H" + df['index'] + "(" + df['n-mer'] + ")"
  
    hor_horindex = dict(zip(df['HOR'], df['horindex']))
    hor_horcolor = dict(zip(df['HOR'], df['hor_color']))
    
    return horclass_color, hor_class, hor_horindex, hor_horcolor


def build_hierarchy(idx, nestlines):
    '''
    idx marks the index of toplayer
    nestlines: 
    [
        [top_rowindex, topline],
        [cover_rowindex, coverline],
        ...
    ]
    '''
    idx_hor_start = 8
    idx_hor_start_idx = 10
    idx_hor_end_idx = 11
    idx_horstv = 24

    sorted_nestlines = sorted(nestlines, key=lambda x: (x[1][idx_hor_start], -(x[1][idx_hor_end_idx] - x[1][idx_hor_start_idx])))

    top_idx = idx
    top_rowidx = sorted_nestlines[0][0]
    top_line = sorted_nestlines[0][1]

    top_mnnum = top_line[idx_hor_end_idx] - top_line[idx_hor_start_idx] + 1
    top_horstv = top_line[idx_horstv]

    root = Node(top_rowidx, top_mnnum, top_idx, top_horstv, None)
    stack = [[root, top_line[idx_hor_start_idx], top_line[idx_hor_end_idx]]]

    for cover in sorted_nestlines[1:]:
        cover_rowidx = cover[0]
        cover_line = cover[1]
        cover_mnnum = cover_line[idx_hor_end_idx] - cover_line[idx_hor_start_idx] + 1
        cover_horstv = cover_line[idx_horstv]
        
        node = Node(cover_rowidx, cover_mnnum, None, cover_horstv, None)
        while stack and (stack[-1][-1]< cover_line[idx_hor_start_idx]):
            stack.pop()
        if stack:
            stack[-1][0].add_child(node)
        stack.append([node, cover_line[idx_hor_start_idx], cover_line[idx_hor_end_idx]])
    return sorted_nestlines, root

def get_hierarchy(root, i=1):
    '''
    for root, hierarchy is top_idx
    for cover, hierarchy is None and need be assigned 
    '''
    if root.hierarchy is not None:
        root_hierarchy = str(root.hierarchy)

    if root.children:
        for child in root.children:
            child.hierarchy = root_hierarchy + '.' + str(i)
            i += 1
            get_hierarchy(child, 1) 
    return root


def assign_bh_layer(root, threshold=0.6):

    hor_stv_ratio = {}
    total_mnnum = root.mnnum

    def get_ratio(root, total_mnnum, hor_stv_ratio):
        if root.children:
            for child in root.children:
                if child.horindex not in hor_stv_ratio:
                    hor_stv_ratio[child.horindex] = child.mnnum / total_mnnum
                else:
                    hor_stv_ratio[child.horindex] += child.mnnum / total_mnnum
                hor_stv_ratio = get_ratio(child, total_mnnum, hor_stv_ratio)
        return hor_stv_ratio

    hor_stv_ratio = get_ratio(root, total_mnnum, hor_stv_ratio)

    candidate_hor_stv = [horstv for horstv, ratio in hor_stv_ratio.items() if ratio > threshold]
    sorted_hor_stv_ratio = {k: v for k, v in sorted(hor_stv_ratio.items(), key=lambda item: item[1], reverse=True)}
    top_hor_stv = list(sorted_hor_stv_ratio.keys())[0]
    # print(len(candidate_hor_stv), top_hor_stv)

    def find_top_hor_stv(node, top_hor_stv):
        if node.horindex == top_hor_stv:
            return True
        for child in node.children:
            if find_top_hor_stv(child, top_hor_stv):
                return True
        return False

    def set_bhlayer(node, layer):
        # node.bhlayer = layer
        for child in node.children:
            child.bhlayer = layer
            set_bhlayer(child, layer)

    def find_bhlayer_with_top(node, top_hor_stv):
        for child in node.children:
            if find_top_hor_stv(child, top_hor_stv):
                if child.horindex == top_hor_stv:
                    child.bhlayer = 1
                    set_bhlayer(child, 0)
                else:
                    child.bhlayer = 0
                    find_bhlayer_with_top(child, top_hor_stv)
            else:
                child.bhlayer = 1
                set_bhlayer(child, 1)

    def get_bhlayer(root, candidate_hor_stv, top_hor_stv):
        if root.horindex == '-':
            # print(root.horindex)
            root.bhlayer = 0
            find_bhlayer_with_top(root, top_hor_stv)
        else:
            if len(candidate_hor_stv) == 0:
                root.bhlayer = 1
                set_bhlayer(root, 0)
            else:
                root.bhlayer = 0
                find_bhlayer_with_top(root, top_hor_stv)

    get_bhlayer(root, candidate_hor_stv, top_hor_stv)
    
    return root


def main(sample, hap, outfile):
    ##load hor stv id, class and color
    horclass_color, hor_class, hor_horindex, hor_horcolor = load_hor_stv()

    ##HiCAT##
    # sample = "CHM13"
    # hap = "-"
    summaryfile = "all.HiCAT.hor.summary.final.xls"
    # outfile = "CHM13.HiCAT.hor.summary.final.hierarchy.bhlayer.xls"

    ##match hor stv id and horclass with color##
    df = pd.read_csv(summaryfile, sep='\t', header=0)
    df['HOR_class'] = df['reorder_hor'].map(hor_class).fillna('-')
    df['color'] = df['HOR_class'].map(horclass_color).fillna('-')
    df['HORstv_index'] = df['reorder_hor'].map(hor_horindex).fillna('-')
    df['HORstv_color'] = df['reorder_hor'].map(hor_horcolor).fillna('-')
    

    def get_rowidx_newinfo(root, rowidx_hierarchy, rowidx_bhlayer):
        rowidx_hierarchy[root.value] = root.hierarchy
        rowidx_bhlayer[root.value] = root.bhlayer

        for child in root.children:
            rowidx_hierarchy, rowidx_bhlayer = get_rowidx_newinfo(child, rowidx_hierarchy, rowidx_bhlayer)

        return rowidx_hierarchy, rowidx_bhlayer
            

    ##single haploid assembly##
    hapdf = df[ (df['sample'] == sample ) & (df['hap'] == hap)].copy()

    updated_subdfs = []

    for ichr, subdf in hapdf.groupby('chromosome'):
        out = []
        subdf = subdf.reset_index(drop=True)
        ori_cols = subdf.columns.tolist()
        # print(ichr)
        # print(subdf.head())

        data = subdf.values.tolist()

        idx_layer = 14

        tophor = []
        '''
        tophor:
        [
            [[i, iline]], 
            [[j, jline], [k, kline]],  
        ]
        '''
        for i, line in enumerate(data):
            if line[idx_layer] == 'top':
                tophor.append([[i, line]])

        pivot = 0
        for i, line in enumerate(data):
            if line[idx_layer] == 'cover':
                if pivot + 1 < len(tophor):
                    while pivot + 1 < len(tophor) and i > tophor[pivot + 1][0][0]:
                        pivot += 1
                    if pivot + 1 < len(tophor):
                        if tophor[pivot][0][0] < i < tophor[pivot + 1][0][0]:
                            tophor[pivot].append([i, line])

        if tophor[-1][0][0] < len(data) - 1:     
            for i in range(tophor[-1][0][0] + 1, len(data)):
                if data[i][idx_layer] == 'cover':
                    tophor[-1].append([i, data[i]])
        

        for i, nestlines in enumerate(tophor):
            if len(nestlines) == 1:
                nestlines[0][1].append(i)
                nestlines[0][1].append(1)
                out.append(nestlines[0][1])
            else:
                sorted_nestlines, root = build_hierarchy(i, nestlines)
                get_hierarchy(root, 1) 
                assign_bh_layer(root, threshold=0.6)

                rowidx_hierarchy = {}
                rowidx_bhlayer = {}
                rowidx_hierarchy, rowidx_bhlayer = get_rowidx_newinfo(root, rowidx_hierarchy, rowidx_bhlayer)

                for rowidx, hierarchy in rowidx_hierarchy.items():
                    for line in sorted_nestlines:
                        if rowidx == line[0]:
                            line[1].append(hierarchy)
                            line[1].append(rowidx_bhlayer[rowidx])
                            out.append(line[1])
        updated_cols = ori_cols + ['hierarchy', 'bhlayer']
        updated_subdf = pd.DataFrame(out, columns=updated_cols)
        updated_subdfs.append(updated_subdf)

    update_hapdf = pd.concat(updated_subdfs, ignore_index=True)
    update_hapdf.to_csv(outfile, sep="\t", header=True, index=False)
    return update_hapdf
    
if __name__ == '__main__':
    sample = sys.argv[1]
    hap=sys.argv[2]
    outfile = f"HiCAT_bhlayer/{sample}_{hap}.HiCAT.hor.summary.final.hierarchy.bhlayer.xls"
    main(sample, hap, outfile)
