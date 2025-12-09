import numpy as np
from collections import defaultdict,Counter
import os
import math
# import networkx as nx
# from networkx.drawing.nx_agraph import write_dot
from subprocess import check_call
from get_continuous_block import regex_mn
import csv


def SingleBlockClusterStat(iblock, seqid_clustid):
    mns = iblock[5].split(',')
    cids = []
    for mn in mns:
        mn = mn.strip()
        cid = seqid_clustid.get(mn,'-')
        cids.append(cid)
    count = Counter(cids)
    return cids, count

def invert_type(mns):
    strands = []
    mns = mns.split(',')
    for mnid in mns:
        ichr, start, end, strand = regex_mn(mnid)
        if strand == "+":
            strands.append(1)
        else:
            strands.append(0)

    p_plus = strands.count(1) / len(strands)
    p_minus = strands.count(0) / len(strands)

    # print(",".join(map(str,strands)))
    def local_invert_region(strands, target):
        max_count = 0
        current_count = 0
        for i in strands:
            if i == target:
                current_count += 1
                max_count = max(max_count, current_count)
            else:
                current_count = 0
        target_total = strands.count(target)
        # print(target, max_count, target_total, len(strands))
        if max_count / target_total > 0.2:
            return True
        else:
            return False

    ty = ""
    if p_plus > 0.8:
        ty = "+"
    elif p_minus > 0.8:
        ty = "-"
    elif p_plus > 0.5 and local_invert_region(strands, 0):
        ty = "+"
    elif p_minus > 0.5 and local_invert_region(strands, 1):
        ty = "-"
    else:
        ty = "mix"
    return strands, ty


def get_nodes(cids):
    nodes = defaultdict(int)
    for i in cids:
        i = (i,)
        nodes[i] += 1
    return nodes


def get_edges(cids, strands, ty):
    edges = defaultdict(int)
    if ty != "mix":
        for i in range(len(cids) - 1):
            strand = strands[i + 1]
            if strand == 1:
                edges[(cids[i], cids[i + 1])] += 1
            else:
                #edges[(cids[i + 1], cids[i])] += 1
                edges[(cids[i], cids[i + 1])] += 1
    else:
        for i in range(len(cids) - 1):
            edges[(cids[i], cids[i + 1])] += 1
    return edges

def get_edge_thr(edges, top_n=3, level="all"):
    bins = [1000, 500, 200,  100, 75, 50, 10, 1, 0]
    def count_ranges(data, bins):
        counts = Counter()
        for value in data:
            for bin_range in bins:
                if value > bin_range:
                    counts[bin_range] += 1
                    break
        return counts

    def find_min_ranges(counts, top_n):
        filtered_counts = {k: v for k, v in counts.items() if k > 10}
        if len(filtered_counts) > 3:
            top_ranges = Counter(filtered_counts).most_common(top_n)
            min_range = min(top_ranges, key=lambda x: x[0])
        elif len(filtered_counts) == 0:
            min_range = list(counts.items())[0]
        else:
            min_range = min(list(filtered_counts.items()), key=lambda x:x[0])
        return min_range

    ec = list(edges.values())
    range_counts = count_ranges(ec, bins)
    if level == "all":
        min_range = find_min_ranges(range_counts, top_n)
        flag_value = min_range
    else:
        flag_value = max(list(range_counts.items()), key=lambda x: x[0])
    print("flag_value:", flag_value)
    if flag_value[0] == 1000:
        thr=100
    elif flag_value[0] == 500:
        thr = 50
    elif flag_value[0] == 200:
        thr = 20
    elif flag_value[0] == 100:
        #thr = 10
        #thr = 15
        thr = 20
    elif flag_value[0] == 75:
        thr = 10
    elif flag_value[0] == 50:
        thr = 9
    elif flag_value[0] == 10:
        thr = 3
    else:
        thr = 2
    return thr

def addNode(G, mn, mncnt):
    lg = 0.01
    if mncnt > 0:
        lg = math.log(mncnt)
    clr = ["red", "#cd5700", "orange", "#ffb300", "yellow", "#ceff1d", "#a7fc00", "#00ff00", "#e0ffff", "#f5fffa"]
    G.add_node(mn, style=f'filled', fillcolor=f'{clr[int(lg)]}', label=f'{mn}[{str(int(mncnt))}]')


def addEdges(G, mn1, mn2, cnt, thr):
    scr = cnt
    thr_wg = [100000000, 1000, 500, 100, 1]
    wgs = [7, 5, 3, 1, 0]
    wg = 3
    while (scr > thr_wg[wg]):
        wg -= 1

    if scr > thr:
        G.add_edge(mn1, mn2, label=f'{scr}', penwidth=f'{wgs[wg]}')


def BuildMonomerGraph(nodes, edges, iblock, nodeThr=2, edgeThr=2):
    mn_set = {tuple(list(x)[:-1]) for x, y in edges.items() if y > nodeThr} | \
             {tuple(list(x)[1:]) for x, y in edges.items() if y > nodeThr}

    graph_label = "#%s %s:%s-%s" % (iblock[0], iblock[1], iblock[2], iblock[3])
    G = nx.DiGraph(label=graph_label)

    for vert in mn_set:
        cnt_mn = 0

        if vert in nodes:
            cnt_mn = nodes[vert]
            lg = 0
            if cnt_mn > 0:
                lg = math.log(cnt_mn)

        clr = ["red", "#cd5700", "orange", "#ffb300", "yellow", "#ceff1d", "#a7fc00", "#00ff00", "#e0ffff", "#f5fffa",
               "#f5fffa", "#f5fffa", "#f5fffa", "#f5fffa"]
        vert = "-".join(list(vert))
        curc = clr[int(lg)]

        G.add_node(vert, style="filled", fillcolor=curc, label=r"%s\n[%s]" % (vert, str(int(cnt_mn))))

    ecnt = 0

    for vt1 in sorted(mn_set):
        for vt2 in sorted(mn_set):
            if list(vt1)[1:] != list(vt2)[:-1]:
                continue
            scr = 0
            if (*vt1, vt2[-1]) in edges:
                scr = edges[(*vt1, vt2[-1])]
            thr_wg = [100000000, 1000, 500, 100, 1]
            wgs = [7, 5, 3, 1, 0]
            wg = 3
            while (scr > thr_wg[wg] and wg > 0):
                wg -= 1

            if scr > edgeThr:
                ecnt += 1
                vrt1 = "-".join(list(vt1))
                vrt2 = "-".join(list(vt2))
                if scr < 1:
                    G.add_edge(vrt1, vrt2, label=str(int(scr)), penwidth=str(wgs[wg]), constraint="false")
                else:
                    G.add_edge(vrt1, vrt2, label=str(int(scr)), penwidth=str(wgs[wg]))
    return G


def drawGraph(G, outprefix, outdir):
    write_dot(G, os.path.join(outdir, outprefix + ".dot"))
    check_call(
        ['dot', '-Tpng', os.path.join(outdir, outprefix + ".dot"), '-o', os.path.join(outdir, outprefix + ".png")])

def splitGraph(G):
    weakly_connected = list(nx.weakly_connected_components(G))
    subgraphs = [G.subgraph(component) for component in weakly_connected]
    return subgraphs

def genCycleInner(G, prefixCycle, usedEdges, cycleList):
    def samecc(c1, c2):
        if len(c1) != len(c2):
            return False
        for j in range(len(c1) - 1):
            if c1[j:-1] + c1[:j] == c2[:-1]:
                return True
        return False

    if len(prefixCycle) > 1 and prefixCycle[0] == prefixCycle[-1]:
        for cc in cycleList:
            if samecc(cc, prefixCycle):
                break
        else:
            cycleList.append(prefixCycle)

    v = prefixCycle[-1]
    for e in G.edges(v):
        if e not in usedEdges:
            genCycleInner(G, prefixCycle + [e[1]], usedEdges | {e}, cycleList)


def genAllCycles(G):
    cycleList = []
    for v in G.nodes():
        genCycleInner(G, [v], set(), cycleList)
    return cycleList

def getCycleCnt(cl, mncen):
    clcnt = 0
    for i in range(len(mncen) - len(cl)):
        if mncen[i:i+len(cl)] == cl:
            clcnt += 1

    ccl = [x + "'" for x in cl[len(cl)::-1]]

    for i in range(len(mncen) - len(cl)):
        if mncen[i:i+len(cl)] == ccl:
            clcnt += 1

    return clcnt

def filterCycles(cycles, hybridSet, mncen, minTrav):
    usedV = set()
    cycles.sort(key=lambda x: -len(x))
    res_cyc = []
    for cl in cycles:
        if len([v for v in cl if v in hybridSet]) > 0:
            continue

        print("cl:", cl, "clcnt:", getCycleCnt(cl, mncen))
        if getCycleCnt(cl, mncen) < minTrav:
            continue

        for v in cl:
            if v not in usedV:
                res_cyc.append(cl)
                break
        usedV |= set(cl)
    return res_cyc

def saveHOR(cycles, outprefix, outdir):
    outfile = os.path.join(outdir, outprefix + ".HORs.tsv")
    with open(outfile, "w") as fw:
        csv_writer = csv.writer(fw, delimiter='\t')
        horid = 1
        for cycle in cycles:
            csv_writer.writerow(["H" + str(horid), ",".join(cycle)])
            horid += 1
    return outfile

def get_iblock_mn_pos(iblock):
    mns = iblock[5].split(',')
    starts = []
    ends = []
    outchr = ""
    for mnid in mns:
        ichr, start, end, strand = regex_mn(mnid)
        starts.append(int(start))
        ends.append(int(end))
        outchr = ichr
    return outchr, starts, ends
 
def decompose_sequence(seq, cycles, outprefix, outdir, outchr, starts, ends):
    outfile = os.path.join(outdir, outprefix + ".HORdecomposition.tsv")
    cycles = [cl[:-1] for cl in cycles]

    def list_in(c, cycle_rotate):
        count_c = len(c)
        for i in cycle_rotate:
            count_i = len(i)
            if count_c == count_i:
                rst = np.zeros(count_c, dtype=bool)
                for j in range(count_c):
                    if c[j] == i[j]:
                        rst[j] = True
                    else:
                        break
                if np.all(rst):
                    return True
        return False

    def sub_in(sub, full):
        count_s = len(sub)
        count_f = len(full)
        if count_f >= count_s:
            for i in range(count_f - count_s + 1):
                rst = np.zeros(count_s, dtype=bool)
                for j in range(count_s):
                    if sub[j] == full[i + j]:
                        rst[j] = True
                    else:
                        break
                if np.all(rst):
                    return True
        return False

    def list_big_than(lhs, rhs):
        count_l = len(lhs)
        count_r = len(rhs)
        if count_l > count_r:
            rst = np.zeros(count_r, dtype=bool)
            for j in range(count_r):
                if lhs[j] == rhs[j]:
                    rst[j] = True
                else:
                    break
            if np.all(rst):
                return True
        return False

    with open(outfile, "w", newline="") as fw:
        csv_writer = csv.writer(fw, delimiter='\t')
        decomposed = []
        index = 0
        while index < len(seq):
            found_subtring = ""
            out_substring = ""
            terminate_flag = False
            found_substring_list = []
            cycle_substring_list = []
            single_base = []
            for cycle in cycles:
                found_substring = ""
                cycle_rotate = [cycle[i:] + cycle[:i] for i in range(len(cycle))]
                if list_in(seq[index:index+len(cycle)], cycle_rotate):
                    decomposed.append(seq[index:index+len(cycle)])
                    csv_writer.writerow([outprefix, index, index+len(cycle), outchr, starts[index], ends[index+len(cycle)-1], ",".join(seq[index:index+len(cycle)])])
                    index += len(cycle)
                    break
                else:
                    i = 1
                    while True:
                        substring = seq[index:index + i]
                        if index + i == len(seq):
                            terminate_flag = True
                            break
                        if sub_in(substring, cycle*2):
                            out_substring = substring
                            i += 1
                        elif len(out_substring) > 1:
                            found_substring_flag = True
                            found_substring_list.append(found_substring_flag)
                            cycle_substring_list.append(out_substring)
                            out_substring=""
                            break
                        else:
                            if len(substring) == 1:
                                single_base.append(substring)
                                break
                            else:
                                single_base.append(out_substring)
                                break
                    if found_substring_list != []:
                        for k in range(len(found_substring_list)):
                            if found_substring_list[k]:
                                if  list_big_than(cycle_substring_list[k], found_subtring):
                                    found_subtring = cycle_substring_list[k]
                        decomposed.append(found_subtring)
                        csv_writer.writerow([outprefix, index, index + len(found_subtring), outchr, starts[index], ends[index + len(found_subtring) -1], ",".join(found_subtring)])
                        index += len(found_subtring)
                        break
                    elif terminate_flag:
                        decomposed.append(substring)
                        csv_writer.writerow([outprefix, index, index + len(substring), outchr, starts[index], ends[index + len(substring) -1], ",".join(substring)])
                        index += len(substring)
                        break
                    else:
                        if len(single_base) == len(cycles):
                            decomposed.append(single_base[0])
                            csv_writer.writerow([outprefix, index, index + len(single_base[0]), outchr, starts[index], ends[index + len(single_base[0]) -1], ",".join(single_base[0])])
                            index += len(single_base[0])
                            break
    return outfile

    
def main_graph(cblocks, seqid_clustid, outdir):
    hormon_outdir = os.path.join(outdir, "graph")
    os.makedirs(hormon_outdir, exist_ok=True)

    #index_blocks = {i: iblock for i, iblock in enumerate(cblocks)}
    #sorted_blocks = dict(sorted(index_blocks.items(), key=lambda item: item[1][6]))

    for i, iblock in enumerate(cblocks):
        #if i != 0:
        #    continue
        cids, count = SingleBlockClusterStat(iblock, seqid_clustid)
        with open(os.path.join(hormon_outdir, str(i) + ".chain.xls"), 'w') as outf1:
            outf1.write("#" + ",".join(map(str, iblock[:5])) + "," + str(iblock[-1]) + "\n")
            outf1.write(",".join(cids) + "\n")

        with open(os.path.join(hormon_outdir, str(i) + ".chain.stat.xls"), 'w') as outf2:
            outf2.write("#" + ",".join(map(str, iblock[:5])) + "," + str(iblock[-1]) + "\n")
            sort_count = dict(sorted(count.items(), key=lambda item: item[1], reverse=True))
            for mnid, c in sort_count.items():
                outf2.write(mnid + "\t" + str(c) + "\n")

        mns = iblock[5]
        strands, ty = invert_type(mns)
        nodes = get_nodes(cids)
        edges = get_edges(cids, strands, ty="none")
        #if ty == "-":
        #    cids = cids[::-1]
        #    print(cids)
        
        if len(edges) > 1:
            edthr = get_edge_thr(edges, 3, "all")
            G = BuildMonomerGraph(nodes, edges, iblock, nodeThr=2, edgeThr=edthr)
            drawGraph(G, str(i), hormon_outdir)

            subgraphs =  splitGraph(G)
            tmp = []
            for j, g in enumerate(subgraphs):
                num_nodes = g.number_of_nodes()
                g_edges = g.edges()
                if num_nodes > 1 and len(g_edges) > 0:
                    drawGraph(g, str(i) + "." + str(j), hormon_outdir)
                    edge_weights = {k: edges[k]for k in edges if k in g_edges}
                    iedthr = get_edge_thr(edge_weights, 1, "sub")
                    g_copy = g.copy(as_view=False)
                    for (u,v),w  in edge_weights.items():
                        if w < iedthr:
                            g_copy.remove_edge(u,v)
                    drawGraph(g_copy, str(i) + "." + str(j) + ".updated", hormon_outdir)
                    cycles = genAllCycles(g_copy)
                    print(str(i) + "." + str(j), cycles)
                    cycles = filterCycles(cycles, {}, cids, 2)
                    print(str(i) + "." + str(j), "filtered", cycles)
                    tmp += cycles
                    print(tmp)

            final_cycles_set = set(tuple(item) for item in tmp)
            final_cycles = [list(item) for item in final_cycles_set]
            print('-------------------------')
            print(final_cycles)
            saveHOR(final_cycles, str(i), hormon_outdir)
            if len(final_cycles) > 0 :
                outchr, starts, ends = get_iblock_mn_pos(iblock)
                decompose_sequence(cids, final_cycles, str(i), hormon_outdir, outchr, starts, ends)
            print("#" + ",".join(map(str, iblock[:5])) + "," + str(iblock[-1]) + " finished.")
        else:
            print("#" + ",".join(map(str, iblock[:5])) + "," + str(iblock[-1]) + " has single monomer. done!")

