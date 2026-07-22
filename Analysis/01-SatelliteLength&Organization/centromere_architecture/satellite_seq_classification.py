# -*- coding: utf-8 -*-
"""
Third classification scheme for the ACROCENTRIC chromosomes
(chr13, chr14, chr15, chr21, chr22):

  * array length threshold = 10kb (fixed here, independent of the main script).
  * HSat1A/HSat1B collapsed to HSat1 (canon_sat, as in the main pipeline).
  * CLASS = the ORDERED sequence of array types along the centromere coordinate,
    e.g.  ASat>BSat>HSat1>HSat3 .
    (Unlike method1 which only counts arrays and ignores order, and method2 which
     counts adjacent transitions.)  Consecutive same-type blocks within <=500kb are
     already merged into one array by segment_arrays; the remaining array order is
     kept verbatim (no second run-length collapse).

Outputs (results/):
  satellite_seq_classification_acro.tsv         per-centromere ordered-sequence class
                                                (rows in the SAME order as the PDFs)
  figures/<chr>_seq_abs.pdf / <chr>_seq_norm.pdf
                                                track plot only (NO heatmap), rows in
                                                class order matching the TSV.

Run:  python satellite_seq_classification.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.collections import PatchCollection

import satellite_classification as sc  # reuse loading / segmentation / colors

ACROS = ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']
MIN_KB = int(os.environ.get('SEQ_MIN_KB', '20'))   # array length threshold (kb)
MIN_LEN = MIN_KB * 1000
TAG = f'_{MIN_KB}kb'           # tag all outputs so threshold variants coexist
OUT = sc.OUT
FIGDIR = sc.FIGDIR


def seq_str(arrays):
    """Ordered array-type sequence, e.g. 'ASat>HSat1>HSat2'."""
    return '>'.join(a[2] for a in arrays) if arrays else '(none)'


def plot_chr(chrom, sub, coord, blocks, colormap, class_colors, xmode='abs'):
    """Track plot only (no heatmap). Rows already in final (class) order."""
    keys = list(sub['key'])
    n = len(keys)
    h = max(4.0, n * 0.14)
    fig = plt.figure(figsize=(16, h))
    gs = fig.add_gridspec(
        1, 5, width_ratios=[0.16, 0.16, 0.16, 1.6, 7.0], wspace=0.03)
    ax_cls = fig.add_subplot(gs[0, 0])
    ax_sp = fig.add_subplot(gs[0, 1])
    ax_pj = fig.add_subplot(gs[0, 2])
    ax_name = fig.add_subplot(gs[0, 3])
    ax_trk = fig.add_subplot(gs[0, 4])

    def strip(ax, values, cmap, title):
        cols = np.array([[sc._hex2rgb(cmap.get(v, '#ffffff')) for v in values]])
        cols = cols.transpose(1, 0, 2)
        ax.imshow(cols, aspect='auto', interpolation='nearest')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(title, fontsize=7, rotation=90, va='bottom')

    strip(ax_cls, list(sub['seq_class']), class_colors, 'class')
    strip(ax_sp, list(sub['superpop']), sc.SUPERPOP_COLORS, 'pop')
    strip(ax_pj, list(sub['project']), sc.PROJECT_COLORS, 'proj')

    # sample names
    ax_name.set_xlim(0, 1); ax_name.set_ylim(n - 0.5, -0.5)
    pitch_pt = (h / n) * 72.0
    fs = max(4.0, min(9.0, pitch_pt * 0.72))
    for r, name in enumerate(sub['sample_hap']):
        ax_name.text(0.0, r, name, fontsize=fs, va='center', ha='left')
    ax_name.axis('off')

    # satellite track (all raw blocks incl rDNA/small), batched
    used = set()
    maxx = 0
    rects = []
    for r, key in enumerate(keys):
        cs, ce, _, _ = coord[key]
        denom = (ce - cs) or 1
        for s, e, sat, col in blocks.get(key, []):
            if xmode == 'norm':
                x = (s - cs) / denom; w = (e - s) / denom
            else:
                x = (s - cs) / 1e6; w = (e - s) / 1e6
            rects.append(Rectangle((x, r - 0.45), w, 0.9,
                                   facecolor=colormap.get(sat, '#888888'),
                                   edgecolor='none'))
            used.add(sat)
            if x + w > maxx:
                maxx = x + w
    ax_trk.add_collection(PatchCollection(rects, match_original=True))
    ax_trk.set_ylim(n - 0.5, -0.5)
    ax_trk.set_xlim(0, (1 if xmode == 'norm' else (maxx * 1.02 if maxx else 1)))
    ax_trk.set_yticks([])
    ax_trk.set_xlabel('position within centromere '
                      + ('(fraction)' if xmode == 'norm' else '(Mb)'))

    sat_handles = [Patch(facecolor=colormap.get(s, '#888888'), label=s)
                   for s in sc.SAT_ORDER if s in used]
    pop_handles = [Patch(facecolor=c, label=k)
                   for k, c in sc.SUPERPOP_COLORS.items()]
    pj_handles = [Patch(facecolor=c, label=k)
                  for k, c in sc.PROJECT_COLORS.items()]
    leg1 = ax_trk.legend(handles=sat_handles, title='satellite', fontsize=7,
                         title_fontsize=8, loc='upper left',
                         bbox_to_anchor=(1.01, 1.0))
    ax_trk.legend(handles=pop_handles + pj_handles, fontsize=6,
                  loc='upper left', bbox_to_anchor=(1.01, 0.5), ncol=2)
    ax_trk.add_artist(leg1)

    fig.suptitle(f"{chrom}  |  n={n}  |  array>={MIN_KB}kb  |  ordered-sequence "
                 f"classes={sub['seq_class'].nunique()}", fontsize=11, y=0.995)
    out = os.path.join(FIGDIR, f"{chrom}_seq{TAG}_{xmode}.pdf")
    fig.savefig(out, dpi=200, bbox_inches='tight')
    plt.close(fig)


def classify_group(chroms, tsv_name, coord, blocks, meta, colormap, selected):
    """Ordered-sequence classification for a set of chromosomes (10kb arrays).
    Writes one TSV and per-chr track PDFs; TSV row order == PDF row order."""
    rows = []
    for key, (cs, ce, proj, chrom) in coord.items():
        if chrom not in chroms or key not in selected:
            continue
        arrays = sc.segment_arrays(blocks.get(key, []), MIN_LEN)
        sample, sample_hap, superpop = meta.get(key, ('NA', key, 'NA'))
        rows.append({
            'key': key, 'chr': chrom, 'project': proj,
            'sample': sample, 'sample_hap': sample_hap, 'superpop': superpop,
            'n_arrays': len(arrays),
            'total_kb': round(sum(e - s for s, e, _ in arrays) / 1000, 1),
            'seq_class': seq_str(arrays),
        })
    df = pd.DataFrame(rows)
    df['chr'] = pd.Categorical(df['chr'], categories=chroms, ordered=True)

    # per-chr class id by frequency: {chr}_{rank}_{count}
    df = sc.assign_classes(df, 'seq_class', 'seq_label')

    # ---- order rows: per chr, most common arrangement first, then by the
    # sequence string, then sample; PDF and TSV share this exact order ----
    freq = df.groupby(['chr', 'seq_class'], observed=True)['key'].transform('size')
    df['_freq'] = freq
    df = df.sort_values(['chr', '_freq', 'seq_class', 'sample_hap'],
                        ascending=[True, False, True, True]).drop(columns='_freq')

    ordered_frames = []
    for chrom, sub in df.groupby('chr', sort=True, observed=True):
        if len(sub) == 0:
            continue
        sub = sub.reset_index(drop=True)
        # distinct-class color strip (tab20)
        uniq = list(dict.fromkeys(sub['seq_class']))
        cmap = plt.get_cmap('tab20')
        class_colors = {s: matplotlib.colors.to_hex(cmap(i % 20))
                        for i, s in enumerate(uniq)}
        plot_chr(chrom, sub, coord, blocks, colormap, class_colors, 'abs')
        plot_chr(chrom, sub, coord, blocks, colormap, class_colors, 'norm')
        sub = sub.copy()
        sub['order'] = range(1, len(sub) + 1)
        ordered_frames.append(sub)
        print(f"  {chrom}: n={len(sub)}  ordered-seq classes={sub['seq_class'].nunique()}")

    out = pd.concat(ordered_frames, ignore_index=True)
    cols = ['order', 'key', 'chr', 'sample_hap', 'project', 'superpop',
            'seq_label', 'seq_class', 'n_arrays', 'total_kb']
    tsv_path = os.path.join(OUT, tsv_name)
    sc.safe_to_csv(out[cols], tsv_path, sep='\t', index=False)
    print('  TSV ->', tsv_path)
    return out


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    coord, blocks = sc.load_blocks()
    # HSat1A/HSat1B -> HSat1 (same as the main pipeline)
    blocks = {k: [(s, e, sc.canon_sat(sat), col) for (s, e, sat, col) in recs]
              for k, recs in blocks.items()}
    meta = sc.load_sample_meta()
    selected = set(sc.select_centromeres(coord, blocks))
    colormap = sc.build_color_map(blocks)

    nonacro = [c for c in sc.CHRS if c not in ACROS]
    print(f'=== acrocentric (chr13/14/15/21/22), array>={MIN_KB}kb ===')
    a = classify_group(ACROS, f'satellite_seq_classification_acro{TAG}.tsv',
                       coord, blocks, meta, colormap, selected)
    print(f'=== non-acrocentric, array>={MIN_KB}kb ===')
    na = classify_group(nonacro, f'satellite_seq_classification_nonacro{TAG}.tsv',
                        coord, blocks, meta, colormap, selected)

    # combined master table: all 24 chromosomes in natural order, each chr's
    # rows kept in its PDF row order (the per-chr `order` column).
    allc = pd.concat([a, na], ignore_index=True)
    allc['chr'] = pd.Categorical(allc['chr'], categories=sc.CHRS, ordered=True)
    allc = allc.sort_values(['chr', 'order']).reset_index(drop=True)
    all_path = os.path.join(OUT, f'satellite_seq_classification_all{TAG}.tsv')
    sc.safe_to_csv(allc, all_path, sep='\t', index=False)
    print('  combined TSV ->', all_path)
    print('\nfigs ->', FIGDIR, f'(<chr>_seq{TAG}_abs.pdf / _seq{TAG}_norm.pdf)')


if __name__ == '__main__':
    main()
