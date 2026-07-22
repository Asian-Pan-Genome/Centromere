# -*- coding: utf-8 -*-
"""
Per-chromosome "representative example" figure: for every Final classification,
show ONE example centromere (its real satellite track in absolute Mb), with the
class size (n = all | n = APG) and two pie charts (superpopulation & project
composition of that whole class).

Representative pick per class:
  1. if the class contains Ref samples -> first available in priority order
     CHM13_-, CN1_Mat, CN1_Pat, HG002_Mat, HG002_Pat, YAO_Mat, YAO_Pat;
  2. else a random APG sample;
  3. else a random HPRC/HGSVC sample.

Row order:  classes WITH APG samples first (by class size desc), then classes
without APG (only Ref/HPRC/HGSVC, by class size desc).

Class name:  CenArch-{rank}  (rank 1 = largest class); the chromosome is shown
once in the figure header (e.g. CEN01), so each row only needs CenArch-{rank}.

Layout: 2026-07 AI spec — A4 portrait, fixed row pitch 0.261in, track height
0.12in, max track width 4.977in.  See results/example_figure/chr4_examples_AI_modifications.md

One PDF per chromosome -> results/example_figure/<chr>_example.pdf
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['pdf.fonttype'] = 42          # editable text in the PDF
# route mathtext (italic n, bold section titles) to Arial too
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.collections import PatchCollection

import satellite_classification as sc

IN = os.path.join(sc.OUT, 'satellite_seq_classification_all_final_fixed.tsv')
OUTDIR = os.path.join(sc.OUT, 'example_figure')
CHRS = sc.CHRS
SEED = 0
REF_PRIORITY = ['CHM13_-', 'CN1_Mat', 'CN1_Pat', 'HG002_Mat', 'HG002_Pat',
                'YAO_Mat', 'YAO_Pat']
SP_ORDER = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
PJ_ORDER = ['APG', 'HGSVC', 'HPRC', 'Ref']

# ── new rank color palette (CenArch-1..10 + others grey) ──────────────────
RANK_COLORS = {
    1:  "#1a3951", 2:  "#087D6F", 3:  "#15959F", 4:  "#F1E4B3",
    5:  "#FFC753", 6:  "#FFB1A8", 7:  "#EC9770", 8:  "#C7402D",
    9:  "#D9430D", 10: "#A60A0A",
}
OTHERS_GRAY = 'grey'

def rank_color(idx):
    return RANK_COLORS.get(idx, OTHERS_GRAY)


# ── AI spec layout constants (all in inches, top-left origin) ──────────────
FIG_W = 8.27       # A4 width
FIG_H = 11.69      # A4 height

# Header (positions from chr4 AI spec — left-edge of text, va='top')
# All Y values adjusted +0.018″ to compensate for matplotlib va='top' vs AI bbox offset
Y_OFFSET     = 0.018  # systematic offset between mpl va='top' and AI text bbox top
CEN_Y        = 0.228 + Y_OFFSET   # CEN{NN} label top
CEN_X        = 0.663              # left edge
N_HEADER_Y   = 0.401 + Y_OFFSET   # n(all) | n(APG) header top
N_HEADER_X   = 0.609              # first italic 'n' left edge
SUPER_HDR_X  = 1.403              # "Super" / "population" left edge
SUPER_HDR_Y1 = 0.278 + Y_OFFSET
SUPER_HDR_Y2 = 0.326 + Y_OFFSET
PROJ_HDR_X   = 1.832              # "Project" left edge
PROJ_HDR_Y   = 0.326 + Y_OFFSET

# Row template (row 1 baseline)
ROW_PITCH      = 0.261   # vertical spacing between rows
TRACK_HEIGHT   = 0.12    # track bar height
TRACK_Y0       = 0.585   # top of first track row
TRACK_X0       = 2.161   # left edge of track
TRACK_MAX_W    = 4.977   # max track width (chr4 longest)

CENARCH_X      = 0.66    # CenArch label x
CENARCH_Y0     = 0.520 + Y_OFFSET   # CenArch label y (row 1)
COLOR_SQ_X     = 0.547   # colour square left edge
COLOR_SQ_SIZE  = 4.4/72  # 4.4 pt -> inches
N_VALUE_Y0     = 0.636 + Y_OFFSET   # n = … y (row 1)
SAMPLE_NAME_Y0 = 0.478 + Y_OFFSET   # sample name y (row 1)
SAMPLE_NAME_X  = 2.161   # sample name x (= track left)

# Pies
PIE_SIZE       = 0.25    # pie axes square side (inches)
SP_PIE_CX      = 1.555   # superpop pie centre x
PJ_PIE_CX      = 1.96    # project pie centre x

# Legend
LEGEND_TITLE_X = 7.17   # bbox_to_anchor; actual text appears ~0.19in to the right
LEGEND_ENTRY_X = 7.36
LEGEND_SWATCH_X = 7.24   # colour swatch left

# X-axis
XAXIS_LABEL_Y  = 2.42    # "Centromere length (Mbp)" y (row 1 base)
TEXT_COLOR      = '#221f1f'


def inch_to_fig(x_inch, y_inch):
    """Convert inch coords (top-left origin) to figure coords (bottom-left)."""
    return x_inch / FIG_W, 1.0 - y_inch / FIG_H


def add_axes_inch(fig, left, top, width, height):
    """Add axes using inch coords (top-left origin for top)."""
    return fig.add_axes([
        left / FIG_W,
        1.0 - (top + height) / FIG_H,
        width / FIG_W,
        height / FIG_H,
    ])


def cen_prefix(chrom):
    c = chrom[3:] if chrom.startswith('chr') else chrom
    if c.isdigit():
        c = f"{int(c):02d}"
    return f"CEN{c}"


def cen_name(chrom, idx):
    return f"CenArch-{idx}"


def disp_hap(name, project):
    """Format sample name for display according to AI spec."""
    if project == 'APG':
        # C038-CHA-S18_Mat -> C038-CHA-S18-01#Mat
        if '_' in name:
            sample, hap = name.rsplit('_', 1)
            return f"{sample}-01#{hap}"
        return name
    elif project == 'Ref':
        if name == 'CHM13_-':
            return 'CHM13'
        # CN1_Mat -> T2T-CN1#Mat
        if '_' in name:
            sample, hap = name.rsplit('_', 1)
            return f"T2T-{sample}#{hap}"
        return f"T2T-{name}"
    else:
        # HGSVC, HPRC: keep original
        return name


def pick_rep(members):
    ref = members[members['project'] == 'Ref']
    if len(ref):
        for name in REF_PRIORITY:
            hit = ref[ref['sample_hap'] == name]
            if len(hit):
                return hit.iloc[0]
        return ref.iloc[0]
    apg = members[members['project'] == 'APG']
    if len(apg):
        return apg.sample(1, random_state=SEED).iloc[0]
    oth = members[members['project'].isin(['HPRC', 'HGSVC'])]
    if len(oth):
        return oth.sample(1, random_state=SEED).iloc[0]
    return members.iloc[0]


def pie(ax, counts, order, cmap):
    """Draw a small pie chart with no border/frame."""
    labels = [k for k in order if counts.get(k, 0) > 0]
    sizes = [counts[k] for k in labels]
    colors = [cmap.get(k, '#cccccc') for k in labels]
    ax.axis('off')
    if not sizes:
        return
    ax.pie(sizes, colors=colors, startangle=90, radius=1.0,
           wedgeprops=dict(edgecolor='none', linewidth=0))
    ax.set_aspect('equal')


def fig_for_chr(chrom, sub, coord, blocks, colormap, outdir, geom):
    # build per-class info, then order rows
    rows = []
    for lab in sub['Final'].unique():
        members = sub[sub['Final'] == lab]
        idx = int(str(lab).split('_')[-2])
        n_all = len(members)
        n_apg = int((members['project'] == 'APG').sum())
        rows.append((lab, idx, n_all, n_apg, n_apg > 0, members))
    # APG-containing first (count desc), then non-APG (count desc)
    rows.sort(key=lambda t: (not t[4], -t[2], t[1]))
    k = len(rows)

    # representative per class
    reps = {lab: pick_rep(members) for lab, *_rest, members in rows}
    # x-axis unified WITHIN this figure = this chromosome's longest example
    maxlen_mb = max((coord[r['key']][1] - coord[r['key']][0])
                    for r in reps.values()) / 1e6 or 1

    # A4 fixed canvas
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    used = set()

    for r, (lab, idx, n_all, n_apg, has_apg, members) in enumerate(rows):
        rep = reps[lab]
        name = cen_name(chrom, idx)
        row_y = r * ROW_PITCH  # offset from row 1 baseline

        # ── colour square (marker) ──
        sq_top = CENARCH_Y0 + row_y
        sq_y_fig = 1.0 - sq_top / FIG_H - COLOR_SQ_SIZE / FIG_H / 2.0
        sq_x_fig = COLOR_SQ_X / FIG_W - COLOR_SQ_SIZE / FIG_W / 2.0
        fig.patches.append(Rectangle(
            (sq_x_fig, sq_y_fig), COLOR_SQ_SIZE / FIG_W, COLOR_SQ_SIZE / FIG_H,
            facecolor=rank_color(idx), edgecolor='none',
            transform=fig.transFigure, clip_on=False))

        # ── CenArch label ──
        fig.text(CENARCH_X / FIG_W, 1.0 - (CENARCH_Y0 + row_y) / FIG_H,
                 name, fontsize=8, va='top', ha='left',
                 fontweight='normal', color=TEXT_COLOR)

        # ── n = … | n = … (both 'n' italic, rest regular) ──
        ny = N_VALUE_Y0 + row_y
        y_fig = 1.0 - ny / FIG_H
        x0 = CENARCH_X / FIG_W
        # italic 'n'
        fig.text(x0, y_fig, 'n', fontsize=6, va='top', ha='left',
                 fontstyle='italic', color=TEXT_COLOR)
        # regular "= NNN  |  "
        x1 = x0 + 0.0075
        fig.text(x1, y_fig, f'= {n_all}  |  ', fontsize=6, va='top',
                 ha='left', fontstyle='normal', color=TEXT_COLOR)
        # italic 'n'
        x2 = x1 + 0.035
        fig.text(x2, y_fig, 'n', fontsize=6, va='top', ha='left',
                 fontstyle='italic', color=TEXT_COLOR)
        # regular "= NNN"
        x3 = x2 + 0.0075
        fig.text(x3, y_fig, f'= {n_apg}', fontsize=6, va='top',
                 ha='left', fontstyle='normal', color=TEXT_COLOR)

        # ── sample name above track ──
        sn_y = SAMPLE_NAME_Y0 + row_y
        fig.text(SAMPLE_NAME_X / FIG_W, 1.0 - sn_y / FIG_H,
                 disp_hap(rep['sample_hap'], rep['project']),
                 fontsize=7, va='top', ha='left', color=TEXT_COLOR)

        # ── superpopulation pie ──
        track_center_y = TRACK_Y0 + row_y + TRACK_HEIGHT / 2.0
        ax_sp = add_axes_inch(fig, SP_PIE_CX - PIE_SIZE / 2.0,
                              track_center_y - PIE_SIZE / 2.0,
                              PIE_SIZE, PIE_SIZE)
        pie(ax_sp, members['superpop'].value_counts().to_dict(),
            SP_ORDER, sc.SUPERPOP_COLORS)

        # ── project pie ──
        ax_pj = add_axes_inch(fig, PJ_PIE_CX - PIE_SIZE / 2.0,
                              track_center_y - PIE_SIZE / 2.0,
                              PIE_SIZE, PIE_SIZE)
        pie(ax_pj, members['project'].value_counts().to_dict(),
            PJ_ORDER, sc.PROJECT_COLORS)

        # ── satellite track ──
        track_top = TRACK_Y0 + row_y
        ax_t = add_axes_inch(fig, TRACK_X0, track_top,
                             TRACK_MAX_W, TRACK_HEIGHT)
        ax_t.set_facecolor('none')
        cs, ce, _, _ = coord[rep['key']]
        cen_mb = (ce - cs) / 1e6
        # scale: TRACK_MAX_W inches = maxlen_mb Mbp
        scale = TRACK_MAX_W / maxlen_mb
        rects = []
        for s, e, sat, col in blocks.get(rep['key'], []):
            x = (s - cs) / 1e6 * scale
            w = max((e - s) / 1e6 * scale, 0.001)  # ensure visible
            rects.append(Rectangle((x, 0.0), w, 1.0,
                                   facecolor=colormap.get(sat, '#888888'),
                                   edgecolor='none'))
            used.add(sat)
        ax_t.add_collection(PatchCollection(rects, match_original=True))
        # black outer border
        ax_t.add_patch(Rectangle((0, 0), cen_mb * scale, 1.0, fill=False,
                                 edgecolor='black', linewidth=0.25))
        ax_t.set_xlim(0, TRACK_MAX_W)
        ax_t.set_ylim(0, 1)
        ax_t.set_yticks([])
        for sp in ['top', 'right', 'left']:
            ax_t.spines[sp].set_visible(False)

        if r != k - 1:
            # not last row: hide bottom spine + ticks
            ax_t.set_xticks([])
            ax_t.spines['bottom'].set_visible(False)
        else:
            # last row: show x-axis
            ax_t.set_xlabel('Centromere length (Mbp)', fontsize=8,
                            color=TEXT_COLOR)
            ax_t.tick_params(labelsize=8, colors=TEXT_COLOR)
            # set nice tick positions
            n_ticks = 5
            tick_step = maxlen_mb / (n_ticks - 1)
            tick_positions = [i * tick_step * scale for i in range(n_ticks)]
            tick_labels = [f"{i * tick_step:.0f}" for i in range(n_ticks)]
            ax_t.set_xticks(tick_positions)
            ax_t.set_xticklabels(tick_labels)
            ax_t.spines['bottom'].set_visible(True)
            ax_t.spines['bottom'].set_color(TEXT_COLOR)

    # ── header texts (4-piece n(all) | n(APG) matching chr4 AI spec) ──
    # CEN{NN}
    fig.text(CEN_X / FIG_W, 1.0 - CEN_Y / FIG_H, cen_prefix(chrom),
             fontsize=8, fontweight='bold', ha='left', va='top',
             color=TEXT_COLOR)
    # n(all) | n(APG): italic 'n' + regular '(all)' + italic '  |  n' + regular '(APG)'
    hy = 1.0 - N_HEADER_Y / FIG_H
    hx = N_HEADER_X / FIG_W
    fig.text(hx, hy, 'n', fontsize=6, fontstyle='italic', ha='left', va='top',
             color=TEXT_COLOR)
    hx += 0.0055                                 # italic 'n' width ≈ 0.046″
    fig.text(hx, hy, '(all)', fontsize=6, fontstyle='normal', ha='left', va='top',
             color=TEXT_COLOR)
    hx += 0.0168                                 # '(all)' width ≈ 0.139″
    fig.text(hx, hy, '  |  n', fontsize=6, fontstyle='italic', ha='left', va='top',
             color=TEXT_COLOR)
    hx += 0.0195                                 # '  |  n' width ≈ 0.161″
    fig.text(hx, hy, '(APG)', fontsize=6, fontstyle='normal', ha='left', va='top',
             color=TEXT_COLOR)
    # Superpopulation header (two lines, left-aligned)
    fig.text(SUPER_HDR_X / FIG_W, 1.0 - SUPER_HDR_Y1 / FIG_H,
             'Super', fontsize=6, ha='left', va='top', color=TEXT_COLOR)
    fig.text(SUPER_HDR_X / FIG_W, 1.0 - SUPER_HDR_Y2 / FIG_H,
             'population', fontsize=6, ha='left', va='top', color=TEXT_COLOR)
    # Project header (left-aligned)
    fig.text(PROJ_HDR_X / FIG_W, 1.0 - PROJ_HDR_Y / FIG_H,
             'Project', fontsize=6, ha='left', va='top', color=TEXT_COLOR)

    # ── combined legend on the right ──
    blank = Patch(facecolor='none', edgecolor='none')
    pj_label = {'APG': 'APGp1', 'HGSVC': 'HGSVC3', 'HPRC': 'HPRCy1', 'Ref': 'Ref'}
    handles, labels, title_rows = [], [], []

    def section(title, items):
        if handles:
            handles.append(blank); labels.append(' ')
        title_rows.append(len(labels))
        handles.append(blank); labels.append(title)
        for h_, l_ in items:
            handles.append(h_); labels.append(l_)

    section('Satellite', [(Patch(facecolor=colormap.get(s, '#888888')), s)
                          for s in sc.SAT_ORDER if s != 'rDNA'])
    section('Superpopulation', [(Patch(facecolor=c), k_)
                                for k_, c in sc.SUPERPOP_COLORS.items()])
    section('Project', [(Patch(facecolor=c), pj_label.get(k_, k_))
                        for k_, c in sc.PROJECT_COLORS.items()])

    leg = fig.legend(handles, labels, loc='center left',
                     bbox_to_anchor=(LEGEND_TITLE_X / FIG_W, 0.5),
                     fontsize=5, frameon=False,
                     handlelength=1.1, labelspacing=0.35)
    leg_texts = leg.get_texts()
    for i in title_rows:
        leg_texts[i].set_fontweight('bold')
        leg_texts[i].set_fontsize(6)
    # Set legend text color
    for txt in leg_texts:
        txt.set_color(TEXT_COLOR)

    out = os.path.join(outdir, f"{chrom}_example.pdf")
    try:
        fig.savefig(out, dpi=200)
    except PermissionError:
        print(f"  [WARN] {os.path.basename(out)} is open/locked -> skipped "
              f"(close it and re-run to refresh)")
    plt.close(fig)
    print(f"  {chrom}: {k} classes -> {os.path.basename(out)}")


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    only = set(sys.argv[1:]) if len(sys.argv) > 1 else None
    df = pd.read_csv(IN, sep='\t', dtype={'Final': str})
    df['chr'] = pd.Categorical(df['chr'], categories=CHRS, ordered=True)
    coord, blocks = sc.load_blocks()
    blocks = {kk: [(s, e, sc.canon_sat(sat), col) for (s, e, sat, col) in recs]
              for kk, recs in blocks.items()}
    colormap = sc.build_color_map(blocks)

    # Fixed A4 canvas — all chromosomes same size
    geom = {}  # kept for compatibility; all layout now uses module constants
    print(f"fixed canvas: {FIG_W:.1f} x {FIG_H:.1f} in (A4) | row_pitch={ROW_PITCH:.3f}in")

    for chrom, sub in df.groupby('chr', sort=True, observed=True):
        if len(sub) == 0 or (only and str(chrom) not in only):
            continue
        fig_for_chr(str(chrom), sub.reset_index(drop=True), coord, blocks,
                    colormap, OUTDIR, geom)
    print('\nfigures ->', OUTDIR, '(<chr>_example.pdf)')


if __name__ == '__main__':
    main()
