# -*- coding: utf-8 -*-
"""
Stacked barplots of complete-centromere counts per chromosome, stacked by the
manual Final classification.  x = chromosome (chr1..22, X, Y), y = number of
sample_hap (complete centromeres).  Segment color: top1..top10 use a custom palette (deep blue → dark red),
rank ≥ 11 merged as grey "others". See RANK_COLORS dict.

Two figures:
  results/class_stacked_barplot_all.pdf   all samples
  results/class_stacked_barplot_APG.pdf   APG project only
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['pdf.fonttype'] = 42      # keep text editable in the PDF
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import satellite_classification as sc

BASE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(sc.OUT, 'satellite_seq_classification_all_final_fixed.tsv')
CHRS = sc.CHRS
# New color palette: top1..top10 + others
RANK_COLORS = {
    1:  "#1a3951",
    2:  "#087D6F",
    3:  "#15959F",
    4:  "#F1E4B3",
    5:  "#FFC753",
    6:  "#FFB1A8",
    7:  "#EC9770",
    8:  "#C7402D",
    9:  "#D9430D",
    10: "#A60A0A",
}
OTHERS_GRAY = 'grey'
BAR_WIDTH = 0.6                 # thinner bars


def rank_color(idx):
    """Return color for a given rank (1-based, 1..10)."""
    return RANK_COLORS.get(idx, OTHERS_GRAY)


FIGSIZE = (8.27, 5.3)          # A4 width, flattened (shorter) height


def stacked_barplot(df, title, out):
    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes([0.095, 0.15, 0.75, 0.74])   # wide, short plot
    maxidx, maxtotal = 0, 0
    has_others = False
    for x, chrom in enumerate(CHRS):
        sub = df[df['chr'] == chrom]
        if len(sub) == 0:
            continue
        order = sorted(sub['Final'].unique(),
                       key=lambda s: int(str(s).split('_')[-2]))
        bottom, others = 0, 0
        for lab in order:
            cnt = int((sub['Final'] == lab).sum())
            if cnt == 0:
                continue
            idx = int(str(lab).split('_')[-2])
            maxidx = max(maxidx, idx)
            if idx <= 10:                            # Arr-1..10 keep their color
                ax.bar(x, cnt, bottom=bottom, width=BAR_WIDTH,
                       color=rank_color(idx), edgecolor='none')
                bottom += cnt
            else:                                    # Arr>=11 -> grey "others"
                others += cnt
        if others > 0:
            ax.bar(x, others, bottom=bottom, width=BAR_WIDTH, color=OTHERS_GRAY,
                   edgecolor='none')
            bottom += others
            has_others = True
        maxtotal = max(maxtotal, bottom)
        ax.text(x, bottom + 5, str(bottom), ha='center', va='bottom',
                fontsize=7)                      # horizontal counts

    ax.set_ylim(0, maxtotal * 1.10)              # small headroom for the counts
    ax.set_xlim(-0.7, len(CHRS) - 0.3)
    ax.set_xticks(range(len(CHRS)))
    ax.set_xticklabels([c.replace('chr', '') for c in CHRS])
    ax.tick_params(axis='both', labelsize=9)
    ax.set_xlabel('chromosome', fontsize=12)
    ax.set_ylabel('The number of complete centromeres', fontsize=12)
    ax.set_title(title, fontsize=12, pad=8)
    ax.spines[['top', 'right']].set_visible(False)

    # legend (right, vertically centered): CenArch-1..10 + grey "others"
    handles = [Patch(facecolor=rank_color(i + 1), label=f'CenArch-{i + 1}')
               for i in range(min(10, maxidx))]
    if has_others:
        handles.append(Patch(facecolor=OTHERS_GRAY, label='others'))
    ax.legend(handles=handles, title='centromere satellite\narrangement',
              fontsize=7.5, title_fontsize=8, ncol=1, loc='center left',
              bbox_to_anchor=(1.02, 0.5), handlelength=1.2, handleheight=1.0,
              labelspacing=0.3)
    fig.savefig(out, dpi=200)            # keep exact A4 canvas (no tight crop)
    plt.close(fig)
    print('->', out)


def main():
    df = pd.read_csv(IN, sep='\t', dtype={'Final': str})
    df['chr'] = pd.Categorical(df['chr'], categories=CHRS, ordered=True)

    stacked_barplot(df, 'Centromere architectures per chromosome (all samples)',
                    os.path.join(sc.OUT, 'class_stacked_barplot_all_v2.pdf'))
    apg = df[df['project'] == 'APG']
    stacked_barplot(apg, 'Centromere architectures per chromosome (APGp1 only)',
                    os.path.join(sc.OUT, 'class_stacked_barplot_APG_v2.pdf'))
    print(f"all samples: {len(df)} centromeres; APG: {len(apg)}")


if __name__ == '__main__':
    main()
