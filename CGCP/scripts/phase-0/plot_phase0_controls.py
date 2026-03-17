#!/usr/bin/env python3
"""
Phase 0 Ground Truth Validation — GraphPad Prism-style figures
BCL-XL/BIM (druggable), MDM2/p53 (druggable), IL-2/IL-2Ra (undruggable)

Run:
    python3 plot_phase0_controls.py

Output:
    ~/projects/rtc-pan-coronavirus/CGCP/00-ground-truth/visuals/
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import warnings
warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════
# PRISM STYLE BASE
# ═══════════════════════════════════════════════════════════════════
plt.rcParams.update({
    # Font
    'font.family':          'sans-serif',
    'font.sans-serif':      ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':            9,

    # Axes
    'axes.linewidth':       0.75,
    'axes.spines.top':      False,
    'axes.spines.right':    False,
    'axes.spines.left':     True,
    'axes.spines.bottom':   True,
    'axes.facecolor':       'white',
    'axes.labelsize':       9,
    'axes.labelpad':        6,       # more breathing room
    'axes.titlesize':       9,
    'axes.titleweight':     'normal',
    'axes.titlepad':        10,

    # Ticks — outward, Prism default
    'xtick.direction':      'out',
    'ytick.direction':      'out',
    'xtick.major.size':     4,
    'ytick.major.size':     4,
    'xtick.major.width':    0.75,
    'ytick.major.width':    0.75,
    'xtick.minor.size':     2.5,
    'ytick.minor.size':     2.5,
    'xtick.labelsize':      8.5,
    'ytick.labelsize':      8.5,
    'xtick.major.pad':      4,
    'ytick.major.pad':      4,

    # No grid
    'axes.grid':            False,

    # Legend
    'legend.fontsize':      8,
    'legend.frameon':       False,
    'legend.handlelength':  1.2,
    'legend.handleheight':  0.9,
    'legend.labelspacing':  0.4,

    # Figure
    'figure.facecolor':     'white',
    'figure.dpi':           150,
    'savefig.dpi':          300,
    'savefig.bbox':         'tight',
    'savefig.facecolor':    'white',
    'savefig.pad_inches':   0.12,    # outer margin

    # Editable text
    'pdf.fonttype':         42,
    'ps.fonttype':          42,
})


# ═══════════════════════════════════════════════════════════════════
# PALETTE
# ═══════════════════════════════════════════════════════════════════
P = {
    'bcl2':   '#2166AC',
    'mdm2':   '#4DAC26',
    'il2':    '#D7191C',
    'bcl2_f': '#92C5DE',
    'mdm2_f': '#A6D96A',
    'il2_f':  '#FDAE61',
    'gray':   '#636363',
    'lgray':  '#BBBBBB',
    'text':   '#1A1A1A',
}

SYSTEMS = ['BCL-XL/BIM', 'MDM2/p53', 'IL-2/IL-2R\u03b1']
COLORS  = [P['bcl2'], P['mdm2'], P['il2']]
FILLS   = [P['bcl2_f'], P['mdm2_f'], P['il2_f']]
LABELS  = ['Druggable', 'Druggable', 'Undruggable']


# ═══════════════════════════════════════════════════════════════════
# PRISM AXIS TRIM — spine ends at last tick, offset from plot
# ═══════════════════════════════════════════════════════════════════
def prism_axes(ax, ymax=None, yticks=None, xticks=None,
               x_label_rot=0):
    """Trim spines to data range. Ticks outward. No grid."""
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    if ymax is not None:
        ax.spines['left'].set_bounds(0, ymax)
        ax.set_ylim(-ymax * 0.01, ymax * 1.01)
    if yticks is not None:
        ax.set_yticks(yticks)
        ax.set_yticklabels([str(t) for t in yticks])
    if xticks is not None:
        ax.set_xticks(xticks)
    if x_label_rot:
        ax.set_xticklabels(ax.get_xticklabels(),
                           rotation=x_label_rot, ha='right')


def prism_axes_h(ax, xmax=None, xticks=None):
    """Horizontal bar variant."""
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    if xmax is not None:
        ax.spines['bottom'].set_bounds(0, xmax)
        ax.set_xlim(-xmax * 0.01, xmax * 1.01)
    if xticks is not None:
        ax.set_xticks(xticks)


# ═══════════════════════════════════════════════════════════════════
# DATA
# ═══════════════════════════════════════════════════════════════════
def load_tsv(path):
    if not os.path.exists(path):
        print(f'  Warning: {path} not found')
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))

def count_res(contacts, key='res_a'):
    c = {}
    for r in contacts:
        c[r[key]] = c.get(r[key], 0) + 1
    return c

def dom_ratio(counts):
    if not counts:
        return 0.0
    v = list(counts.values())
    return round(max(v) / sum(v), 4)

BASE = os.path.expanduser(
    '~/projects/rtc-pan-coronavirus/CGCP/00-ground-truth')
OUT  = os.path.join(BASE, 'visuals')
os.makedirs(OUT, exist_ok=True)

bc = count_res(load_tsv(os.path.join(
    BASE, 'positive-controls/bcl2-bax/bcl2_contacts.tsv')), 'res_a')
mb = count_res(load_tsv(os.path.join(
    BASE, 'positive-controls/mdm2-p53/mdm2_contacts.tsv')), 'res_b')
ic = count_res(load_tsv(os.path.join(
    BASE, 'negative-controls/il2_contacts.tsv')), 'res_a')


# ═══════════════════════════════════════════════════════════════════
# FIG 1 — Dominance Ratio
# Wide figure → x-labels never rotate into bars
# ═══════════════════════════════════════════════════════════════════
def fig1():
    # Wider figure gives x-labels room without rotation
    fig, ax = plt.subplots(figsize=(4.4, 3.6))
    fig.subplots_adjust(bottom=0.18, top=0.88, left=0.17, right=0.88)

    dr  = [dom_ratio(bc), dom_ratio(mb), dom_ratio(ic)]
    x   = np.arange(3)
    w   = 0.40       # narrow Prism bars

    for i, (xi, val) in enumerate(zip(x, dr)):
        ax.bar(xi, val, width=w,
               facecolor=FILLS[i], edgecolor=COLORS[i],
               linewidth=0.9, zorder=3,
               clip_on=False)
        # Value label: fixed 8% of ymax above bar top
        ax.text(xi, val + 0.013,
                f'{val:.3f}',
                ha='center', va='bottom',
                fontsize=8, color=P['text'],
                zorder=5)

    # Threshold line — drawn before bars so it sits behind
    ax.axhline(0.15,
               color=P['gray'], linewidth=0.8,
               linestyle=(0, (5, 3)),
               zorder=1)
    # Threshold label — right of plot area, not on top of bars
    ax.text(2.24, 0.148,
            'threshold = 0.15',
            fontsize=7, color=P['gray'],
            va='top', ha='left',
            clip_on=False)

    # X-axis
    ax.set_xticks(x)
    ax.set_xticklabels(SYSTEMS, fontsize=8.5)
    ax.set_xlim(-0.6, 2.6)

    # Y-axis
    yticks = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
    ax.set_ylabel('Dominance ratio', fontsize=9)
    prism_axes(ax, ymax=0.30, yticks=yticks)

    # Panel label — outside axes, top-left
    ax.text(-0.15, 1.06, 'A',
            transform=ax.transAxes,
            fontsize=12, fontweight='bold', color=P['text'])

    p = os.path.join(OUT, 'Fig1_DominanceRatio.png')
    fig.savefig(p)
    plt.close()
    print(f'  Saved: {p}')


# ═══════════════════════════════════════════════════════════════════
# FIG 2 — Interface profile: residues + avg contacts
# Extra bottom margin keeps x-labels from touching bars
# ═══════════════════════════════════════════════════════════════════
def fig2():
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.8))
    fig.subplots_adjust(wspace=0.50,
                        bottom=0.22, top=0.88,
                        left=0.12, right=0.96)

    x = np.arange(3)
    w = 0.40

    for ax, vals, ymax, yticks, ylabel, panel in [
        (axes[0],
         [len(bc), len(mb), len(ic)],
         32, [0, 8, 16, 24, 32],
         'Interface residues (chain A)', 'A'),
        (axes[1],
         [round(sum(bc.values()) / len(bc), 1) if bc else 0,
          round(sum(mb.values()) / len(mb), 1) if mb else 0,
          round(sum(ic.values()) / len(ic), 1) if ic else 0],
         4.0, [0, 1.0, 2.0, 3.0, 4.0],
         'Avg contacts per residue', 'B'),
    ]:
        for i, (xi, val) in enumerate(zip(x, vals)):
            ax.bar(xi, val, width=w,
                   facecolor=FILLS[i], edgecolor=COLORS[i],
                   linewidth=0.9, zorder=3, clip_on=False)
            # Value label gap = 2% of ymax
            gap = ymax * 0.02
            fmt = str(val) if isinstance(val, int) else f'{val:.1f}'
            ax.text(xi, val + gap, fmt,
                    ha='center', va='bottom',
                    fontsize=8, color=P['text'], zorder=5)

        ax.set_xticks(x)
        ax.set_xticklabels(SYSTEMS, fontsize=8.5)
        ax.set_xlim(-0.6, 2.6)
        ax.set_ylabel(ylabel, fontsize=9)
        prism_axes(ax, ymax=ymax, yticks=yticks)
        ax.text(-0.20, 1.06, panel,
                transform=ax.transAxes,
                fontsize=12, fontweight='bold', color=P['text'])

    # Shared legend below both panels
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=FILLS[i], edgecolor=COLORS[i],
                     linewidth=0.9,
                     label=f'{SYSTEMS[i]}  ({LABELS[i]})')
               for i in range(3)]
    fig.legend(handles=handles,
               loc='lower center', ncol=3,
               fontsize=7.5, frameon=False,
               bbox_to_anchor=(0.54, 0.00))

    p = os.path.join(OUT, 'Fig2_InterfaceProfile.png')
    fig.savefig(p)
    plt.close()
    print(f'  Saved: {p}')


# ═══════════════════════════════════════════════════════════════════
# FIG 3 — Hotspot Detection (horizontal bars)
# Right-side labels use fixed offset, no overlap with bars
# ═══════════════════════════════════════════════════════════════════
def fig3():
    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.2))
    fig.subplots_adjust(wspace=0.60,
                        bottom=0.20, top=0.85,
                        left=0.14, right=0.96)

    pairs = [
        (axes[0], bc, ['PHE105', 'VAL126', 'PHE146'],
         P['bcl2'], P['bcl2_f'], 'BCL-XL / BIM', '3FDL', 'A'),
        (axes[1], mb, ['PHE19',  'TRP23',  'LEU26'],
         P['mdm2'], P['mdm2_f'], 'MDM2 / p53',   '1YCR', 'B'),
    ]

    for ax, counts, hotspots, edge, fill, system, pdb, panel in pairs:
        vals  = [counts.get(h, 0) for h in hotspots]
        y     = np.arange(len(hotspots))
        xmax  = max(vals)
        # x-axis limit: bar + label region; label region = 45% of xmax
        xlim  = xmax * 1.55

        for i, (yi, val) in enumerate(zip(y, vals)):
            ax.barh(yi, val, height=0.40,
                    facecolor=fill, edgecolor=edge,
                    linewidth=0.9, zorder=3, clip_on=False)
            # Label: fixed 8% of xmax right of bar end
            ax.text(val + xmax * 0.08, yi,
                    f'{val} contacts  \u2713',
                    va='center', ha='left',
                    fontsize=8.5, color=edge, zorder=5)

        # Tick spacing
        step   = max(1, int(xmax / 4))
        xticks = list(range(0, xmax + 1, step))

        ax.set_yticks(y)
        ax.set_yticklabels(hotspots, fontsize=9.5)
        ax.set_xlabel('Number of contacts', fontsize=9)
        ax.set_xlim(-xmax * 0.03, xlim)
        prism_axes_h(ax, xmax=xmax, xticks=xticks)
        ax.invert_yaxis()

        # Detection badge — bottom-right, clear of bars
        n = sum(1 for v in vals if v > 0)
        ax.text(0.97, 0.04,
                f'{n}/{len(hotspots)} detected',
                transform=ax.transAxes,
                ha='right', va='bottom',
                fontsize=7.5, color=edge,
                bbox=dict(boxstyle='round,pad=0.28',
                          facecolor=fill + 'CC',
                          edgecolor=edge,
                          linewidth=0.6))

        # Title
        ax.set_title(f'{panel}   {system}',
                     loc='left', fontsize=9.5, pad=8)

        # Subtitle below x-axis label — extra bottom margin handles this
        ax.text(0.5, -0.30,
                f'PDB: {pdb}  \u2014  known PPI anchor residues',
                transform=ax.transAxes,
                ha='center', fontsize=7.5,
                color=P['gray'], style='italic')

        ax.text(-0.20, 1.08, panel,
                transform=ax.transAxes,
                fontsize=12, fontweight='bold', color=P['text'])

    p = os.path.join(OUT, 'Fig3_HotspotDetection.png')
    fig.savefig(p)
    plt.close()
    print(f'  Saved: {p}')


# ═══════════════════════════════════════════════════════════════════
# FIG 4 — Contact Distribution Profiles
# Hotspot labels rotated; mean label in right margin not over bars
# ═══════════════════════════════════════════════════════════════════
def fig4():
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 4.0))
    fig.subplots_adjust(wspace=0.48,
                        bottom=0.14, top=0.84,
                        left=0.08, right=0.97)

    datasets = [
        (bc, P['bcl2'], P['bcl2_f'],
         ['PHE105', 'VAL126', 'PHE146'],
         'BCL-XL / BIM', 'A', 'Druggable'),
        (mb, P['mdm2'], P['mdm2_f'],
         ['PHE19', 'TRP23', 'LEU26'],
         'MDM2 / p53',   'B', 'Druggable'),
        (ic, P['il2'],  P['il2_f'],
         [],
         'IL-2 / IL-2R\u03b1', 'C', 'Undruggable'),
    ]

    for ax, (counts, edge, fill, hotspots, name, panel, label) \
            in zip(axes, datasets):
        if not counts:
            continue

        items  = sorted(counts.items(), key=lambda z: z[1], reverse=True)
        res    = [r for r, _ in items]
        vals   = np.array([v for _, v in items])
        x      = np.arange(len(vals))
        vmax   = int(vals.max())

        # ── Bars ──────────────────────────────────────────────────
        for i, (xi, val, r) in enumerate(zip(x, vals, res)):
            fc = edge if r in hotspots else fill
            ax.bar(xi, val, width=0.70,
                   facecolor=fc, edgecolor='none',
                   linewidth=0, zorder=3)

        # ── Hotspot residue labels ─────────────────────────────────
        # Stagger: odd-index hotspots sit a bit higher to avoid crowding
        for idx, (r, val) in enumerate(items):
            if r in hotspots:
                extra = vmax * 0.10 if idx % 2 == 1 else 0
                ax.text(idx,
                        val + vmax * 0.12 + extra,
                        r,
                        ha='center', va='bottom',
                        fontsize=6.5, rotation=40,
                        rotation_mode='anchor',
                        color=P['text'], zorder=6,
                        clip_on=False)

        # ── Mean dashed line ─────────────────────────────────────
        mean_v = vals.mean()
        ax.axhline(mean_v,
                   color=P['lgray'], linewidth=0.8,
                   linestyle=(0, (5, 3)),
                   zorder=2)
        # Mean label: right edge of axes, in top margin — never over bars
        ax.text(1.01, mean_v / (vmax * 1.65),
                f'\u03bc = {mean_v:.1f}',
                transform=ax.get_yaxis_transform(),
                ha='left', va='center',
                fontsize=7, color=P['gray'],
                clip_on=False)

        # ── Dominance ratio — top-right corner ────────────────────
        dr = dom_ratio(counts)
        ax.text(0.97, 0.95,
                f'ratio = {dr:.3f}',
                transform=ax.transAxes,
                ha='right', va='top',
                fontsize=8, color=edge,
                fontweight='bold')

        # ── Axes ─────────────────────────────────────────────────
        # ymax: tall enough for rotated labels above bars
        ymax   = int(np.ceil(vmax * 1.65))
        step   = max(1, ymax // 5)
        yticks = list(range(0, ymax + 1, step))

        ax.set_xlim(-0.8, len(vals) - 0.2)
        ax.set_ylim(0, ymax)
        ax.set_xticks([])
        prism_axes(ax, ymax=max(yticks), yticks=yticks)

        ax.set_xlabel('Residues (ranked by contacts)',
                      fontsize=8.5, labelpad=6)
        ax.set_ylabel('Number of contacts', fontsize=8.5)
        ax.set_title(f'{panel}   {name}\n({label})',
                     loc='left', fontsize=9,
                     pad=8, linespacing=1.55)
        ax.text(-0.14, 1.07, panel,
                transform=ax.transAxes,
                fontsize=12, fontweight='bold',
                color=P['text'])

    p = os.path.join(OUT, 'Fig4_ContactDistributions.png')
    fig.savefig(p)
    plt.close()
    print(f'  Saved: {p}')


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print('Phase 0 validation \u2014 Prism-style figures\n')
    fig1()
    fig2()
    fig3()
    fig4()
    print(f'\nAll figures saved to:\n  {OUT}')
