"""
CGCP Prism Style Module
Shared figure aesthetics for all Phase 2 scripts.
Import at top of any script:
    from prism_style import apply_prism, prism_axes, COLORS

GraphPad Prism aesthetic:
  - White background, no gridlines
  - Thick outward ticks
  - Bold axis labels and tick labels
  - Narrow spaced bars
  - Minimal axes (bottom + left only, ending at last tick)
  - Arial/Helvetica font
  - Flat muted colors
  - No chartjunk
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# ── Color palette ─────────────────────────────────────────────
COLORS = {
    # Element colors
    'red':        '#D7191C',
    'red_fill':   '#FDAE61',
    'blue':       '#2166AC',
    'blue_fill':  '#92C5DE',
    'green':      '#1A7D2E',
    'green_fill': '#A8D5B0',
    'orange':     '#E6550D',
    'orange_fill':'#FDD0A2',
    'purple':     '#6A3D9A',
    'purple_fill':'#CAB2D6',
    'gray':       '#636363',
    'gray_fill':  '#CCCCCC',
    'black':      '#1A1A1A',
    'white':      '#FFFFFF',
    # Decision colors
    'anchor':     '#D7191C',
    'anchor_f':   '#FDAE61',
    'include':    '#1A7D2E',
    'include_f':  '#A8D5B0',
    'secondary':  '#2166AC',
    'secondary_f':'#92C5DE',
    'exclude':    '#636363',
    'exclude_f':  '#CCCCCC',
    # Feature colors
    'aromatic':      '#D7191C',
    'aromatic_f':    '#FDAE61',
    'hydrophobic':   '#636363',
    'hydrophobic_f': '#CCCCCC',
    'charged_pos':   '#2166AC',
    'charged_pos_f': '#92C5DE',
    'charged_neg':   '#1A7D2E',
    'charged_neg_f': '#A8D5B0',
    'hbond':         '#8B4513',
    'hbond_f':       '#DEB887',
}

CLUSTER_COLORS = [
    '#D7191C', '#2166AC', '#4DAC26', '#8B4513',
    '#9467BD', '#E6550D', '#1F78B4', '#33A02C',
]
CLUSTER_FILLS = [
    '#FDAE61', '#92C5DE', '#A6D96A', '#DEB887',
    '#C5B0D5', '#FDD0A2', '#A6CEE3', '#B2DF8A',
]


def apply_prism():
    """Apply Prism-style rcParams globally."""
    plt.rcParams.update({
        # Font
        'font.family':          'sans-serif',
        'font.sans-serif':      ['Arial', 'Helvetica',
                                 'DejaVu Sans'],
        'font.size':            10,
        'font.weight':          'normal',

        # Axes
        'axes.linewidth':       1.2,
        'axes.spines.top':      False,
        'axes.spines.right':    False,
        'axes.facecolor':       'white',
        'axes.labelsize':       11,
        'axes.labelweight':     'bold',
        'axes.labelpad':        8,
        'axes.titlepad':        10,
        'axes.titlesize':       10,
        'axes.titleweight':     'normal',

        # Ticks — outward, thick
        'xtick.direction':      'out',
        'ytick.direction':      'out',
        'xtick.major.size':     6,
        'ytick.major.size':     6,
        'xtick.minor.size':     3,
        'ytick.minor.size':     3,
        'xtick.major.width':    1.2,
        'ytick.major.width':    1.2,
        'xtick.minor.width':    0.8,
        'ytick.minor.width':    0.8,
        'xtick.major.pad':      5,
        'ytick.major.pad':      5,
        'xtick.labelsize':      9,
        'ytick.labelsize':      9,

        # Grid — none
        'axes.grid':            False,

        # Legend
        'legend.fontsize':      9,
        'legend.frameon':       True,
        'legend.facecolor':     'white',
        'legend.edgecolor':     '#CCCCCC',
        'legend.framealpha':    0.9,
        'legend.borderpad':     0.5,
        'legend.labelspacing':  0.4,

        # Figure
        'figure.facecolor':     'white',
        'figure.dpi':           100,

        # Save
        'savefig.dpi':          300,
        'savefig.bbox':         'tight',
        'savefig.facecolor':    'white',
        'savefig.pad_inches':   0.2,
        'pdf.fonttype':         42,
        'ps.fonttype':          42,
    })


def prism_axes(ax,
               ymin=0,
               ymax=None,
               yticks=None,
               xmin=None,
               xmax=None,
               spine_offset=6):
    """
    Apply Prism-style spine formatting to an axis.

    Parameters
    ----------
    ax         : matplotlib Axes
    ymin       : float, lower y spine bound (default 0)
    ymax       : float, upper y spine bound
    yticks     : list, explicit y tick positions
    xmin       : float, lower x spine bound
    xmax       : float, upper x spine bound
    spine_offset: int, pixels to offset spines outward
    """
    ax.spines['left'].set_position(
        ('outward', spine_offset))
    ax.spines['bottom'].set_position(
        ('outward', spine_offset))
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)

    if ymax is not None:
        ax.spines['left'].set_bounds(ymin, ymax)
        ax.set_ylim(ymin - ymax * 0.04,
                    ymax * 1.12)
    if yticks is not None:
        ax.set_yticks(yticks)

    if xmin is not None and xmax is not None:
        ax.spines['bottom'].set_bounds(xmin, xmax)


def set_xticklabels_vertical(ax, labels,
                              fontsize=8,
                              fontweight='normal'):
    """
    Set x tick labels vertically to avoid overlap.
    """
    ax.set_xticklabels(
        labels,
        fontsize=fontsize,
        fontweight=fontweight,
        rotation=90,
        ha='center',
        va='top',
    )
    ax.tick_params(axis='x',
                   pad=4,
                   length=6,
                   width=1.2)


def bar_with_label(ax, x, y,
                   facecolor, edgecolor,
                   width=0.45,
                   label_offset=None,
                   fontsize=8,
                   fontweight='bold',
                   fmt='{:.3f}',
                   clip_on=False):
    """
    Draw a bar and add a value label above it.
    """
    ax.bar(x, y, width=width,
           facecolor=facecolor,
           edgecolor=edgecolor,
           linewidth=1.2,
           zorder=3,
           clip_on=clip_on)

    if label_offset is None:
        label_offset = max(y * 0.02, 0.015)

    ax.text(x, y + label_offset,
            fmt.format(y),
            ha='center', va='bottom',
            fontsize=fontsize,
            fontweight=fontweight,
            color=edgecolor,
            clip_on=clip_on)


def add_significance_bracket(ax, x1, x2, y,
                               label='*',
                               color='black',
                               linewidth=1.0,
                               fontsize=10):
    """Add a significance bracket between two bars."""
    h = y * 0.03
    ax.plot([x1, x1, x2, x2],
            [y, y + h, y + h, y],
            color=color,
            linewidth=linewidth)
    ax.text((x1 + x2) / 2, y + h * 1.2,
            label,
            ha='center', va='bottom',
            fontsize=fontsize,
            color=color)


def make_legend(ax, items,
                loc='upper right',
                fontsize=9):
    """
    Make a clean legend.
    items = list of (facecolor, edgecolor, label) tuples
    """
    handles = [
        mpatches.Patch(
            facecolor=fc,
            edgecolor=ec,
            linewidth=1.0,
            label=lbl)
        for fc, ec, lbl in items
    ]
    ax.legend(handles=handles,
              loc=loc,
              fontsize=fontsize,
              frameon=True,
              facecolor='white',
              edgecolor='#CCCCCC',
              framealpha=0.9)
    return handles


def panel_title(ax, letter, title,
                subtitle=None,
                fontsize=10):
    """
    Add a panel title in Prism style.
    letter: 'A', 'B', 'C' etc.
    """
    full = f"{letter}   {title}"
    if subtitle:
        full += f"\n{subtitle}"
    ax.set_title(full,
                 loc='left',
                 fontsize=fontsize,
                 fontweight='normal',
                 pad=8,
                 linespacing=1.5)
