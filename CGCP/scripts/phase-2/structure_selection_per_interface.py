#!/usr/bin/env python3
"""
CGCP Structure Selection — One publication figure per interface.
Shows BOTH resolution and R-free; annotates the actual selection
criterion used (resolution, R-free, interface completeness, AF3).

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/structure_selection_per_interface.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines  as mlines
import numpy as np
import os

# ── Publication theme (theme_classic equivalent) ─────────────
RC = {
    'font.family':       'DejaVu Sans',
    'font.size':         12,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.grid':         False,
    'axes.linewidth':    1.3,
    'xtick.direction':   'out',
    'ytick.direction':   'out',
    'xtick.major.width': 1.3,
    'ytick.major.width': 1.3,
    'xtick.major.size':  6,
    'ytick.major.size':  6,
    'figure.facecolor':  'white',
    'axes.facecolor':    'white',
    'pdf.fonttype':      42,
    'ps.fonttype':       42,
}
plt.rcParams.update(RC)

# ── Colour palette ────────────────────────────────────────────
C = {
    'sel_e':  '#1A5276', 'sel_f':  '#D6EAF8',   # selected
    'alt_e':  '#6C7A89', 'alt_f':  '#ECF0F1',   # alternative
    'fail_e': '#922B21', 'fail_f': '#FADBD8',   # failed criterion
    'af3_e':  '#B7770D', 'af3_f':  '#FEF9E7',   # AlphaFold3
    'res':    '#1A5276',                         # resolution bar label
    'rfree':  '#C0392B',                         # R-free
    'note':   '#2C3E50',
}

# ── Interface data ────────────────────────────────────────────
# criterion: 'resolution' | 'rfree' | 'completeness' | 'af3'
INTERFACES = {

  'NSP12-NSP7': {
    'title':     'NSP12–NSP7 Interface',
    'selected':  '7C2K',
    'criterion': 'completeness',
    'crit_note': 'Selected for interface completeness:\n7C2K resolves NSP12+NSP7+NSP8 ternary complex\n(7BV2/6NUR lack NSP7 chain at interface)',
    'structs': [
      {'id':'7BV2', 'res':2.50,'method':'Cryo-EM','rfree':None,
       'selected':False,'note':'No NSP7 at interface'},
      {'id':'7C2K', 'res':2.93,'method':'Cryo-EM','rfree':None,
       'selected':True, 'note':'NSP12+NSP7+NSP8 ternary'},
      {'id':'6NUR', 'res':3.10,'method':'Cryo-EM','rfree':None,
       'selected':False,'note':'Lower resolution'},
    ]
  },

  'NSP12-NSP8': {
    'title':     'NSP12–NSP8 Interface',
    'selected':  '7C2K',
    'criterion': 'completeness',
    'crit_note': 'Selected for interface completeness:\n7C2K contains two NSP8 chains (B+D)\n(7BV2/6NUR have single NSP8 chain only)',
    'structs': [
      {'id':'7BV2', 'res':2.50,'method':'Cryo-EM','rfree':None,
       'selected':False,'note':'Single NSP8 chain only'},
      {'id':'7C2K', 'res':2.93,'method':'Cryo-EM','rfree':None,
       'selected':True, 'note':'Two NSP8 chains (B+D)'},
      {'id':'6NUR', 'res':3.10,'method':'Cryo-EM','rfree':None,
       'selected':False,'note':'Lower resolution'},
    ]
  },

  'NSP9-NSP12': {
    'title':     'NSP9–NSP12 Interface',
    'selected':  '8SQK',
    'criterion': 'completeness',
    'crit_note': 'Only structure capturing the NiRAN–NSP9\nco-crystal interface (unique structure)',
    'structs': [
      {'id':'8SQK', 'res':3.01,'method':'Cryo-EM','rfree':None,
       'selected':True,'note':'Unique NiRAN–NSP9 complex'},
    ]
  },

  'NSP10-NSP16': {
    'title':     'NSP10–NSP16 Interface',
    'selected':  '6W4H',
    'criterion': 'rfree',
    'crit_note': 'Selected by R-free:\n6W4H R-free=0.070 (exceptionally low)\nDespite 6WKQ/6WVN at higher resolution',
    'structs': [
      {'id':'6W4H', 'res':1.80,'method':'X-ray','rfree':0.070,
       'selected':True, 'note':'Best R-free (0.070)'},
      {'id':'6WKQ', 'res':1.98,'method':'X-ray','rfree':0.097,
       'selected':False,'note':'Higher R-free'},
      {'id':'6WVN', 'res':2.00,'method':'X-ray','rfree':0.097,
       'selected':False,'note':'Higher R-free'},
    ]
  },

  'NSP7-NSP8': {
    'title':     'NSP7–NSP8 Interface',
    'selected':  'AF3 ModeB',
    'criterion': 'af3',
    'crit_note': 'AF3 ModeB selected — crystal structures failed\ngeometry criterion: PHE92 (anchor) is 35 A\nfrom NSP7 partner in all PDB structures (ModeA).\nAF3 ModeB: PHE92 = 6.78 A (druggable).',
    'structs': [
      {'id':'7BV2\n(ModeA)', 'res':2.50,'method':'Cryo-EM','rfree':None,
       'selected':False,'failed':True,'note':'PHE92 = 35 A from partner'},
      {'id':'6NUR\n(ModeA)', 'res':3.10,'method':'Cryo-EM','rfree':None,
       'selected':False,'failed':True,'note':'PHE92 = 35 A from partner'},
      {'id':'AF3\nModeB',    'res':None,'method':'AlphaFold3','rfree':None,
       'selected':True,'note':'PHE92 = 6.78 A (druggable)'},
    ]
  },

  'NSP10-NSP14': {
    'title':     'NSP10–NSP14 Interface',
    'selected':  '7DIY',
    'criterion': 'completeness',
    'crit_note': 'Only structure available for the\nNSP10–NSP14 exonuclease proofreading complex',
    'structs': [
      {'id':'7DIY', 'res':2.69,'method':'X-ray','rfree':0.264,
       'selected':True,'note':'Unique exonuclease complex'},
    ]
  },

  'NSP13-Helicase': {
    'title':     'NSP13 Helicase Homodimer Interface',
    'selected':  '7NIO',
    'criterion': 'resolution',
    'crit_note': 'Selected by resolution:\n7NIO (X-ray, 2.20 A) vs 7CXM (Cryo-EM, 2.90 A)\nX-ray preferred for higher accuracy at interface',
    'structs': [
      {'id':'7NIO', 'res':2.20,'method':'X-ray', 'rfree':0.287,
       'selected':True, 'note':'Best resolution, X-ray'},
      {'id':'7CXM', 'res':2.90,'method':'Cryo-EM','rfree':None,
       'selected':False,'note':'Lower resolution'},
    ]
  },

  'NSP12-NSP13': {
    'title':     'NSP12–NSP13 Interface',
    'selected':  '7RDY',
    'criterion': 'completeness',
    'crit_note': 'Selected for interface completeness:\n7RDY explicitly captures polymerase–helicase\njunction residues (NSP12 C-tail + NSP13 N-term)\n6XEZ lower resolution; 7CXM lacks junction',
    'structs': [
      {'id':'7CXM', 'res':2.90,'method':'Cryo-EM','rfree':None,
       'selected':False,'note':'Lacks junction residues'},
      {'id':'7RDY', 'res':3.10,'method':'Cryo-EM','rfree':None,
       'selected':True, 'note':'Captures interface residues'},
      {'id':'6XEZ', 'res':3.50,'method':'Cryo-EM','rfree':None,
       'selected':False,'note':'Lowest resolution'},
    ]
  },

  'NSP15': {
    'title':     'NSP15 Endoribonuclease Homodimer Interface',
    'selected':  '9HH5',
    'criterion': 'rfree',
    'crit_note': 'Selected by R-free:\n9HH5 R-free=0.023 (exceptionally low vs 0.178-0.219\nfor higher-resolution alternatives)\nLow R-free indicates superior model quality',
    'structs': [
      {'id':'6WLC', 'res':1.82,'method':'X-ray', 'rfree':0.195,'selected':False,'note':'High R-free'},
      {'id':'6WXC', 'res':1.85,'method':'X-ray', 'rfree':0.194,'selected':False,'note':'High R-free'},
      {'id':'6W01', 'res':1.90,'method':'X-ray', 'rfree':0.185,'selected':False,'note':'High R-free'},
      {'id':'9HH5', 'res':2.08,'method':'X-ray', 'rfree':0.023,'selected':True, 'note':'Best R-free (0.023)'},
      {'id':'6VWW', 'res':2.20,'method':'X-ray', 'rfree':0.178,'selected':False,'note':'High R-free'},
      {'id':'2H85', 'res':2.60,'method':'X-ray', 'rfree':0.219,'selected':False,'note':'High R-free'},
      {'id':'7TJ2', 'res':3.20,'method':'Cryo-EM','rfree':None, 'selected':False,'note':'Lowest resolution'},
    ]
  },
}

CRIT_COLORS = {
    'resolution':   '#1A5276',
    'rfree':        '#922B21',
    'completeness': '#1E8449',
    'af3':          '#B7770D',
}
CRIT_LABELS = {
    'resolution':   'Selection criterion: Resolution',
    'rfree':        'Selection criterion: R-free',
    'completeness': 'Selection criterion: Interface completeness',
    'af3':          'Selection criterion: AF3 (crystal geometry failure)',
}


def draw_figure(iface_key, data):
    structs   = data['structs']
    n         = len(structs)
    criterion = data['criterion']
    crit_col  = CRIT_COLORS[criterion]

    has_rfree = any(s.get('rfree') is not None for s in structs)
    has_af3   = any(s['method'] == 'AlphaFold3' for s in structs)

    # Figure size: wider for many structures
    fw = max(5.5, n * 1.6 + 2.0)
    fig, ax1 = plt.subplots(figsize=(fw, 5.8))
    fig.subplots_adjust(left=0.13, right=0.87 if has_rfree else 0.95,
                        top=0.82, bottom=0.20)

    xs = np.arange(n)
    bw = 0.52

    # ── Resolution bars ───────────────────────────────────────
    for i, s in enumerate(structs):
        if s.get('failed'):
            ec, fc = C['fail_e'], C['fail_f']
        elif s.get('selected'):
            ec = crit_col
            fc = C['sel_f']
        else:
            ec, fc = C['alt_e'], C['alt_f']

        lw = 2.0 if s.get('selected') else 1.1

        if s['res'] is not None:
            ax1.bar(i, s['res'], width=bw, facecolor=fc,
                    edgecolor=ec, linewidth=lw,
                    zorder=3, clip_on=False)
            ax1.text(i, s['res'] + 0.05, f"{s['res']:.2f} Å",
                     ha='center', va='bottom',
                     fontsize=9,
                     fontweight='bold' if s.get('selected') else 'normal',
                     color=ec, clip_on=False)
        else:
            # AF3 placeholder hatched bar
            ax1.bar(i, 3.8, width=bw,
                    facecolor=C['af3_f'], edgecolor=C['af3_e'],
                    linewidth=2.0, hatch='///', alpha=0.8,
                    zorder=3, clip_on=False)
            ax1.text(i, 1.9, 'Predicted\nmodel (AF3)',
                     ha='center', va='center',
                     fontsize=8.5, color=C['af3_e'],
                     fontweight='bold', clip_on=False)

        # Failed X overlay
        if s.get('failed'):
            h = s['res'] * 0.5 if s['res'] else 1.9
            ax1.text(i, h, 'X', ha='center', va='center',
                     fontsize=20, color=C['fail_e'],
                     fontweight='bold', alpha=0.50, clip_on=False)

        # Selected annotation bracket
        if s.get('selected') and s['method'] != 'AlphaFold3':
            top = s['res'] + 0.60
            ax1.annotate(
                'SELECTED',
                xy=(i, s['res'] + 0.08),
                xytext=(i, top),
                ha='center', va='bottom',
                fontsize=8.5, fontweight='bold',
                color=crit_col,
                arrowprops=dict(arrowstyle='-[,widthB=1.0,lengthB=0.3',
                                color=crit_col, lw=1.2),
                annotation_clip=False)

        if s.get('selected') and s['method'] == 'AlphaFold3':
            ax1.text(i, 4.15, 'SELECTED', ha='center', va='bottom',
                     fontsize=8.5, fontweight='bold',
                     color=crit_col, clip_on=False)

    # ── R-free twin axis ──────────────────────────────────────
    if has_rfree:
        ax2 = ax1.twinx()
        ax2.spines['right'].set_visible(True)
        ax2.spines['right'].set_linewidth(1.1)
        ax2.spines['top'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.set_ylim(0, 0.45)
        ax2.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        ax2.tick_params(axis='y', labelsize=10,
                        colors=C['rfree'],
                        direction='out', length=5, width=1.1)
        ax2.set_ylabel('R-free', fontsize=12,
                       color=C['rfree'], fontweight='bold',
                       labelpad=6)

        for i, s in enumerate(structs):
            if s.get('rfree') is None:
                continue
            sel = s.get('selected')
            ax2.scatter([i], [s['rfree']],
                        marker='D',
                        s=80 if sel else 55,
                        color=C['rfree'],
                        zorder=7,
                        edgecolors='white' if not sel else crit_col,
                        linewidths=1.5 if sel else 0.8)
            offset = 0.22 if i < n - 1 else -0.22
            ax2.text(i + offset, s['rfree'] + 0.010,
                     f"{s['rfree']:.3f}",
                     fontsize=8, color=C['rfree'],
                     va='bottom',
                     fontweight='bold' if sel else 'normal')

    # ── Axis styling ──────────────────────────────────────────
    ax1.set_xticks(xs)
    ax1.set_xticklabels([s['id'] for s in structs],
                        fontsize=10, ha='center')
    ax1.tick_params(axis='x', length=0, pad=7)
    ax1.tick_params(axis='y', direction='out',
                    length=6, width=1.3, labelsize=10)
    ax1.set_ylim(0, 4.9)
    ax1.set_yticks([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
    ax1.set_ylabel('Resolution (Å)', fontsize=12, fontweight='bold')
    ax1.set_xlabel('PDB / Structure ID', fontsize=12, fontweight='bold')
    ax1.spines['left'].set_position(('outward', 5))
    ax1.spines['bottom'].set_position(('outward', 5))
    ax1.set_xlim(-0.6, n - 0.4)

    # "lower is better" arrow
    ax1.annotate('', xy=(-.55, 1.3), xytext=(-.55, 2.8),
                 arrowprops=dict(arrowstyle='->',
                                 color='#AAAAAA', lw=1.0),
                 annotation_clip=False)
    ax1.text(-.55, 2.05, 'Better', fontsize=8,
             color='#AAAAAA', ha='center', va='center',
             rotation=90, clip_on=False)

    # ── Method labels under bars ──────────────────────────────
    for i, s in enumerate(structs):
        col = (C['fail_e'] if s.get('failed')
               else C['af3_e'] if s['method'] == 'AlphaFold3'
               else C['sel_e'] if s.get('selected')
               else C['alt_e'])
        short = {'X-ray':'X-ray','Cryo-EM':'Cryo-EM','AlphaFold3':'AF3'}
        ax1.text(i, -0.38, f"[{short[s['method']]}]",
                 ha='center', va='top',
                 fontsize=8, color=col,
                 fontstyle='italic',
                 clip_on=False)

    # ── Criterion badge ───────────────────────────────────────
    badge = CRIT_LABELS[criterion]
    ax1.text(0.5, 1.04, badge,
             transform=ax1.transAxes,
             ha='center', va='bottom',
             fontsize=9.5, fontweight='bold',
             color='white',
             bbox=dict(boxstyle='round,pad=0.4',
                       facecolor=crit_col,
                       edgecolor='none'))

    # ── Selection rationale note ──────────────────────────────
    ax1.text(1.02 if has_rfree else 0.98,
             0.98, data['crit_note'],
             transform=ax1.transAxes,
             ha='left' if has_rfree else 'right',
             va='top',
             fontsize=8.5, color=C['note'],
             bbox=dict(boxstyle='round,pad=0.45',
                       facecolor='#FDFEFE',
                       edgecolor='#BFC9CA',
                       linewidth=0.8),
             clip_on=False)

    # ── Legend ────────────────────────────────────────────────
    handles = [
        mpatches.Patch(facecolor=C['sel_f'],
                       edgecolor=crit_col, linewidth=2.0,
                       label='Selected structure'),
        mpatches.Patch(facecolor=C['alt_f'],
                       edgecolor=C['alt_e'], linewidth=1.1,
                       label='Candidate (not selected)'),
    ]
    if any(s.get('failed') for s in structs):
        handles.append(mpatches.Patch(
            facecolor=C['fail_f'], edgecolor=C['fail_e'],
            linewidth=1.1, label='Failed interface criterion'))
    if has_af3:
        handles.append(mpatches.Patch(
            facecolor=C['af3_f'], edgecolor=C['af3_e'],
            linewidth=1.5, hatch='///',
            label='AlphaFold3 (predicted model)'))
    if has_rfree:
        handles.append(mlines.Line2D(
            [], [], marker='D', color=C['rfree'],
            linestyle='None', markersize=8,
            label='R-free value'))

    ax1.legend(handles=handles, loc='upper left',
               fontsize=8.5, frameon=True,
               facecolor='white', edgecolor='#BFC9CA',
               framealpha=1.0,
               bbox_to_anchor=(0.01, 0.99),
               borderpad=0.6, handlelength=1.4)

    # ── Title ─────────────────────────────────────────────────
    fig.suptitle(data['title'],
                 fontsize=14, fontweight='bold',
                 x=0.5, y=0.97)

    return fig


# ── Save all figures ──────────────────────────────────────────
OUT = os.path.expanduser(
    '~/projects/rtc-pan-coronavirus/CGCP/02-deep-dive/structure_selection')
os.makedirs(OUT, exist_ok=True)

for iface_key, data in INTERFACES.items():
    fig = draw_figure(iface_key, data)
    for ext in ['png', 'pdf']:
        path = os.path.join(OUT,
            f'Fig_StructureSelection_{iface_key}.{ext}')
        fig.savefig(path, dpi=300,
                    bbox_inches='tight',
                    facecolor='white')
    plt.close(fig)
    print(f"Saved: {iface_key}")

print(f"\nAll figures saved to:\n{OUT}")
