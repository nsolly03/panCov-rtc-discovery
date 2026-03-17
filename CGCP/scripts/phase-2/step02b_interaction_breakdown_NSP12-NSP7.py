#!/usr/bin/env python3
"""
CGCP Phase 2 Step 2b - Interaction Breakdown: NSP12-NSP7
3-panel figure:
  A: Stacked bars — SB / Hbond / Hydrophobic per residue
  B: Contact partner map — NSP12 residue contacts which NSP7 residue
  C: Residue chemistry classification

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step02b_interaction_breakdown_NSP12-NSP7.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')
from Bio import PDB

plt.rcParams.update({
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':          9,
    'axes.linewidth':     0.75,
    'axes.spines.top':    False,
    'axes.spines.right':  False,
    'axes.facecolor':     'white',
    'axes.labelsize':     9,
    'axes.labelpad':      6,
    'axes.titlesize':     9,
    'axes.titlepad':      8,
    'xtick.direction':    'out',
    'ytick.direction':    'out',
    'xtick.major.size':   4,
    'ytick.major.size':   4,
    'xtick.major.width':  0.75,
    'ytick.major.width':  0.75,
    'xtick.major.pad':    4,
    'ytick.major.pad':    4,
    'xtick.labelsize':    8,
    'ytick.labelsize':    8,
    'axes.grid':          False,
    'legend.fontsize':    8,
    'legend.frameon':     False,
    'figure.facecolor':   'white',
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.facecolor':  'white',
    'savefig.pad_inches': 0.15,
    'pdf.fonttype':       42,
})

# Prism palette
P = {
    'sb':      '#7B2D8B',   # purple  - salt bridge
    'sb_f':    '#C994C7',
    'hb':      '#2166AC',   # blue    - H-bond
    'hb_f':    '#92C5DE',
    'hy':      '#E6550D',   # orange  - hydrophobic
    'hy_f':    '#FDAE6B',
    'anchor':  '#D7191C',   # red     - PHE440
    'nsp7':    '#4DAC26',   # green   - NSP7
    'gray':    '#636363',
    'text':    '#1A1A1A',
}

# Residue chemistry classification
CHEM = {
    'ALA': 'hydrophobic', 'VAL': 'hydrophobic', 'ILE': 'hydrophobic',
    'LEU': 'hydrophobic', 'MET': 'hydrophobic', 'PRO': 'hydrophobic',
    'PHE': 'aromatic',    'TYR': 'aromatic',    'TRP': 'aromatic',
    'GLY': 'hydrophobic',
    'SER': 'polar',       'THR': 'polar',       'CYS': 'polar',
    'ASN': 'polar',       'GLN': 'polar',
    'LYS': 'charged_pos', 'ARG': 'charged_pos', 'HIS': 'charged_pos',
    'ASP': 'charged_neg', 'GLU': 'charged_neg',
}
CHEM_COLOR = {
    'aromatic':    '#D7191C',
    'hydrophobic': '#E6550D',
    'polar':       '#2166AC',
    'charged_pos': '#7B2D8B',
    'charged_neg': '#4DAC26',
}
CHEM_LABEL = {
    'aromatic':    'Aromatic (pi-stack)',
    'hydrophobic': 'Hydrophobic',
    'polar':       'Polar (H-bond)',
    'charged_pos': 'Charged+ (SB donor)',
    'charged_neg': 'Charged- (SB acceptor)',
}

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
VAL_DIR = os.path.join(BASE, '02-validation/NSP12-NSP7')
PDB_FILE = os.path.join(BASE,
    '03-virtual-screening/NSP12-NSP7_3/receptor_NSP12-NSP7_3.pdb')
OUT_DIR = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP7/step-02-contacts')
os.makedirs(OUT_DIR, exist_ok=True)


# ── Load CSV data ────────────────────────────────────────────
def load_csv():
    rows = []
    with open(os.path.join(VAL_DIR,
              'composite_ranking_NSP12-NSP7_3.csv')) as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


# ── Compute contact partners from PDB ───────────────────────
def compute_partners(cutoff=5.0):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', PDB_FILE)[0]
    chain_a = structure['A']
    chain_c = structure['C']

    res_a = [r for r in chain_a if r.get_id()[0] == ' ']
    res_c = [r for r in chain_c if r.get_id()[0] == ' ']

    # Build residue name lookup
    aa_lookup = {}
    for r in res_a:
        aa_lookup[('A', r.get_id()[1])] = r.get_resname()
    for r in res_c:
        aa_lookup[('C', r.get_id()[1])] = r.get_resname()

    partners = {}   # {(chain, pos): [(partner_chain, partner_pos, partner_aa, min_dist)]}

    for ra in res_a:
        pos_a = ra.get_id()[1]
        key   = ('A', pos_a)
        partners[key] = []
        for rc in res_c:
            pos_c = rc.get_id()[1]
            min_d = float('inf')
            for aa in ra.get_atoms():
                for ac in rc.get_atoms():
                    d = np.linalg.norm(aa.coord - ac.coord)
                    if d < min_d:
                        min_d = d
            if min_d <= cutoff:
                partners[key].append({
                    'partner_chain': 'C',
                    'partner_pos':   pos_c,
                    'partner_aa':    rc.get_resname(),
                    'min_dist':      round(float(min_d), 2),
                })
        # Sort by distance
        partners[key].sort(key=lambda x: x['min_dist'])

    return partners, aa_lookup


# ── PANEL A: Stacked interaction bars ───────────────────────
def panel_a(ax, rows):
    # NSP12 top 10 by composite
    nsp12 = [r for r in rows if r['chain'] == 'NSP12']
    nsp12.sort(key=lambda x: -float(x['composite']))
    nsp12 = nsp12[:10]

    labels = []
    sb_vals, hb_vals, hy_vals = [], [], []

    for r in nsp12:
        aa  = r['residue_aa'].split('-')[-1]
        aa  = ''.join([c for c in aa if not c.isdigit()])
        aa  = aa.replace('NSP12','').replace('NSP7','').strip()[:3]
        pos = r['position']
        labels.append(f"{aa}{pos}")
        sb_vals.append(int(float(r['sb_loss'])))
        hb_vals.append(int(float(r['hb_loss'])))
        hy_vals.append(int(float(r['hy_loss'])))

    x = np.arange(len(nsp12))
    w = 0.42

    # Stacked bars
    b1 = ax.bar(x, hy_vals, width=w,
                facecolor=P['hy_f'], edgecolor=P['hy'],
                linewidth=0.9, zorder=3, label='Hydrophobic')
    b2 = ax.bar(x, hb_vals, width=w,
                bottom=hy_vals,
                facecolor=P['hb_f'], edgecolor=P['hb'],
                linewidth=0.9, zorder=3, label='H-bond')
    b3 = ax.bar(x, sb_vals, width=w,
                bottom=[h+b for h,b in zip(hy_vals, hb_vals)],
                facecolor=P['sb_f'], edgecolor=P['sb'],
                linewidth=0.9, zorder=3, label='Salt bridge')

    # Total label on top
    for xi, hy, hb, sb in zip(x, hy_vals, hb_vals, sb_vals):
        total = hy + hb + sb
        if total > 0:
            ax.text(xi, total + 0.15, str(total),
                    ha='center', va='bottom',
                    fontsize=7, color=P['text'],
                    fontweight='bold', clip_on=False)

    # Anchor marker
    ax.text(0, -1.2, 'anchor', ha='center',
            fontsize=6.5, color=P['anchor'],
            fontstyle='italic', clip_on=False)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.5,
                       rotation=90, ha='center', va='top')
    ax.tick_params(axis='x', pad=2, length=3, width=0.75)
    ax.set_xlim(-0.6, len(nsp12) - 0.4)
    ax.set_ylabel('Number of contacts', fontsize=9)

    ymax = max([h+b+s for h,b,s in
                zip(hy_vals,hb_vals,sb_vals)]) + 2
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['left'].set_bounds(0, ymax)
    ax.set_ylim(-1.5, ymax * 1.1)
    ax.set_yticks(range(0, int(ymax)+2, 2))

    ax.set_title(
        'A   Interaction type breakdown — NSP12 residues\n'
        '(hydrophobic / H-bond / salt bridge contacts at 5.0 A)',
        loc='left', fontsize=8.5, pad=6, linespacing=1.5)

    return ax


# ── PANEL B: Contact partner map ────────────────────────────
def panel_b(ax, rows, partners, aa_lookup):
    # Top 8 NSP12 hotspots
    nsp12 = [r for r in rows if r['chain'] == 'NSP12']
    nsp12.sort(key=lambda x: -float(x['composite']))
    nsp12 = nsp12[:8]

    # Collect all NSP7 partners
    nsp7_partners = set()
    for r in nsp12:
        pos = int(r['position'])
        key = ('A', pos)
        for p in partners.get(key, []):
            nsp7_partners.add(
                (p['partner_pos'],
                 p['partner_aa'][:3]))
    nsp7_list = sorted(nsp7_partners, key=lambda x: x[0])

    # Build matrix
    nsp12_labels = []
    for r in nsp12:
        aa  = r['residue_aa'].split('-')[-1]
        aa  = ''.join([c for c in aa if not c.isdigit()])
        aa  = aa.replace('NSP12','').replace('NSP7','').strip()[:3]
        nsp12_labels.append(f"{aa}{r['position']}")

    nsp7_labels = [f"{aa}{pos}" for pos, aa in nsp7_list]

    matrix = np.zeros((len(nsp12), len(nsp7_list)))
    for i, r in enumerate(nsp12):
        pos = int(r['position'])
        key = ('A', pos)
        for p in partners.get(key, []):
            for j, (p7pos, _) in enumerate(nsp7_list):
                if p['partner_pos'] == p7pos:
                    matrix[i, j] = p['min_dist']

    # Plot as bubble chart
    for i in range(len(nsp12)):
        for j in range(len(nsp7_list)):
            d = matrix[i, j]
            if d > 0:
                size = max(20, 300 * (1 - (d - 2.5) / 3.0))
                is_anchor = (nsp12[i]['position'] == '440')
                color = P['anchor'] if is_anchor else P['hb']
                ax.scatter(j, i, s=size,
                           color=color, alpha=0.7,
                           edgecolors='white',
                           linewidths=0.5, zorder=3)
                ax.text(j, i, f'{d:.1f}',
                        ha='center', va='center',
                        fontsize=5.5, color='white',
                        fontweight='bold')

    ax.set_xticks(range(len(nsp7_labels)))
    ax.set_xticklabels(nsp7_labels, fontsize=7.5,
                       rotation=90, ha='center', va='top')
    ax.set_yticks(range(len(nsp12_labels)))
    ax.set_yticklabels(nsp12_labels, fontsize=8)
    ax.set_xlim(-0.6, len(nsp7_labels) - 0.4)
    ax.set_ylim(-0.6, len(nsp12) - 0.4)
    ax.invert_yaxis()

    ax.set_xlabel('NSP7 partner residue (chain C)', fontsize=9)
    ax.set_ylabel('NSP12 residue (chain A)', fontsize=9)

    # Light grid lines
    for i in range(len(nsp12)):
        ax.axhline(i, color='#EEEEEE', linewidth=0.5, zorder=1)
    for j in range(len(nsp7_labels)):
        ax.axvline(j, color='#EEEEEE', linewidth=0.5, zorder=1)

    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))

    ax.set_title(
        'B   Contact partner map — NSP12 vs NSP7\n'
        '(bubble size = proximity; label = min dist A)',
        loc='left', fontsize=8.5, pad=6, linespacing=1.5)

    return ax


# ── PANEL C: Residue chemistry ───────────────────────────────
def panel_c(ax, rows, aa_lookup):
    all_res = [r for r in rows
               if r['chain'] in ('NSP12', 'NSP7')]
    all_res.sort(key=lambda x: -float(x['composite']))
    all_res = all_res[:16]

    labels, colors, vals = [], [], []
    for r in all_res:
        aa  = r['residue_aa'].split('-')[-1]
        aa  = ''.join([c for c in aa if not c.isdigit()])
        aa  = aa.replace('NSP12','').replace('NSP7','').strip()[:3]
        pos = r['position']
        chain = r['chain']
        labels.append(f"{chain[:4]}\n{aa}{pos}")
        chem  = CHEM.get(aa.upper(), 'hydrophobic')
        colors.append(CHEM_COLOR[chem])
        vals.append(float(r['composite']))

    x = np.arange(len(all_res))
    w = 0.42
    for xi, val, col, lbl in zip(x, vals, colors, labels):
        ax.bar(xi, val, width=w,
               facecolor=col + '55',
               edgecolor=col,
               linewidth=0.9, zorder=3,
               clip_on=False)
        ax.scatter(xi, val, color=col,
                   s=16, zorder=5, clip_on=False)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6.5,
                       rotation=90, ha='center', va='top')
    ax.tick_params(axis='x', pad=2, length=3, width=0.75)
    ax.set_xlim(-0.6, len(all_res) - 0.4)
    ax.set_ylabel('CGCP composite score', fontsize=9)

    ymax = 1.25
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['left'].set_bounds(0, 1.0)
    ax.set_ylim(-0.08, ymax)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

    ax.set_title(
        'C   Residue chemistry — all interface residues\n'
        '(color = biophysical class; NSP12 + NSP7 combined)',
        loc='left', fontsize=8.5, pad=6, linespacing=1.5)

    # Chemistry legend
    handles = [mpatches.Patch(
                   facecolor=CHEM_COLOR[k]+'55',
                   edgecolor=CHEM_COLOR[k],
                   linewidth=0.9,
                   label=CHEM_LABEL[k])
               for k in CHEM_COLOR]
    ax.legend(handles=handles, loc='upper right',
              fontsize=7, frameon=True,
              facecolor='white', edgecolor='#CCCCCC',
              framealpha=0.9)

    return ax


# ── MAIN FIGURE ──────────────────────────────────────────────
def make_figure(rows, partners, aa_lookup):
    fig = plt.figure(figsize=(16, 13))
    gs  = fig.add_gridspec(2, 2,
                           hspace=0.65, wspace=0.45,
                           left=0.08, right=0.97,
                           top=0.94, bottom=0.12)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])   # full width bottom

    panel_a(ax_a, rows)
    panel_b(ax_b, rows, partners, aa_lookup)
    panel_c(ax_c, rows, aa_lookup)

    figpath = os.path.join(OUT_DIR,
        'Fig_Step02b_InteractionBreakdown.png')
    fig.savefig(figpath)
    plt.close()
    print(f"  Figure: {figpath}")


# ── MAIN ─────────────────────────────────────────────────────
if __name__ == '__main__':
    print('CGCP Phase 2 Step 2b - Interaction Breakdown')
    print('Computing contact partners from PDB...')

    rows = load_csv()
    partners, aa_lookup = compute_partners(cutoff=5.0)

    print(f"  Partners computed for {len(partners)} NSP12 residues")
    print('Building figure...')
    make_figure(rows, partners, aa_lookup)

    # Print partner summary
    print('\nContact partner summary (top NSP12 residues):')
    nsp12 = [r for r in rows if r['chain'] == 'NSP12']
    nsp12.sort(key=lambda x: -float(x['composite']))
    for r in nsp12[:6]:
        pos = int(r['position'])
        aa  = r['residue_aa'].split('-')[-1]
        aa  = ''.join([c for c in aa if not c.isdigit()])
        aa  = aa.replace('NSP12','').strip()[:3]
        key = ('A', pos)
        pts = partners.get(key, [])
        pt_str = ', '.join(
            [f"{p['partner_aa'][:3]}{p['partner_pos']}({p['min_dist']}A)"
             for p in pts[:4]])
        print(f"  {aa}{pos:<6} contacts: {pt_str}")

    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Phase 2 Step 3 - feature classification")
