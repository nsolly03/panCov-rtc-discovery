#!/usr/bin/env python3
"""
CGCP Phase 2 Step 2 — Raw Contact Mapping: NSP12-NSP7
Reads existing validation outputs — no recomputation.
Reformats into CGCP Step 2 standard TSV + JSON + Prism figure.

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step02_contact_mapping_NSP12-NSP7.py
"""

import os
import json
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

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
    'xtick.major.pad':    5,
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
    'savefig.pad_inches': 0.12,
    'pdf.fonttype':       42,
})

P = {
    'nsp12':   '#2166AC',
    'nsp12_f': '#92C5DE',
    'nsp7':    '#4DAC26',
    'nsp7_f':  '#A6D96A',
    'anchor':  '#D7191C',
    'anchor_f':'#FDAE61',
    'gray':    '#636363',
    'text':    '#1A1A1A',
}

# ── Paths ────────────────────────────────────────────────────
BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
VAL_DIR = os.path.join(BASE, '02-validation/NSP12-NSP7')
OUT_DIR = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP7/step-02-contacts')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Read existing data ───────────────────────────────────────
def load_data():
    # Interface analysis JSON
    with open(os.path.join(VAL_DIR, 'interface_analysis_3.json')) as f:
        iface = json.load(f)

    # Composite ranking CSV
    rows = []
    with open(os.path.join(VAL_DIR,
              'composite_ranking_NSP12-NSP7_3.csv')) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    # Conservation CSVs
    cons_nsp12 = {}
    with open(os.path.join(VAL_DIR, 'conservation_NSP12.csv')) as f:
        reader = csv.DictReader(f)
        for row in reader:
            cons_nsp12[int(row['position'])] = {
                'aa':   row['aa_SARS2'],
                'cons': float(row['conservation']),
            }

    cons_nsp7 = {}
    with open(os.path.join(VAL_DIR, 'conservation_NSP7.csv')) as f:
        reader = csv.DictReader(f)
        for row in reader:
            cons_nsp7[int(row['position'])] = {
                'aa':   row['aa_SARS2'],
                'cons': float(row['conservation']),
            }

    return iface, rows, cons_nsp12, cons_nsp7

# ── Build contact table ──────────────────────────────────────
def build_contact_table(iface, rows, cons_nsp12, cons_nsp7):
    contacts = []

    for row in rows:
        chain   = row['chain']
        pos     = int(row['position'])
        aa      = row['residue_aa'].split('-')[-1].replace(
                      chain, '').strip()
        bsa     = float(row['bsa'])
        ct      = float(row['contact_score'])
        comp    = float(row['composite'])
        is_sb   = row['primary_sb'] == 'True'
        hcore   = row['hcore'] == 'True'

        # Conservation from CSV
        if chain == 'NSP12':
            cons = cons_nsp12.get(pos, {}).get('cons', 0.0)
        else:
            cons = cons_nsp7.get(pos, {}).get('cons', 0.0)

        # Is it a known hotspot?
        hs_list = (iface['hotspots_NSP12']
                   if chain == 'NSP12'
                   else iface['hotspots_NSP7'])
        is_hotspot = pos in hs_list

        contacts.append({
            'residue':    f'{chain}-{aa}{pos}',
            'chain':      chain,
            'position':   pos,
            'aa':         aa[:3] if len(aa) >= 3 else aa,
            'bsa':        round(bsa, 2),
            'contacts':   int(ct),
            'conservation': round(cons, 3),
            'composite':  round(comp, 4),
            'salt_bridge': is_sb,
            'hcore':      hcore,
            'is_hotspot': is_hotspot,
        })

    # Sort by composite score descending
    contacts.sort(key=lambda x: -x['composite'])
    return contacts

# ── Save TSV ─────────────────────────────────────────────────
def save_tsv(contacts):
    path = os.path.join(OUT_DIR, 'contact_map_NSP12-NSP7.tsv')
    fields = ['residue', 'chain', 'position', 'aa',
              'bsa', 'contacts', 'conservation',
              'composite', 'salt_bridge', 'hcore', 'is_hotspot']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(contacts)
    print(f"  TSV: {path}")

# ── Save JSON ────────────────────────────────────────────────
def save_json(contacts, iface):
    nsp12 = [c for c in contacts if c['chain'] == 'NSP12']
    nsp7  = [c for c in contacts if c['chain'] == 'NSP7']
    out = {
        'interface':      'NSP12-NSP7',
        'date':           '2026-03-18',
        'method':         'BSA + contact score from BioPython '
                          '(02-validation pipeline)',
        'cutoff_A':       5.0,
        'n_total':        len(contacts),
        'n_nsp12':        len(nsp12),
        'n_nsp7':         len(nsp7),
        'top_anchor':     'PHE440',
        'anchor_composite': next(
            c['composite'] for c in contacts
            if c['position'] == 440
        ),
        'salt_bridges':   iface['salt_bridges_all'],
        'hotspots_nsp12': iface['hotspots_NSP12'],
        'hotspots_nsp7':  iface['hotspots_NSP7'],
        'contacts':       contacts,
    }
    path = os.path.join(OUT_DIR, 'contact_map_NSP12-NSP7.json')
    with open(path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"  JSON: {path}")
    return out

# ── Prism figure — 2 panels ──────────────────────────────────
def prism_axes(ax, ymax=None, yticks=None):
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    if ymax is not None:
        ax.spines['left'].set_bounds(0, ymax)
        ax.set_ylim(-ymax * 0.02, ymax * 1.05)
    if yticks is not None:
        ax.set_yticks(yticks)

def make_figure(contacts):
    # Top 10 NSP12 and top 8 NSP7 by composite
    nsp12 = [c for c in contacts
             if c['chain'] == 'NSP12'][:10]
    nsp7  = [c for c in contacts
             if c['chain'] == 'NSP7'][:8]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.0))
    fig.subplots_adjust(wspace=0.48, bottom=0.20,
                        top=0.88, left=0.09, right=0.97)

    for ax, group, color_f, color_e, title_lbl, panel in [
        (axes[0], nsp12, P['nsp12_f'], P['nsp12'],
         'NSP12 (chain A)', 'A'),
        (axes[1], nsp7,  P['nsp7_f'],  P['nsp7'],
         'NSP7 (chain C)',  'B'),
    ]:
        labels = [f"{c['aa']}{c['position']}" for c in group]
        vals   = [c['composite'] for c in group]
        x      = np.arange(len(group))
        w      = 0.50

        for i, (xi, val, c) in enumerate(zip(x, vals, group)):
            # Anchor gets red
            fc = P['anchor_f'] if c['position'] == 440 \
                 and c['chain'] == 'NSP12' else color_f
            ec = P['anchor']   if c['position'] == 440 \
                 and c['chain'] == 'NSP12' else color_e
            ax.bar(xi, val, width=w,
                   facecolor=fc, edgecolor=ec,
                   linewidth=0.9, zorder=3,
                   clip_on=False)
            # Individual data point dot
            ax.scatter(xi, val, color=ec,
                       s=18, zorder=5,
                       clip_on=False)
            ax.text(xi, val + 0.01,
                    f'{val:.3f}',
                    ha='center', va='bottom',
                    fontsize=6.5, color=ec,
                    clip_on=False)

            # Salt bridge marker
            if c['salt_bridge']:
                ax.text(xi, -0.04, '★',
                        ha='center', fontsize=7,
                        color='purple', clip_on=False)

        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7.5,
                           rotation=90, ha='center')
        ax.set_ylabel('CGCP composite score', fontsize=9)
        ax.set_xlim(-0.55, len(group) - 0.45)

        ymax = max(vals) * 1.3
        step = 0.2
        yticks = [round(i * step, 1)
                  for i in range(int(ymax / step) + 2)]
        prism_axes(ax, ymax=max(yticks), yticks=yticks)

        ax.set_title(
            f'{panel}   Interface contact ranking — {title_lbl}\n'
            f'(BSA + contact score + conservation, 5.0 Å cutoff)',
            loc='left', fontsize=8.5, pad=6,
            linespacing=1.5,
        )

    # Legend
    handles = [
        mpatches.Patch(facecolor=P['anchor_f'],
                       edgecolor=P['anchor'],
                       linewidth=0.9,
                       label='PHE440 — primary anchor'),
        mpatches.Patch(facecolor=P['nsp12_f'],
                       edgecolor=P['nsp12'],
                       linewidth=0.9,
                       label='NSP12 interface residues'),
        mpatches.Patch(facecolor=P['nsp7_f'],
                       edgecolor=P['nsp7'],
                       linewidth=0.9,
                       label='NSP7 interface residues'),
        plt.Line2D([0], [0], marker='*', color='w',
                   markerfacecolor='purple',
                   markersize=9,
                   label='Salt bridge partner'),
    ]
    fig.legend(handles=handles, loc='lower center',
               ncol=4, fontsize=7.5, frameon=False,
               bbox_to_anchor=(0.53, 0.00))

    path = os.path.join(OUT_DIR,
                        'Fig_Step02_ContactMap_Prism.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")

# ── Print summary ────────────────────────────────────────────
def print_summary(contacts, out):
    print('\n' + '=' * 60)
    print('STEP 2 CONTACT MAP SUMMARY — NSP12-NSP7')
    print('=' * 60)
    print(f"  Total interface residues: {out['n_total']}")
    print(f"  NSP12 residues: {out['n_nsp12']}")
    print(f"  NSP7  residues: {out['n_nsp7']}")
    print(f"\n  Top 5 by composite score:")
    print(f"  {'Residue':<20} {'BSA':>6} "
          f"{'Contacts':>9} {'Cons':>6} {'Composite':>10}")
    print(f"  {'-'*55}")
    for c in contacts[:5]:
        sb = ' ★SB' if c['salt_bridge'] else ''
        print(f"  {c['residue']:<20} "
              f"{c['bsa']:>6.1f} "
              f"{c['contacts']:>9} "
              f"{c['conservation']:>6.3f} "
              f"{c['composite']:>10.4f}{sb}")
    print(f"\n  Salt bridges:")
    for sb in out['salt_bridges']:
        print(f"    {sb}")
    print('=' * 60)

# ── Main ─────────────────────────────────────────────────────
if __name__ == '__main__':
    print('CGCP Phase 2 Step 2 — Contact Mapping: NSP12-NSP7')
    print('Reading from existing validation outputs...\n')

    iface, rows, cons_nsp12, cons_nsp7 = load_data()
    contacts = build_contact_table(
        iface, rows, cons_nsp12, cons_nsp7)

    save_tsv(contacts)
    out = save_json(contacts, iface)
    make_figure(contacts)
    print_summary(contacts, out)

    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Phase 2 Step 3 — feature classification")
