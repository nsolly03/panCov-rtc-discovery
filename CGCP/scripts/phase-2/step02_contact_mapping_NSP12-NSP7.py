#!/usr/bin/env python3
"""
CGCP Phase 2 Step 2 - Raw Contact Mapping: NSP12-NSP7
Reads existing validation outputs - no recomputation.
Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step02_contact_mapping_NSP12-NSP7.py
"""

import os, json, csv
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
    'savefig.pad_inches': 0.12,
    'pdf.fonttype':       42,
})

P = {
    'nsp12':   '#2166AC', 'nsp12_f': '#92C5DE',
    'nsp7':    '#4DAC26', 'nsp7_f':  '#A6D96A',
    'anchor':  '#D7191C', 'anchor_f':'#FDAE61',
    'gray':    '#636363', 'text':    '#1A1A1A',
}

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
VAL_DIR = os.path.join(BASE, '02-validation/NSP12-NSP7')
OUT_DIR = os.path.join(BASE, 'CGCP/02-deep-dive/NSP12-NSP7/step-02-contacts')
os.makedirs(OUT_DIR, exist_ok=True)


def load_data():
    with open(os.path.join(VAL_DIR, 'interface_analysis_3.json')) as f:
        iface = json.load(f)
    rows = []
    with open(os.path.join(VAL_DIR, 'composite_ranking_NSP12-NSP7_3.csv')) as f:
        for row in csv.DictReader(f):
            rows.append(row)
    cons_nsp12 = {}
    with open(os.path.join(VAL_DIR, 'conservation_NSP12.csv')) as f:
        for row in csv.DictReader(f):
            cons_nsp12[int(row['position'])] = float(row['conservation'])
    cons_nsp7 = {}
    with open(os.path.join(VAL_DIR, 'conservation_NSP7.csv')) as f:
        for row in csv.DictReader(f):
            cons_nsp7[int(row['position'])] = float(row['conservation'])
    return iface, rows, cons_nsp12, cons_nsp7


def build_contact_table(iface, rows, cons_nsp12, cons_nsp7):
    contacts = []
    for row in rows:
        chain = row['chain']
        pos   = int(row['position'])
        aa    = row['residue_aa'].split('-')[-1]
        aa    = ''.join([c for c in aa if not c.isdigit()])
        aa    = aa.replace('NSP12','').replace('NSP7','').strip()[:3]
        cons  = cons_nsp12.get(pos, 0.0) if chain == 'NSP12' \
                else cons_nsp7.get(pos, 0.0)
        hs    = iface['hotspots_NSP12'] if chain == 'NSP12' \
                else iface['hotspots_NSP7']
        contacts.append({
            'residue':      f'{chain}-{aa}{pos}',
            'chain':        chain,
            'position':     pos,
            'aa':           aa,
            'bsa':          round(float(row['bsa']), 2),
            'contacts':     int(float(row['contact_score'])),
            'conservation': round(cons, 3),
            'composite':    round(float(row['composite']), 4),
            'salt_bridge':  row['primary_sb'] == 'True',
            'hcore':        row['hcore'] == 'True',
            'is_hotspot':   pos in hs,
        })
    contacts.sort(key=lambda x: -x['composite'])
    return contacts


def save_tsv(contacts):
    fpath = os.path.join(OUT_DIR, 'contact_map_NSP12-NSP7.tsv')
    fields = ['residue','chain','position','aa','bsa','contacts',
              'conservation','composite','salt_bridge','hcore','is_hotspot']
    with open(fpath, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(contacts)
    print(f"  TSV: {fpath}")


def save_json(contacts, iface):
    nsp12 = [c for c in contacts if c['chain'] == 'NSP12']
    nsp7  = [c for c in contacts if c['chain'] == 'NSP7']
    out = {
        'interface':   'NSP12-NSP7',
        'date':        '2026-03-18',
        'n_total':     len(contacts),
        'n_nsp12':     len(nsp12),
        'n_nsp7':      len(nsp7),
        'top_anchor':  'PHE440',
        'salt_bridges': iface['salt_bridges_all'],
        'hotspots_nsp12': iface['hotspots_NSP12'],
        'hotspots_nsp7':  iface['hotspots_NSP7'],
        'contacts':    contacts,
    }
    fpath = os.path.join(OUT_DIR, 'contact_map_NSP12-NSP7.json')
    with open(fpath, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"  JSON: {fpath}")
    return out


def prism_axes(ax, ymax=None, yticks=None):
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    if ymax is not None:
        ax.spines['left'].set_bounds(0, ymax)
        ax.set_ylim(-ymax * 0.02, ymax * 1.05)
    if yticks is not None:
        ax.set_yticks(yticks)


def make_figure(contacts):
    nsp12 = [c for c in contacts if c['chain'] == 'NSP12'][:10]
    nsp7  = [c for c in contacts if c['chain'] == 'NSP7'][:8]

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 6.5))
    fig.subplots_adjust(wspace=0.45, bottom=0.25,
                        top=0.88, left=0.08, right=0.97)

    groups = [
        (axes[0], nsp12, P['nsp12_f'], P['nsp12'], 'NSP12 (chain A)', 'A'),
        (axes[1], nsp7,  P['nsp7_f'],  P['nsp7'],  'NSP7 (chain C)',  'B'),
    ]

    for ax, group, color_f, color_e, title_lbl, panel in groups:
        labels = [c['aa'][:3] + str(c['position']) for c in group]
        vals   = [c['composite'] for c in group]
        x      = np.arange(len(group))
        w      = 0.32

        for xi, val, c in zip(x, vals, group):
            is_anchor = (c['position'] == 440 and c['chain'] == 'NSP12')
            fc = P['anchor_f'] if is_anchor else color_f
            ec = P['anchor']   if is_anchor else color_e
            ax.bar(xi, val, width=w, facecolor=fc, edgecolor=ec,
                   linewidth=0.9, zorder=3, clip_on=False)
            ax.scatter(xi, val, color=ec, s=16, zorder=5, clip_on=False)
            ax.text(xi, val + 0.012, f'{val:.3f}',
                    ha='center', va='bottom', fontsize=6.5,
                    color=ec, clip_on=False)
            if c['salt_bridge']:
                ax.text(xi, -0.055, '*', ha='center',
                        fontsize=9, color='purple', clip_on=False)

        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7.5,
                           rotation=90, ha='center', va='top')
        ax.tick_params(axis='x', pad=2, length=3, width=0.75)
        ax.set_xlim(-0.6, len(group) - 0.4)
        ax.set_ylabel('CGCP composite score', fontsize=9)

        ymax   = round(max(vals) * 1.35 + 0.1, 1)
        step   = 0.2
        yticks = [round(i * step, 1) for i in range(int(ymax/step) + 1)]
        prism_axes(ax, ymax=max(yticks), yticks=yticks)

        panel_title = panel + '   Contact ranking - ' + title_lbl
        ax.set_title(panel_title, loc='left', fontsize=8.5, pad=6)

    handles = [
        mpatches.Patch(facecolor=P['anchor_f'], edgecolor=P['anchor'],
                       linewidth=0.9, label='PHE440 - primary anchor'),
        mpatches.Patch(facecolor=P['nsp12_f'], edgecolor=P['nsp12'],
                       linewidth=0.9, label='NSP12 interface residues'),
        mpatches.Patch(facecolor=P['nsp7_f'], edgecolor=P['nsp7'],
                       linewidth=0.9, label='NSP7 interface residues'),
        plt.Line2D([0],[0], marker='*', color='w',
                   markerfacecolor='purple', markersize=9,
                   label='Salt bridge partner'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=4,
               fontsize=7.5, frameon=False, bbox_to_anchor=(0.53, 0.00))

    figpath = os.path.join(OUT_DIR, 'Fig_Step02_ContactMap_Prism.png')
    fig.savefig(figpath)
    plt.close()
    print(f"  Figure: {figpath}")


def print_summary(contacts, out):
    print('\n' + '='*60)
    print('STEP 2 CONTACT MAP SUMMARY - NSP12-NSP7')
    print('='*60)
    print(f"  Total: {out['n_total']}  NSP12: {out['n_nsp12']}  NSP7: {out['n_nsp7']}")
    print(f"\n  {'Residue':<22} {'BSA':>6} {'Contacts':>9} {'Cons':>6} {'Composite':>10}")
    print(f"  {'-'*57}")
    for c in contacts[:5]:
        sb = ' *SB' if c['salt_bridge'] else ''
        print(f"  {c['residue']:<22} {c['bsa']:>6.1f} {c['contacts']:>9}"
              f" {c['conservation']:>6.3f} {c['composite']:>10.4f}{sb}")
    print('='*60)


if __name__ == '__main__':
    print('CGCP Phase 2 Step 2 - Contact Mapping: NSP12-NSP7')
    iface, rows, cons_nsp12, cons_nsp7 = load_data()
    contacts = build_contact_table(iface, rows, cons_nsp12, cons_nsp7)
    save_tsv(contacts)
    out = save_json(contacts, iface)
    make_figure(contacts)
    print_summary(contacts, out)
    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Phase 2 Step 3 - feature classification")
