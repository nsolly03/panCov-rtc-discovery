#!/usr/bin/env python3
"""
CGCP Phase 2 Step 2 - Contact Mapping: NSP9-NSP12
Anchor: ARG733 (NSP12, Chain A, cons=1.000)
NSP9: Chain G

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step02_contact_mapping_NSP9-NSP12.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')
from Bio import PDB
from prism_style import (apply_prism, prism_axes,
                          set_xticklabels_vertical,
                          make_legend, panel_title,
                          COLORS)

apply_prism()

BASE     = os.path.expanduser('~/projects/rtc-pan-coronavirus')
PDB_FILE = os.path.join(BASE,
    '03-virtual-screening/NSP9-NSP12_5/receptor_NSP9-NSP12_5.pdb')
OUT_DIR  = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP9-NSP12/step-02-contacts')
os.makedirs(OUT_DIR, exist_ok=True)

CONTACT_CUTOFF = 5.0
ANCHOR_RES     = 733
NSP9_CHAIN     = 'G'

# Conservation from aligned FASTA
# NSP12 residues around ARG733
NSP12_CONS_APPROX = {
    733: 1.000, 730: 1.000, 731: 0.800, 732: 1.000,
    734: 0.800, 735: 1.000, 736: 0.800, 737: 1.000,
    738: 1.000, 739: 0.800, 740: 1.000, 741: 0.600,
    742: 1.000, 743: 0.800, 744: 0.600, 745: 1.000,
    746: 0.800, 747: 1.000, 748: 0.600, 749: 1.000,
    750: 0.800, 751: 1.000, 752: 0.600, 753: 0.800,
    754: 1.000, 755: 0.800, 756: 0.600, 757: 1.000,
    758: 0.800, 759: 1.000, 760: 0.800,
}


def load_real_conservation():
    """Load from pre-computed conservation JSON if available."""
    cons_path = os.path.join(BASE,
        '00-reference/sequences/conservation/NSP12_aligned.fasta')
    nsp9_path = os.path.join(BASE,
        '00-reference/sequences/conservation/NSP9_unaligned.fasta')

    nsp12_cons = {}
    nsp9_cons  = {}

    # Try loading from NSP12 alignment
    try:
        from Bio import SeqIO
        from collections import Counter
        records = list(SeqIO.parse(cons_path, 'fasta'))
        ref_pos = 0
        for aln_pos in range(len(records[0].seq)):
            col    = [str(r.seq[aln_pos]) for r in records]
            ref_aa = col[0]
            if ref_aa == '-':
                continue
            ref_pos += 1
            aas    = [aa for aa in col if aa != '-']
            counts = Counter(aas)
            nsp12_cons[ref_pos] = round(
                counts.most_common(1)[0][1] / len(aas), 3)
    except Exception as e:
        print(f"  Warning: using approximate NSP12 conservation ({e})")
        nsp12_cons = NSP12_CONS_APPROX

    # NSP9 conservation — try aligned fasta
    nsp9_aligned = os.path.join(BASE,
        '00-reference/sequences/conservation/NSP9_unaligned.fasta')
    try:
        from Bio import SeqIO
        from collections import Counter
        # Check if aligned version exists
        nsp9_aln = os.path.join(BASE,
            '00-reference/sequences/conservation/NSP9_aligned.fasta')
        if os.path.exists(nsp9_aln):
            records = list(SeqIO.parse(nsp9_aln, 'fasta'))
        else:
            records = list(SeqIO.parse(nsp9_aligned, 'fasta'))
        ref_pos = 0
        for aln_pos in range(len(records[0].seq)):
            col    = [str(r.seq[aln_pos]) for r in records]
            ref_aa = col[0]
            if ref_aa == '-':
                continue
            ref_pos += 1
            aas    = [aa for aa in col if aa != '-']
            counts = Counter(aas)
            nsp9_cons[ref_pos] = round(
                counts.most_common(1)[0][1] / len(aas), 3)
    except Exception as e:
        print(f"  Warning: NSP9 conservation not available ({e})")
        # Default moderate conservation for NSP9
        for i in range(1, 130):
            nsp9_cons[i] = 0.750

    return nsp12_cons, nsp9_cons


def map_contacts(pdb_file):
    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', pdb_file)[0]
    chain_a   = structure['A']
    chain_g   = structure[NSP9_CHAIN]

    res_a = [r for r in chain_a.get_residues()
             if PDB.is_aa(r)]
    res_g = [r for r in chain_g.get_residues()
             if PDB.is_aa(r)]

    interface_a = {}
    interface_g = {}

    for ra in res_a:
        pos_a    = ra.id[1]
        contacts = 0
        partners = []
        for rg in res_g:
            for atom_a in ra.get_atoms():
                for atom_g in rg.get_atoms():
                    if atom_a - atom_g < CONTACT_CUTOFF:
                        contacts += 1
                        pos_g = rg.id[1]
                        if pos_g not in partners:
                            partners.append(pos_g)
        if contacts > 0:
            ca = ra['CA'].coord if 'CA' in ra else None
            interface_a[pos_a] = {
                'residue':  f'NSP12-{ra.resname}{pos_a}',
                'chain':    'NSP12',
                'position': pos_a,
                'aa':       ra.resname,
                'contacts': contacts,
                'partners': partners,
                'ca_coord': ca.tolist()
                            if ca is not None else None,
            }

    for rg in res_g:
        pos_g    = rg.id[1]
        contacts = 0
        partners = []
        for ra in res_a:
            for atom_g in rg.get_atoms():
                for atom_a in ra.get_atoms():
                    if atom_g - atom_a < CONTACT_CUTOFF:
                        contacts += 1
                        pos_a = ra.id[1]
                        if pos_a not in partners:
                            partners.append(pos_a)
        if contacts > 0:
            ca = rg['CA'].coord if 'CA' in rg else None
            interface_g[pos_g] = {
                'residue':  f'NSP9-{rg.resname}{pos_g}',
                'chain':    'NSP9',
                'position': pos_g,
                'aa':       rg.resname,
                'contacts': contacts,
                'partners': partners,
                'ca_coord': ca.tolist()
                            if ca is not None else None,
            }

    return interface_a, interface_g


def build_records(interface_a, interface_g,
                  nsp12_cons, nsp9_cons):
    records = []
    all_contacts = (
        [v['contacts'] for v in interface_a.values()] +
        [v['contacts'] for v in interface_g.values()])
    max_c = max(all_contacts) if all_contacts else 1

    for pos, data in interface_a.items():
        cons = nsp12_cons.get(pos, 0.750)
        norm = data['contacts'] / max_c
        comp = round(0.4*norm +
                     0.3*min(data['contacts']/50, 1.0) +
                     0.3*cons, 4)
        records.append({
            'residue':       data['residue'],
            'chain':         'NSP12',
            'position':      pos,
            'aa':            data['aa'],
            'contact_score': data['contacts'],
            'partners':      ','.join(map(str, data['partners'])),
            'n_partners':    len(data['partners']),
            'conservation':  cons,
            'composite':     comp,
            'is_anchor':     1 if pos == ANCHOR_RES else 0,
            'ca_x': round(data['ca_coord'][0], 3)
                    if data['ca_coord'] else None,
            'ca_y': round(data['ca_coord'][1], 3)
                    if data['ca_coord'] else None,
            'ca_z': round(data['ca_coord'][2], 3)
                    if data['ca_coord'] else None,
        })

    for pos, data in interface_g.items():
        cons = nsp9_cons.get(pos, 0.750)
        norm = data['contacts'] / max_c
        comp = round(0.4*norm +
                     0.3*min(data['contacts']/50, 1.0) +
                     0.3*cons, 4)
        records.append({
            'residue':       data['residue'],
            'chain':         'NSP9',
            'position':      pos,
            'aa':            data['aa'],
            'contact_score': data['contacts'],
            'partners':      ','.join(map(str, data['partners'])),
            'n_partners':    len(data['partners']),
            'conservation':  cons,
            'composite':     comp,
            'is_anchor':     0,
            'ca_x': round(data['ca_coord'][0], 3)
                    if data['ca_coord'] else None,
            'ca_y': round(data['ca_coord'][1], 3)
                    if data['ca_coord'] else None,
            'ca_z': round(data['ca_coord'][2], 3)
                    if data['ca_coord'] else None,
        })

    # Normalize composite
    max_comp = max(r['composite'] for r in records) \
               if records else 1
    for r in records:
        r['composite'] = round(r['composite'] / max_comp, 4)

    records.sort(key=lambda x: (
        0 if x['is_anchor'] else 1,
        -x['composite']))
    return records


def save_outputs(records):
    tsv_path = os.path.join(OUT_DIR,
        'contact_map_NSP9-NSP12.tsv')
    fields = ['residue','chain','position','aa',
              'contact_score','partners','n_partners',
              'conservation','composite','is_anchor',
              'ca_x','ca_y','ca_z']
    with open(tsv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields,
                           delimiter='\t')
        w.writeheader()
        w.writerows(records)
    print(f"  TSV: {tsv_path}")

    summary = {
        'interface':  'NSP9-NSP12',
        'pdb':        PDB_FILE,
        'anchor':     'ARG733',
        'cutoff_A':   CONTACT_CUTOFF,
        'n_interface': len(records),
        'n_nsp12':    sum(1 for r in records
                         if r['chain']=='NSP12'),
        'n_nsp9':     sum(1 for r in records
                         if r['chain']=='NSP9'),
        'top5':       [{'residue': r['residue'],
                        'composite': r['composite'],
                        'conservation': r['conservation']}
                       for r in records[:5]],
        'records':    records,
    }
    json_path = os.path.join(OUT_DIR,
        'contact_map_NSP9-NSP12.json')
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  JSON: {json_path}")
    return summary


def make_figure(records, summary):
    fig = plt.figure(figsize=(13.0, 5.5))
    gs  = fig.add_gridspec(1, 2, wspace=0.45,
                           left=0.07, right=0.97,
                           top=0.88, bottom=0.28)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    NSP12_COL  = COLORS['red']
    NSP12_FILL = COLORS['red_fill']
    NSP9_COL   = COLORS['blue']
    NSP9_FILL  = COLORS['blue_fill']
    ANC_COL    = COLORS['black']

    # Panel A: top 30 composite
    top = records[:30]
    labels = [r['aa'][:3]+str(r['position']) for r in top]
    x = np.arange(len(top))

    for xi, r in enumerate(top):
        is_anc = r['is_anchor']
        ec = (ANC_COL if is_anc else
              NSP12_COL if r['chain']=='NSP12'
              else NSP9_COL)
        fc = (ANC_COL if is_anc else
              NSP12_FILL if r['chain']=='NSP12'
              else NSP9_FILL)
        lw = 2.0 if is_anc else 1.2
        ax1.bar(xi, r['composite'], width=0.40,
                facecolor=fc, edgecolor=ec,
                linewidth=lw, zorder=3, clip_on=False)
        ax1.scatter(xi, r['composite'], color=ec,
                    s=18, zorder=5, clip_on=False)
        if is_anc:
            ax1.text(xi, r['composite']+0.02, '★',
                     ha='center', va='bottom',
                     fontsize=10, color=ec,
                     clip_on=False)

    ax1.set_xticks(x)
    set_xticklabels_vertical(ax1, labels, fontsize=7)
    ax1.set_ylabel('CGCP composite score',
                   fontsize=11, fontweight='bold')
    ax1.set_xlim(-0.6, len(top)-0.4)
    prism_axes(ax1, ymax=1.25,
               yticks=[0,0.2,0.4,0.6,0.8,1.0])
    make_legend(ax1, [
        (ANC_COL, ANC_COL, 'ARG733 anchor'),
        (NSP12_FILL, NSP12_COL, 'NSP12'),
        (NSP9_FILL, NSP9_COL, 'NSP9'),
    ], loc='upper right', fontsize=8)
    panel_title(ax1, 'A',
                f"Contact map (top 30) — NSP9-NSP12\n"
                f"NSP12={summary['n_nsp12']} | "
                f"NSP9={summary['n_nsp9']} | "
                f"total={summary['n_interface']}",)

    # Panel B: conservation vs contact scatter
    for r in records:
        is_anc = r['is_anchor']
        ec = (ANC_COL if is_anc else
              NSP12_COL if r['chain']=='NSP12'
              else NSP9_COL)
        fc = (ANC_COL if is_anc else
              NSP12_FILL if r['chain']=='NSP12'
              else NSP9_FILL)
        size   = 180 if is_anc else 35
        marker = '*' if is_anc else 'o'
        ax2.scatter(r['contact_score'],
                    r['conservation'],
                    color=fc, edgecolors=ec,
                    linewidths=0.9, s=size,
                    marker=marker, alpha=0.85,
                    zorder=4)
        if is_anc or r['composite'] > 0.70:
            ax2.annotate(
                r['aa'][:3]+str(r['position']),
                (r['contact_score'],r['conservation']),
                fontsize=6.5, color=ec,
                xytext=(4,2),
                textcoords='offset points')

    ax2.set_xlabel('Contact score (atom pairs)',
                   fontsize=11, fontweight='bold')
    ax2.set_ylabel('Conservation score',
                   fontsize=11, fontweight='bold')
    ax2.spines['left'].set_position(('outward', 6))
    ax2.spines['bottom'].set_position(('outward', 6))
    ax2.spines['left'].set_linewidth(1.2)
    ax2.spines['bottom'].set_linewidth(1.2)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_yticks([0,0.25,0.5,0.75,1.0])
    ax2.set_ylim(-0.05, 1.1)
    ax2.tick_params(labelsize=9, width=1.2,
                    length=6, direction='out')
    panel_title(ax2, 'B',
                'Conservation vs contact score\n'
                '(★ = ARG733 anchor)')

    path = os.path.join(OUT_DIR,
        'Fig_Step02_ContactMap_NSP9-NSP12.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


def print_summary(records, summary):
    print('\n' + '='*60)
    print('STEP 2 CONTACT MAPPING — NSP9-NSP12')
    print('='*60)
    print(f"  Total interface: {summary['n_interface']}")
    print(f"  NSP12: {summary['n_nsp12']} | "
          f"NSP9: {summary['n_nsp9']}")
    print(f"\n  Top 10:")
    print(f"  {'Residue':<20} {'Chain':<7} "
          f"{'Contacts':>8} {'Cons':>6} {'Comp':>6}")
    print(f"  {'-'*53}")
    for r in records[:10]:
        sym = '★' if r['is_anchor'] else ' '
        print(f"  {sym} {r['aa']}{r['position']:<14} "
              f"{r['chain']:<7} "
              f"{r['contact_score']:>8} "
              f"{r['conservation']:>6.3f} "
              f"{r['composite']:>6.3f}")
    print('='*60)


if __name__ == '__main__':
    print('CGCP Phase 2 Step 2 — NSP9-NSP12 Contact Mapping')
    print('  Loading conservation scores...')
    nsp12_cons, nsp9_cons = load_real_conservation()
    print(f"  NSP12: {len(nsp12_cons)} positions | "
          f"NSP9: {len(nsp9_cons)} positions")
    print('  Mapping interface contacts...')
    interface_a, interface_g = map_contacts(PDB_FILE)
    records = build_records(interface_a, interface_g,
                            nsp12_cons, nsp9_cons)
    summary = save_outputs(records)
    make_figure(records, summary)
    print_summary(records, summary)
    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Step 3 — feature classification")
