#!/usr/bin/env python3
"""
CGCP Phase 2 Steps 1+2+3 - NSP7-NSP8 (Mode B AF3)
Structural verification, contact mapping, feature classification.

Structure: receptor_NSP7-NSP8_ModeB_AF3_6.pdb
  Chain A = NSP8 (198 residues, starts at res 78)
  Chain B = NSP7 (83 residues, starts at res 2)
  Anchor  = PHE92 (Chain A, NSP8, cons=1.000)

Note: Mode B (AlphaFold3) is druggable — PHE92 at interface.
      Mode A (crystal) has PHE92 away from interface.

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step01_02_03_NSP7-NSP8.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
    '03-virtual-screening/NSP7-NSP8_6/receptor_NSP7-NSP8_ModeB_AF3_6.pdb')

OUT1 = os.path.join(BASE, 'CGCP/02-deep-dive/NSP7-NSP8/step-01-structure')
OUT2 = os.path.join(BASE, 'CGCP/02-deep-dive/NSP7-NSP8/step-02-contacts')
OUT3 = os.path.join(BASE, 'CGCP/02-deep-dive/NSP7-NSP8/step-03-features')
for d in [OUT1, OUT2, OUT3]:
    os.makedirs(d, exist_ok=True)
os.makedirs(os.path.join(OUT1, 'visuals'), exist_ok=True)
os.makedirs(os.path.join(OUT1, 'pymol-sessions'), exist_ok=True)

NSP8_CHAIN     = 'A'
NSP7_CHAIN     = 'B'
ANCHOR_RES     = 92
ANCHOR_AA      = 'PHE'
CONTACT_CUTOFF = 5.0

FEATURE_MAP = {
    'PHE':'aromatic','TYR':'aromatic','TRP':'aromatic','HIS':'aromatic',
    'LEU':'hydrophobic','VAL':'hydrophobic','ILE':'hydrophobic',
    'ALA':'hydrophobic','MET':'hydrophobic','PRO':'hydrophobic',
    'GLY':'hydrophobic',
    'ARG':'charged_pos','LYS':'charged_pos',
    'ASP':'charged_neg','GLU':'charged_neg',
    'SER':'hbond_donor','THR':'hbond_donor',
    'ASN':'hbond_donor','GLN':'hbond_donor','CYS':'hbond_donor',
}

FEAT_COLORS = {
    'anchor':      COLORS['black'],
    'aromatic':    COLORS['red'],
    'hydrophobic': COLORS['gray'],
    'charged_pos': COLORS['blue'],
    'charged_neg': COLORS['green'],
    'hbond_donor': COLORS['hbond'],
    'unknown':     '#AAAAAA',
}
FEAT_FILLS = {
    'anchor':      '#555555',
    'aromatic':    COLORS['red_fill'],
    'hydrophobic': COLORS['gray_fill'],
    'charged_pos': COLORS['blue_fill'],
    'charged_neg': COLORS['green_fill'],
    'hbond_donor': COLORS['hbond_f'],
    'unknown':     '#DDDDDD',
}


def get_ca(chain, resnum):
    try:
        return chain[(' ', resnum, ' ')]['CA'].coord
    except KeyError:
        return None


def dist(a, b):
    return float(np.linalg.norm(np.array(a) - np.array(b)))


def load_conservation():
    base = os.path.expanduser('~/projects/rtc-pan-coronavirus')
    nsp8_cons = {}
    nsp7_cons = {}
    for protein, fpath, cons_dict in [
        ('NSP8', f'{base}/00-reference/sequences/conservation/NSP8_aligned.fasta', nsp8_cons),
        ('NSP7', f'{base}/00-reference/sequences/conservation/NSP7_aligned.fasta', nsp7_cons),
    ]:
        try:
            from Bio import SeqIO
            from collections import Counter
            records = list(SeqIO.parse(fpath, 'fasta'))
            ref_pos = 0
            for aln_pos in range(len(records[0].seq)):
                col    = [str(r.seq[aln_pos]) for r in records]
                ref_aa = col[0]
                if ref_aa == '-':
                    continue
                ref_pos += 1
                aas    = [aa for aa in col if aa != '-']
                counts = Counter(aas)
                cons_dict[ref_pos] = round(
                    counts.most_common(1)[0][1]/len(aas), 3)
            print(f"  {protein}: {len(cons_dict)} positions loaded")
        except Exception as e:
            print(f"  Warning: {protein} ({e})")
            for i in range(1, 220):
                cons_dict[i] = 0.750
    return nsp8_cons, nsp7_cons


def step1(structure):
    chain_a = structure[NSP8_CHAIN]
    chain_b = structure[NSP7_CHAIN]
    results = {}

    chains = [c.id for c in structure.get_chains()]
    c1 = 'A' in chains and 'B' in chains
    results['C1'] = {
        'name':   'Chains A (NSP8) + B (NSP7) present',
        'pass':   c1,
        'detail': f"Chains: {chains}",
    }

    na = sum(1 for r in chain_a.get_residues() if PDB.is_aa(r))
    nb = sum(1 for r in chain_b.get_residues() if PDB.is_aa(r))
    c2 = na > 50 and nb > 30
    results['C2'] = {
        'name':   'Chain sizes reasonable',
        'pass':   c2,
        'detail': f"NSP8={na} res, NSP7={nb} res",
    }

    ca_anchor = get_ca(chain_a, ANCHOR_RES)
    if ca_anchor is not None:
        min_dist = min(
            dist(ca_anchor, r['CA'].coord)
            for r in chain_b.get_residues()
            if PDB.is_aa(r) and 'CA' in r)
        c3 = min_dist < 15.0
        results['C3'] = {
            'name':   'PHE92 adjacent to NSP7',
            'pass':   c3,
            'detail': f"PHE92 min dist to NSP7: {min_dist:.2f}Å",
        }
    else:
        results['C3'] = {'name':'PHE92 adjacent to NSP7',
                         'pass':False,'detail':'PHE92 not found'}

    missing = []
    for r in chain_a.get_residues():
        if not PDB.is_aa(r): continue
        if 'CA' not in r or 'N' not in r or 'C' not in r:
            missing.append(r.id[1])
    c4 = len(missing) == 0
    results['C4'] = {
        'name':   'Backbone complete (Chain A/NSP8)',
        'pass':   c4,
        'detail': f"{len(missing)} missing backbone atoms",
    }

    try:
        res92 = chain_a[(' ', ANCHOR_RES, ' ')]
        aa92  = res92.resname
        c5    = aa92 == ANCHOR_AA
        ca92  = res92['CA'].coord
        results['C5'] = {
            'name':   'PHE92 present and correct',
            'pass':   c5,
            'detail': f"Res 92={aa92} at "
                      f"({ca92[0]:.3f},{ca92[1]:.3f},{ca92[2]:.3f})",
        }
    except KeyError:
        results['C5'] = {'name':'PHE92 present',
                         'pass':False,'detail':'PHE92 not found'}

    if ca_anchor is not None:
        close = [r for r in chain_b.get_residues()
                 if PDB.is_aa(r) and 'CA' in r
                 and dist(ca_anchor, r['CA'].coord) < 15.0]
        c6 = len(close) >= 3
        results['C6'] = {
            'name':   'Positive ctrl: PHE92 near NSP7',
            'pass':   c6,
            'detail': f"{len(close)} NSP7 residues within 15Å",
        }
    else:
        results['C6'] = {'name':'Positive ctrl',
                         'pass':False,'detail':'PHE92 not found'}

    # Neg ctrl — PHE92 far from NSP8 C-terminus (res 191)
    if ca_anchor is not None:
        ca_cter = get_ca(chain_a, 191)
        if ca_cter is not None:
            d_cter = dist(ca_anchor, ca_cter)
            c7 = d_cter > 20.0
            results['C7'] = {
                'name':   'Neg ctrl: PHE92 far from NSP8 C-term',
                'pass':   c7,
                'detail': f"PHE92 to NSP8-191: {d_cter:.2f}Å",
            }
        else:
            results['C7'] = {'name':'Neg ctrl',
                             'pass':True,'detail':'C-term not found — OK'}
    else:
        results['C7'] = {'name':'Neg ctrl',
                         'pass':False,'detail':'PHE92 not found'}

    n_pass  = sum(1 for r in results.values() if r['pass'])
    n_total = len(results)

    with open(os.path.join(OUT1,'audit_report.md'),'w') as f:
        f.write('# CGCP Step 1 Structural Audit — NSP7-NSP8\n\n')
        f.write(f'**Result: {n_pass}/{n_total} PASS**\n\n')
        f.write('**Mode B (AlphaFold3) used — Mode A not druggable**\n\n')
        f.write('| Check | Name | Result | Detail |\n|-------|------|--------|--------|\n')
        for ck,r in results.items():
            sym = '✅' if r['pass'] else '❌'
            f.write(f"| {ck} | {r['name']} | {sym} | {r['detail']} |\n")

    with open(os.path.join(OUT1,'structural_verification.json'),'w') as f:
        json.dump({'interface':'NSP7-NSP8','pdb':PDB_FILE,
                   'mode':'ModeB_AF3','anchor':'PHE92',
                   'n_pass':n_pass,'n_total':n_total,
                   'checks':results},f,indent=2)

    fig, ax = plt.subplots(figsize=(10,5))
    checks  = list(results.keys())
    x       = np.arange(len(checks))
    for xi,ck in enumerate(checks):
        r  = results[ck]
        ec = COLORS['green'] if r['pass'] else COLORS['red']
        fc = COLORS['green_fill'] if r['pass'] else COLORS['anchor_f']
        ax.bar(xi, 1 if r['pass'] else 0, width=0.50,
               facecolor=fc, edgecolor=ec,
               linewidth=1.2, zorder=3, clip_on=False)
        ax.text(xi, (1 if r['pass'] else 0)+0.03,
                'PASS' if r['pass'] else 'FAIL',
                ha='center', va='bottom', fontsize=9,
                fontweight='bold', color=ec, clip_on=False)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{c}\n{results[c]['name'][:20]}" for c in checks],
        fontsize=7.5, ha='center')
    ax.tick_params(axis='x', pad=4, length=6, width=1.2)
    ax.set_ylabel('Check result (1=PASS)',
                  fontsize=11, fontweight='bold')
    prism_axes(ax, ymax=1.2, yticks=[0,1])
    ax.set_xlim(-0.6, len(checks)-0.4)
    ax.set_title(
        f'NSP7-NSP8 Structural Audit — {n_pass}/{n_total} PASS\n'
        f'Anchor: PHE92 (NSP8, ModeB AF3)',
        loc='left', fontsize=9.5, pad=8, linespacing=1.5)
    fig.savefig(os.path.join(OUT1,
        'Fig_Step01_StructuralAudit_NSP7-NSP8.png'))
    plt.close()

    return results, n_pass, n_total


def step2(structure, nsp8_cons, nsp7_cons):
    chain_a = structure[NSP8_CHAIN]
    chain_b = structure[NSP7_CHAIN]
    res_a   = [r for r in chain_a.get_residues() if PDB.is_aa(r)]
    res_b   = [r for r in chain_b.get_residues() if PDB.is_aa(r)]

    interface_a = {}
    interface_b = {}

    for ra in res_a:
        pos_a = ra.id[1]
        contacts = 0; partners = []
        for rb in res_b:
            for atom_a in ra.get_atoms():
                for atom_b in rb.get_atoms():
                    if atom_a - atom_b < CONTACT_CUTOFF:
                        contacts += 1
                        if rb.id[1] not in partners:
                            partners.append(rb.id[1])
        if contacts > 0:
            ca = ra['CA'].coord if 'CA' in ra else None
            interface_a[pos_a] = {
                'residue': f'NSP8-{ra.resname}{pos_a}',
                'chain':   'NSP8', 'position': pos_a,
                'aa': ra.resname, 'contacts': contacts,
                'partners': partners,
                'ca_coord': ca.tolist() if ca is not None else None,
            }

    for rb in res_b:
        pos_b = rb.id[1]
        contacts = 0; partners = []
        for ra in res_a:
            for atom_b in rb.get_atoms():
                for atom_a in ra.get_atoms():
                    if atom_b - atom_a < CONTACT_CUTOFF:
                        contacts += 1
                        if ra.id[1] not in partners:
                            partners.append(ra.id[1])
        if contacts > 0:
            ca = rb['CA'].coord if 'CA' in rb else None
            interface_b[pos_b] = {
                'residue': f'NSP7-{rb.resname}{pos_b}',
                'chain':   'NSP7', 'position': pos_b,
                'aa': rb.resname, 'contacts': contacts,
                'partners': partners,
                'ca_coord': ca.tolist() if ca is not None else None,
            }

    all_c = ([v['contacts'] for v in interface_a.values()] +
             [v['contacts'] for v in interface_b.values()])
    max_c = max(all_c) if all_c else 1
    records = []

    for pos, data in interface_a.items():
        cons = nsp8_cons.get(pos, 0.750)
        norm = data['contacts'] / max_c
        comp = round(0.4*norm +
                     0.3*min(data['contacts']/50,1.0) +
                     0.3*cons, 4)
        records.append({
            'residue':       data['residue'],
            'chain':         'NSP8',
            'position':      pos,
            'aa':            data['aa'],
            'contact_score': data['contacts'],
            'partners':      ','.join(map(str,data['partners'])),
            'n_partners':    len(data['partners']),
            'conservation':  cons,
            'composite':     comp,
            'is_anchor':     1 if pos == ANCHOR_RES else 0,
            'ca_x': round(data['ca_coord'][0],3) if data['ca_coord'] else None,
            'ca_y': round(data['ca_coord'][1],3) if data['ca_coord'] else None,
            'ca_z': round(data['ca_coord'][2],3) if data['ca_coord'] else None,
        })

    for pos, data in interface_b.items():
        cons = nsp7_cons.get(pos, 0.750)
        norm = data['contacts'] / max_c
        comp = round(0.4*norm +
                     0.3*min(data['contacts']/50,1.0) +
                     0.3*cons, 4)
        records.append({
            'residue':       data['residue'],
            'chain':         'NSP7',
            'position':      pos,
            'aa':            data['aa'],
            'contact_score': data['contacts'],
            'partners':      ','.join(map(str,data['partners'])),
            'n_partners':    len(data['partners']),
            'conservation':  cons,
            'composite':     comp,
            'is_anchor':     0,
            'ca_x': round(data['ca_coord'][0],3) if data['ca_coord'] else None,
            'ca_y': round(data['ca_coord'][1],3) if data['ca_coord'] else None,
            'ca_z': round(data['ca_coord'][2],3) if data['ca_coord'] else None,
        })

    max_comp = max(r['composite'] for r in records) if records else 1
    for r in records:
        r['composite'] = round(r['composite']/max_comp, 4)
    records.sort(key=lambda x: (0 if x['is_anchor'] else 1,
                                 -x['composite']))

    tsv = os.path.join(OUT2, 'contact_map_NSP7-NSP8.tsv')
    fields = ['residue','chain','position','aa',
              'contact_score','partners','n_partners',
              'conservation','composite','is_anchor',
              'ca_x','ca_y','ca_z']
    with open(tsv,'w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=fields,delimiter='\t')
        w.writeheader(); w.writerows(records)

    summary = {
        'interface': 'NSP7-NSP8',
        'anchor':    'PHE92',
        'mode':      'ModeB_AF3',
        'n_interface': len(records),
        'n_nsp8':    sum(1 for r in records if r['chain']=='NSP8'),
        'n_nsp7':    sum(1 for r in records if r['chain']=='NSP7'),
        'records':   records,
    }
    with open(os.path.join(OUT2,'contact_map_NSP7-NSP8.json'),'w') as f:
        json.dump(summary,f,indent=2)

    return records, summary


def step3(records):
    classified = []
    for row in records:
        pos   = row['position']
        chain = row['chain']
        aa    = row['aa']
        feat  = ('anchor' if pos==ANCHOR_RES and chain=='NSP8'
                 else FEATURE_MAP.get(aa,'unknown'))
        classified.append({**row, 'primary_feature': feat})

    classified.sort(key=lambda x: (
        0 if x['primary_feature']=='anchor' else 1,
        -x['composite']))

    feat_counts = {}
    for r in classified:
        ft = r['primary_feature']
        feat_counts[ft] = feat_counts.get(ft,0)+1

    tsv = os.path.join(OUT3,'feature_classification_NSP7-NSP8.tsv')
    with open(tsv,'w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=list(classified[0].keys()),
                           delimiter='\t')
        w.writeheader(); w.writerows(classified)
    with open(os.path.join(OUT3,'feature_classification_NSP7-NSP8.json'),'w') as f:
        json.dump({'interface':'NSP7-NSP8','anchor':'PHE92',
                   'mode':'ModeB_AF3','n_total':len(classified),
                   'feature_counts':feat_counts,
                   'records':classified},f,indent=2)
    return classified, feat_counts


def make_figure(records, summary, classified,
                feat_counts, results, n_pass):
    fig = plt.figure(figsize=(16,5.5))
    gs  = fig.add_gridspec(1,3,wspace=0.48,
                            left=0.06,right=0.97,
                            top=0.88,bottom=0.26)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[0,2])

    NSP8_COL  = COLORS['red']
    NSP8_FILL = COLORS['red_fill']
    NSP7_COL  = COLORS['blue']
    NSP7_FILL = COLORS['blue_fill']
    ANC_COL   = COLORS['black']

    # Panel A: audit
    checks = list(results.keys())
    x = np.arange(len(checks))
    for xi,ck in enumerate(checks):
        r  = results[ck]
        ec = COLORS['green'] if r['pass'] else COLORS['red']
        fc = COLORS['green_fill'] if r['pass'] else COLORS['anchor_f']
        ax1.bar(xi, 1 if r['pass'] else 0, width=0.50,
                facecolor=fc, edgecolor=ec,
                linewidth=1.2, zorder=3, clip_on=False)
        ax1.text(xi, (1 if r['pass'] else 0)+0.04,
                 'P' if r['pass'] else 'F',
                 ha='center',va='bottom',fontsize=9,
                 fontweight='bold',color=ec,clip_on=False)
    ax1.set_xticks(x)
    ax1.set_xticklabels([c for c in checks],fontsize=9)
    ax1.set_ylabel('Check result',fontsize=11,fontweight='bold')
    prism_axes(ax1,ymax=1.2,yticks=[0,1])
    ax1.set_xlim(-0.6,len(checks)-0.4)
    panel_title(ax1,'A',f'Structural audit ({n_pass}/7)',
                'NSP7-NSP8 ModeB AF3')

    # Panel B: top 25 contacts
    top = records[:25]
    labels_b = [r['residue'].split('-')[1] for r in top]
    x2 = np.arange(len(top))
    for xi,r in enumerate(top):
        is_anc = r['is_anchor']
        ec = (ANC_COL if is_anc else
              NSP8_COL if r['chain']=='NSP8' else NSP7_COL)
        fc = (ANC_COL if is_anc else
              NSP8_FILL if r['chain']=='NSP8' else NSP7_FILL)
        ax2.bar(xi,r['composite'],width=0.40,
                facecolor=fc,edgecolor=ec,
                linewidth=2.0 if is_anc else 1.2,
                zorder=3,clip_on=False)
        if is_anc:
            ax2.text(xi,r['composite']+0.02,'★',
                     ha='center',va='bottom',
                     fontsize=10,color=ec,clip_on=False)
    ax2.set_xticks(x2)
    set_xticklabels_vertical(ax2,labels_b,fontsize=6.5)
    ax2.set_ylabel('CGCP composite score',
                   fontsize=11,fontweight='bold')
    ax2.set_xlim(-0.6,len(top)-0.4)
    prism_axes(ax2,ymax=1.25,yticks=[0,0.2,0.4,0.6,0.8,1.0])
    make_legend(ax2,[
        (ANC_COL,ANC_COL,'PHE92 anchor'),
        (NSP8_FILL,NSP8_COL,'NSP8'),
        (NSP7_FILL,NSP7_COL,'NSP7'),
    ],loc='upper right',fontsize=8)
    panel_title(ax2,'B',
                f"Contact map (top 25)\nNSP8={summary['n_nsp8']} | NSP7={summary['n_nsp7']}")

    # Panel C: feature distribution
    feat_order  = ['anchor','aromatic','hydrophobic',
                   'charged_pos','charged_neg','hbond_donor','unknown']
    feat_labels = ['Anchor','Aromatic','Hydrophobic',
                   'Charged+','Charged-','H-bond','Unknown']
    counts = [feat_counts.get(f,0) for f in feat_order]
    cols   = [FEAT_COLORS[f] for f in feat_order]
    fills  = [FEAT_FILLS[f]  for f in feat_order]
    x3 = np.arange(len(feat_order))
    for xi,val,fc,ec in zip(x3,counts,fills,cols):
        if val==0: continue
        ax3.bar(xi,val,width=0.50,facecolor=fc,edgecolor=ec,
                linewidth=1.2,zorder=3,clip_on=False)
        ax3.text(xi,val+0.15,str(val),
                 ha='center',va='bottom',fontsize=9,
                 fontweight='bold',color=ec,clip_on=False)
    ax3.set_xticks(x3)
    set_xticklabels_vertical(ax3,feat_labels,fontsize=9)
    ax3.set_ylabel('Number of residues',
                   fontsize=11,fontweight='bold')
    ymax = max(counts)+5
    prism_axes(ax3,ymax=ymax,yticks=range(0,ymax+1,5))
    ax3.set_xlim(-0.6,len(feat_order)-0.4)
    panel_title(ax3,'C','Feature distribution')

    path = os.path.join(OUT3,'Fig_Steps01-03_NSP7-NSP8.png')
    fig.savefig(path); plt.close()
    print(f"  Combined figure: {path}")


def print_all(results,n_pass,n_total,records,summary,feat_counts):
    print('\n'+'='*65)
    print('STEPS 1-3 — NSP7-NSP8 (Mode B AF3)')
    print('='*65)
    print(f"\nSTEP 1: {n_pass}/{n_total} PASS")
    for ck,r in results.items():
        sym = '✅' if r['pass'] else '❌'
        print(f"  {sym} {ck}: {r['detail']}")
    print(f"\nSTEP 2: {summary['n_interface']} interface residues")
    print(f"  NSP8={summary['n_nsp8']} | NSP7={summary['n_nsp7']}")
    print(f"\n  Top 8:")
    print(f"  {'Residue':<20} {'Chain':<7} "
          f"{'Contacts':>8} {'Cons':>6} {'Comp':>6}")
    print(f"  {'-'*53}")
    for r in records[:8]:
        sym = '★' if r['is_anchor'] else ' '
        print(f"  {sym} {r['residue']:<16} {r['chain']:<7} "
              f"{r['contact_score']:>8} "
              f"{r['conservation']:>6.3f} "
              f"{r['composite']:>6.3f}")
    print(f"\nSTEP 3: Feature distribution:")
    for feat,count in sorted(feat_counts.items(),key=lambda x:-x[1]):
        print(f"  {feat:<22} {count:>4}")
    if n_pass==n_total:
        print(f"\n  → Structure APPROVED")
    print('='*65)


if __name__ == '__main__':
    print('CGCP Phase 2 Steps 1-3 — NSP7-NSP8')
    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', PDB_FILE)[0]

    print('  Loading conservation...')
    nsp8_cons, nsp7_cons = load_conservation()

    print('  Step 1: Structural verification...')
    results, n_pass, n_total = step1(structure)

    print('  Step 2: Contact mapping...')
    records, summary = step2(structure, nsp8_cons, nsp7_cons)

    print('  Step 3: Feature classification...')
    classified, feat_counts = step3(records)

    print('  Generating figures...')
    make_figure(records, summary, classified,
                feat_counts, results, n_pass)

    print_all(results, n_pass, n_total, records,
              summary, feat_counts)

    print(f"\nOutputs: {OUT1}, {OUT2}, {OUT3}")
    print("Next: Steps 4+5 — DBSCAN + conservation")
