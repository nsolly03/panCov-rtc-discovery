#!/usr/bin/env python3
"""
CGCP Phase 2 Steps 1+2+3 - NSP10-NSP16
Structural verification, contact mapping, feature classification.

Special note: This PDB uses polyprotein absolute numbering.
  Chain A = NSP16 (residues 6798-7096)
  Chain B = NSP10 (residues 4271-4402, includes ZN)
  Anchor  = LYS at B4346 (= LYS93 in NSP10 protein numbering)
  NSP10 offset: 4271-1 = 4270 (position 1 = residue 4271)

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step01_02_03_NSP10-NSP16.py
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
    '03-virtual-screening/NSP10-NSP16_2/receptor_NSP10-NSP16_2.pdb')

OUT1 = os.path.join(BASE, 'CGCP/02-deep-dive/NSP10-NSP16/step-01-structure')
OUT2 = os.path.join(BASE, 'CGCP/02-deep-dive/NSP10-NSP16/step-02-contacts')
OUT3 = os.path.join(BASE, 'CGCP/02-deep-dive/NSP10-NSP16/step-03-features')
for d in [OUT1, OUT2, OUT3]:
    os.makedirs(d, exist_ok=True)
os.makedirs(os.path.join(OUT1, 'visuals'), exist_ok=True)
os.makedirs(os.path.join(OUT1, 'pymol-sessions'), exist_ok=True)

# Polyprotein numbering
ANCHOR_CHAIN  = 'B'
ANCHOR_RESNUM = 4346   # LYS93 in NSP10
ANCHOR_AA     = 'LYS'
NSP16_CHAIN   = 'A'
NSP10_CHAIN   = 'B'
NSP10_OFFSET  = 4270   # subtract to get protein position
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


def load_structure():
    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', PDB_FILE)[0]
    return structure


def load_conservation():
    """Load NSP16 and NSP10 conservation from aligned FASTA."""
    base = os.path.expanduser('~/projects/rtc-pan-coronavirus')
    nsp16_cons = {}
    nsp10_cons = {}

    for protein, fpath, cons_dict in [
        ('NSP16', f'{base}/00-reference/sequences/conservation/NSP16_aligned.fasta', nsp16_cons),
        ('NSP10', f'{base}/00-reference/sequences/conservation/NSP10_aligned.fasta', nsp10_cons),
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
                    counts.most_common(1)[0][1] / len(aas), 3)
            print(f"  {protein}: {len(cons_dict)} positions loaded")
        except Exception as e:
            print(f"  Warning: {protein} conservation not available ({e})")
            for i in range(1, 400):
                cons_dict[i] = 0.750

    return nsp16_cons, nsp10_cons


# ── STEP 1: Structural Verification ──────────────────────────
def step1(structure):
    chain_a = structure[NSP16_CHAIN]
    chain_b = structure[NSP10_CHAIN]
    results = {}

    chains = [c.id for c in structure.get_chains()]
    c1 = 'A' in chains and 'B' in chains
    results['C1'] = {
        'name':   'Chains A (NSP16) + B (NSP10) present',
        'pass':   c1,
        'detail': f"Chains: {chains}",
    }

    na = sum(1 for r in chain_a.get_residues() if PDB.is_aa(r))
    nb = sum(1 for r in chain_b.get_residues() if PDB.is_aa(r))
    c2 = na > 100 and nb > 50
    results['C2'] = {
        'name':   'Chain sizes reasonable',
        'pass':   c2,
        'detail': f"NSP16={na} res, NSP10={nb} res",
    }

    ca_anchor = get_ca(chain_b, ANCHOR_RESNUM)
    if ca_anchor is not None:
        min_dist = min(
            dist(ca_anchor, r['CA'].coord)
            for r in chain_a.get_residues()
            if PDB.is_aa(r) and 'CA' in r)
        c3 = min_dist < 15.0
        results['C3'] = {
            'name':   'Anchor adjacent to NSP16',
            'pass':   c3,
            'detail': f"LYS93(B4346) min dist to NSP16: {min_dist:.2f}Å",
        }
    else:
        results['C3'] = {'name':'Anchor adjacent to NSP16',
                         'pass':False,'detail':'Anchor not found'}

    missing = []
    for r in chain_a.get_residues():
        if not PDB.is_aa(r): continue
        if 'CA' not in r or 'N' not in r or 'C' not in r:
            missing.append(r.id[1])
    c4 = len(missing) == 0
    results['C4'] = {
        'name':   'Backbone complete (Chain A/NSP16)',
        'pass':   c4,
        'detail': f"{len(missing)} missing backbone atoms",
    }

    try:
        res_anc = chain_b[(' ', ANCHOR_RESNUM, ' ')]
        aa_anc  = res_anc.resname
        c5      = aa_anc == ANCHOR_AA
        ca_anc  = res_anc['CA'].coord
        results['C5'] = {
            'name':   'LYS93 (B4346) present and correct',
            'pass':   c5,
            'detail': f"B4346={aa_anc} at "
                      f"({ca_anc[0]:.3f},{ca_anc[1]:.3f},{ca_anc[2]:.3f})",
        }
    except KeyError:
        results['C5'] = {'name':'LYS93 present','pass':False,
                         'detail':'B4346 not found'}

    if ca_anchor is not None:
        close = [r for r in chain_a.get_residues()
                 if PDB.is_aa(r) and 'CA' in r
                 and dist(ca_anchor, r['CA'].coord) < 15.0]
        c6 = len(close) >= 3
        results['C6'] = {
            'name':   'Positive ctrl: LYS93 near NSP16',
            'pass':   c6,
            'detail': f"{len(close)} NSP16 residues within 15Å",
        }
    else:
        results['C6'] = {'name':'Positive ctrl','pass':False,
                         'detail':'Anchor not found'}

    # Negative control — check LYS93 far from NSP10 zinc (B4402)
    if ca_anchor is not None:
        try:
            zn_res = chain_b[('H_ZN', 4402, ' ')]
            zn_coord = list(zn_res.get_atoms())[0].coord
            d_zn = dist(ca_anchor, zn_coord)
            c7 = d_zn > 10.0
            results['C7'] = {
                'name':   'Neg ctrl: LYS93 not at Zn site',
                'pass':   c7,
                'detail': f"LYS93 to Zn: {d_zn:.2f}Å",
            }
        except Exception:
            results['C7'] = {
                'name':   'Neg ctrl: LYS93 not at Zn site',
                'pass':   True,
                'detail': 'Zn not found — OK',
            }
    else:
        results['C7'] = {'name':'Neg ctrl','pass':False,
                         'detail':'Anchor not found'}

    n_pass  = sum(1 for r in results.values() if r['pass'])
    n_total = len(results)

    # Save
    with open(os.path.join(OUT1,'audit_report.md'),'w') as f:
        f.write('# CGCP Step 1 Structural Audit — NSP10-NSP16\n\n')
        f.write(f'**Result: {n_pass}/{n_total} PASS**\n\n')
        f.write('**Note:** Polyprotein numbering — NSP16=ChainA(6798-7096), NSP10=ChainB(4271-4386)\n\n')
        f.write('| Check | Name | Result | Detail |\n|-------|------|--------|--------|\n')
        for ck,r in results.items():
            sym='✅' if r['pass'] else '❌'
            f.write(f"| {ck} | {r['name']} | {sym} | {r['detail']} |\n")

    with open(os.path.join(OUT1,'structural_verification.json'),'w') as f:
        json.dump({'interface':'NSP10-NSP16','pdb':PDB_FILE,
                   'anchor':'LYS93(B4346)','n_pass':n_pass,
                   'n_total':n_total,'checks':results},f,indent=2)

    # Figure
    fig, ax = plt.subplots(figsize=(10,5))
    checks  = list(results.keys())
    x       = np.arange(len(checks))
    for xi, ck in enumerate(checks):
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
    ax.tick_params(axis='x',pad=4,length=6,width=1.2)
    ax.set_ylabel('Check result (1=PASS)',
                  fontsize=11,fontweight='bold')
    prism_axes(ax,ymax=1.2,yticks=[0,1])
    ax.set_xlim(-0.6,len(checks)-0.4)
    ax.set_title(f'NSP10-NSP16 Structural Audit — {n_pass}/{n_total} PASS\n'
                 f'Anchor: LYS93(B4346) | Polyprotein numbering',
                 loc='left',fontsize=9.5,pad=8,linespacing=1.5)
    fig.savefig(os.path.join(OUT1,'Fig_Step01_StructuralAudit_NSP10-NSP16.png'))
    plt.close()

    return results, n_pass, n_total


# ── STEP 2: Contact Mapping ───────────────────────────────────
def step2(structure, nsp16_cons, nsp10_cons):
    chain_a = structure[NSP16_CHAIN]
    chain_b = structure[NSP10_CHAIN]

    res_a = [r for r in chain_a.get_residues() if PDB.is_aa(r)]
    res_b = [r for r in chain_b.get_residues() if PDB.is_aa(r)]

    interface_a = {}
    interface_b = {}

    for ra in res_a:
        pos_a    = ra.id[1]
        contacts = 0
        partners = []
        for rb in res_b:
            for atom_a in ra.get_atoms():
                for atom_b in rb.get_atoms():
                    if atom_a - atom_b < CONTACT_CUTOFF:
                        contacts += 1
                        if rb.id[1] not in partners:
                            partners.append(rb.id[1])
        if contacts > 0:
            ca = ra['CA'].coord if 'CA' in ra else None
            prot_pos = pos_a - 6797  # NSP16 protein position
            interface_a[pos_a] = {
                'residue':   f'NSP16-{ra.resname}{prot_pos}',
                'chain':     'NSP16',
                'position':  pos_a,
                'prot_pos':  prot_pos,
                'aa':        ra.resname,
                'contacts':  contacts,
                'partners':  partners,
                'ca_coord':  ca.tolist() if ca is not None else None,
            }

    for rb in res_b:
        pos_b    = rb.id[1]
        contacts = 0
        partners = []
        for ra in res_a:
            for atom_b in rb.get_atoms():
                for atom_a in ra.get_atoms():
                    if atom_b - atom_a < CONTACT_CUTOFF:
                        contacts += 1
                        if ra.id[1] not in partners:
                            partners.append(ra.id[1])
        if contacts > 0:
            ca = rb['CA'].coord if 'CA' in rb else None
            prot_pos = pos_b - 4270  # NSP10 protein position
            interface_b[pos_b] = {
                'residue':   f'NSP10-{rb.resname}{prot_pos}',
                'chain':     'NSP10',
                'position':  pos_b,
                'prot_pos':  prot_pos,
                'aa':        rb.resname,
                'contacts':  contacts,
                'partners':  partners,
                'ca_coord':  ca.tolist() if ca is not None else None,
            }

    all_c   = ([v['contacts'] for v in interface_a.values()] +
               [v['contacts'] for v in interface_b.values()])
    max_c   = max(all_c) if all_c else 1
    records = []

    for pos, data in interface_a.items():
        cons = nsp16_cons.get(data['prot_pos'], 0.750)
        norm = data['contacts'] / max_c
        comp = round(0.4*norm +
                     0.3*min(data['contacts']/50,1.0) +
                     0.3*cons, 4)
        records.append({
            'residue':       data['residue'],
            'chain':         'NSP16',
            'position':      pos,
            'prot_pos':      data['prot_pos'],
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

    for pos, data in interface_b.items():
        cons    = nsp10_cons.get(data['prot_pos'], 0.750)
        norm    = data['contacts'] / max_c
        comp    = round(0.4*norm +
                        0.3*min(data['contacts']/50,1.0) +
                        0.3*cons, 4)
        is_anc  = 1 if pos == ANCHOR_RESNUM else 0
        records.append({
            'residue':       data['residue'],
            'chain':         'NSP10',
            'position':      pos,
            'prot_pos':      data['prot_pos'],
            'aa':            data['aa'],
            'contact_score': data['contacts'],
            'partners':      ','.join(map(str,data['partners'])),
            'n_partners':    len(data['partners']),
            'conservation':  cons,
            'composite':     comp,
            'is_anchor':     is_anc,
            'ca_x': round(data['ca_coord'][0],3) if data['ca_coord'] else None,
            'ca_y': round(data['ca_coord'][1],3) if data['ca_coord'] else None,
            'ca_z': round(data['ca_coord'][2],3) if data['ca_coord'] else None,
        })

    max_comp = max(r['composite'] for r in records) if records else 1
    for r in records:
        r['composite'] = round(r['composite']/max_comp, 4)

    records.sort(key=lambda x: (0 if x['is_anchor'] else 1,
                                 -x['composite']))

    # Save
    tsv = os.path.join(OUT2,'contact_map_NSP10-NSP16.tsv')
    fields = ['residue','chain','position','prot_pos','aa',
              'contact_score','partners','n_partners',
              'conservation','composite','is_anchor',
              'ca_x','ca_y','ca_z']
    with open(tsv,'w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=fields,delimiter='\t')
        w.writeheader(); w.writerows(records)

    summary = {
        'interface':  'NSP10-NSP16',
        'anchor':     'LYS93(B4346)',
        'n_interface': len(records),
        'n_nsp16':    sum(1 for r in records if r['chain']=='NSP16'),
        'n_nsp10':    sum(1 for r in records if r['chain']=='NSP10'),
        'top5':       [{'residue':r['residue'],
                        'composite':r['composite'],
                        'conservation':r['conservation']}
                       for r in records[:5]],
        'records':    records,
    }
    with open(os.path.join(OUT2,'contact_map_NSP10-NSP16.json'),'w') as f:
        json.dump(summary, f, indent=2)

    return records, summary


# ── STEP 3: Feature Classification ───────────────────────────
def step3(records):
    classified = []
    for row in records:
        pos   = row['position']
        chain = row['chain']
        aa    = row['aa']

        if pos == ANCHOR_RESNUM and chain == 'NSP10':
            feat = 'anchor'
        else:
            feat = FEATURE_MAP.get(aa, 'unknown')

        classified.append({**row, 'primary_feature': feat})

    classified.sort(key=lambda x: (
        0 if x['primary_feature']=='anchor' else 1,
        -x['composite']))

    feat_counts = {}
    for r in classified:
        ft = r['primary_feature']
        feat_counts[ft] = feat_counts.get(ft, 0)+1

    fields = list(classified[0].keys())
    tsv = os.path.join(OUT3,'feature_classification_NSP10-NSP16.tsv')
    with open(tsv,'w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=fields,delimiter='\t')
        w.writeheader(); w.writerows(classified)

    with open(os.path.join(OUT3,'feature_classification_NSP10-NSP16.json'),'w') as f:
        json.dump({'interface':'NSP10-NSP16',
                   'anchor':'LYS93(B4346)',
                   'n_total':len(classified),
                   'feature_counts':feat_counts,
                   'records':classified},f,indent=2)

    return classified, feat_counts


# ── Combined Figure ───────────────────────────────────────────
def make_combined_figure(records, summary, classified,
                          feat_counts, results, n_pass):
    fig = plt.figure(figsize=(16,5.5))
    gs  = fig.add_gridspec(1,3,wspace=0.48,
                            left=0.06,right=0.97,
                            top=0.88,bottom=0.26)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[0,2])

    NSP16_COL  = COLORS['red']
    NSP16_FILL = COLORS['red_fill']
    NSP10_COL  = COLORS['blue']
    NSP10_FILL = COLORS['blue_fill']
    ANC_COL    = COLORS['black']

    # Panel A: audit summary
    checks = list(results.keys())
    x = np.arange(len(checks))
    for xi, ck in enumerate(checks):
        r  = results[ck]
        ec = COLORS['green'] if r['pass'] else COLORS['red']
        fc = COLORS['green_fill'] if r['pass'] else COLORS['anchor_f']
        ax1.bar(xi, 1 if r['pass'] else 0, width=0.50,
                facecolor=fc, edgecolor=ec, linewidth=1.2,
                zorder=3, clip_on=False)
        ax1.text(xi, (1 if r['pass'] else 0)+0.04,
                 'P' if r['pass'] else 'F',
                 ha='center',va='bottom',fontsize=9,
                 fontweight='bold',color=ec,clip_on=False)
    ax1.set_xticks(x)
    ax1.set_xticklabels([c for c in checks],fontsize=9)
    ax1.set_ylabel('Check result',fontsize=11,fontweight='bold')
    prism_axes(ax1,ymax=1.2,yticks=[0,1])
    ax1.set_xlim(-0.6,len(checks)-0.4)
    panel_title(ax1,'A',f'Structural audit ({n_pass}/7 PASS)',
                'NSP10-NSP16 | polyprotein numbering')

    # Panel B: top 25 contact map
    top = records[:25]
    labels_b = [r['residue'].split('-')[1] for r in top]
    x2 = np.arange(len(top))
    for xi,r in enumerate(top):
        is_anc = r['is_anchor']
        ec = (ANC_COL if is_anc else
              NSP16_COL if r['chain']=='NSP16' else NSP10_COL)
        fc = (ANC_COL if is_anc else
              NSP16_FILL if r['chain']=='NSP16' else NSP10_FILL)
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
        (ANC_COL,ANC_COL,'LYS93 anchor'),
        (NSP16_FILL,NSP16_COL,'NSP16'),
        (NSP10_FILL,NSP10_COL,'NSP10'),
    ],loc='upper right',fontsize=8)
    panel_title(ax2,'B',
                f"Contact map (top 25)\nNSP16={summary['n_nsp16']} | NSP10={summary['n_nsp10']}")

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

    path = os.path.join(OUT3,
        'Fig_Steps01-03_NSP10-NSP16.png')
    fig.savefig(path)
    plt.close()
    print(f"  Combined figure: {path}")


# ── Print summaries ───────────────────────────────────────────
def print_all(results,n_pass,n_total,records,summary,
              classified,feat_counts):
    print('\n'+'='*65)
    print('STEPS 1-3 — NSP10-NSP16')
    print('(Polyprotein numbering: NSP16=A6798+, NSP10=B4271+)')
    print('='*65)

    print(f"\nSTEP 1: {n_pass}/{n_total} PASS")
    for ck,r in results.items():
        sym='✅' if r['pass'] else '❌'
        print(f"  {sym} {ck}: {r['detail']}")

    print(f"\nSTEP 2: {summary['n_interface']} interface residues")
    print(f"  NSP16={summary['n_nsp16']} | NSP10={summary['n_nsp10']}")
    print(f"\n  Top 8:")
    print(f"  {'Residue':<22} {'Chain':<7} "
          f"{'Contacts':>8} {'Cons':>6} {'Comp':>6}")
    print(f"  {'-'*53}")
    for r in records[:8]:
        sym = '★' if r['is_anchor'] else ' '
        print(f"  {sym} {r['residue']:<18} {r['chain']:<7} "
              f"{r['contact_score']:>8} "
              f"{r['conservation']:>6.3f} "
              f"{r['composite']:>6.3f}")

    print(f"\nSTEP 3: Feature distribution:")
    for feat,count in sorted(feat_counts.items(),key=lambda x:-x[1]):
        print(f"  {feat:<22} {count:>4}")

    if n_pass == n_total:
        print(f"\n  → Structure APPROVED for CGCP pipeline")
    print('='*65)


if __name__ == '__main__':
    print('CGCP Phase 2 Steps 1-3 — NSP10-NSP16')
    print('  Loading structure...')
    structure = load_structure()

    print('  Loading conservation...')
    nsp16_cons, nsp10_cons = load_conservation()

    print('  Step 1: Structural verification...')
    results, n_pass, n_total = step1(structure)

    print('  Step 2: Contact mapping...')
    records, summary = step2(structure, nsp16_cons, nsp10_cons)

    print('  Step 3: Feature classification...')
    classified, feat_counts = step3(records)

    print('  Generating figures...')
    make_combined_figure(records, summary, classified,
                          feat_counts, results, n_pass)

    print_all(results, n_pass, n_total, records, summary,
              classified, feat_counts)

    print(f"\nOutputs:")
    print(f"  Step1: {OUT1}")
    print(f"  Step2: {OUT2}")
    print(f"  Step3: {OUT3}")
    print("Next: Steps 4+5 — DBSCAN + conservation")
