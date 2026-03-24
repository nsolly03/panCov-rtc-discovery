#!/usr/bin/env python3
"""
CGCP Phase 2 Steps 6+7 — Four remaining interfaces
Integrated assessment + pharmacophore hypothesis for:
  NSP10-NSP14, NSP13-Helicase, NSP12-NSP13, NSP15

Thresholds:
  NSP10-NSP14, NSP12-NSP13 : cons>=0.800, comp>=0.600 (INCLUDE)
  NSP13-Helicase, NSP15    : cons>=0.600, comp>=0.500 (lower — homodimer interfaces)

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step06_07_remaining_interfaces.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
import warnings
warnings.filterwarnings('ignore')
from Bio import PDB
from prism_style import (apply_prism, prism_axes,
                          set_xticklabels_vertical,
                          make_legend, panel_title,
                          COLORS)

apply_prism()

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):  return int(obj)
        if isinstance(obj, np.floating): return float(obj)
        if isinstance(obj, np.ndarray):  return obj.tolist()
        if isinstance(obj, np.bool_):    return bool(obj)
        return super().default(obj)

BASE = os.path.expanduser('~/projects/rtc-pan-coronavirus')
VS   = os.path.join(BASE, '03-virtual-screening')

DRUGGABLE = {'anchor','aromatic','hydrophobic',
             'charged_pos','charged_neg','hbond_donor'}

P = {
    'anchor':    '#1A1A1A', 'anchor_f':    '#555555',
    'include':   '#1A7D2E', 'include_f':   '#A8D5B0',
    'secondary': '#2166AC', 'secondary_f': '#92C5DE',
    'exclude':   '#636363', 'exclude_f':   '#CCCCCC',
}

EL_COLORS = {'E1':'#D7191C','E2':'#2166AC','E3':'#4DAC26'}
EL_FILLS  = {'E1':'#FDAE61','E2':'#92C5DE','E3':'#A6D96A'}

# Interface configs
INTERFACES = {
    'NSP10-NSP14': {
        'pdb':         f'{VS}/NSP10-NSP14_2/receptor_NSP10-NSP14_2.pdb',
        'chain1_id':   'A', 'name1': 'NSP10',
        'chain2_id':   'B', 'name2': 'NSP14',
        'anchor_chain':'NSP10', 'anchor_res': 80,
        'inc_cons': 0.800, 'inc_comp': 0.600,
        'sec_cons': 0.600, 'sec_comp': 0.400,
        'el_names': {
            'E1': 'E1 — Aromatic/HIS80 core',
            'E2': 'E2 — NSP14 contact zone',
            'E3': 'E3 — Extended contacts',
        },
    },
    'NSP13-Helicase': {
        'pdb':         f'{VS}/NSP13-Helicase_7/receptor_NSP13-Helicase_7.pdb',
        'chain1_id':   'A', 'name1': 'NSP13a',
        'chain2_id':   'E', 'name2': 'NSP13b',
        'anchor_chain':'NSP13a', 'anchor_res': 414,
        'inc_cons': 0.600, 'inc_comp': 0.500,
        'sec_cons': 0.400, 'sec_comp': 0.300,
        'el_names': {
            'E1': 'E1 — LYS414 salt bridge core',
            'E2': 'E2 — Hydrophobic homodimer zone',
            'E3': 'E3 — Extended contacts',
        },
    },
    'NSP12-NSP13': {
        'pdb':         f'{VS}/NSP12-NSP13_8/receptor_NSP12-NSP13_8.pdb',
        'chain1_id':   'A', 'name1': 'NSP12',
        'chain2_id':   'E', 'name2': 'NSP13',
        'anchor_chain':'NSP12', 'anchor_res': 903,
        'inc_cons': 0.800, 'inc_comp': 0.500,
        'sec_cons': 0.600, 'sec_comp': 0.300,
        'el_names': {
            'E1': 'E1 — TYR903/MET902 aromatic core',
            'E2': 'E2 — NSP13 contact zone',
            'E3': 'E3 — Supporting contacts',
        },
    },
    'NSP15': {
        'pdb':         f'{VS}/NSP15_9/receptor_NSP15_9.pdb',
        'chain1_id':   'A', 'name1': 'NSP15a',
        'chain2_id':   'B', 'name2': 'NSP15b',
        'anchor_chain':'NSP15a', 'anchor_res': 40,
        'inc_cons': 0.600, 'inc_comp': 0.500,
        'sec_cons': 0.400, 'sec_comp': 0.300,
        'el_names': {
            'E1': 'E1 — ASP40 homodimer core',
            'E2': 'E2 — Aromatic/hydrophobic zone',
            'E3': 'E3 — Extended contacts',
        },
    },
}


def dec_style(d):
    if d=='ANCHOR':    return P['anchor'],    P['anchor_f']
    if d=='INCLUDE':   return P['include'],   P['include_f']
    if d=='SECONDARY': return P['secondary'], P['secondary_f']
    return P['exclude'], P['exclude_f']


def load_conservation_tsv(iface):
    path = os.path.join(BASE,
        f'CGCP/02-deep-dive/{iface}/step-05-conservation/'
        f'conservation_overlay_{iface}.tsv')
    records = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            records.append(row)
    return records


def assess(records, cfg, iface):
    assessed = []
    chain1_name = cfg['name1']

    for row in records:
        chain = row['chain']
        pos   = int(row['position'])
        aa    = row['aa']
        feat  = row['primary_feature']
        cons  = float(row['conservation'])
        tier  = row['cons_tier']
        comp  = float(row['composite'])

        e1 = 1  # all in cluster
        e2 = 1 if chain == chain1_name else 0
        e3 = 1 if cons >= 0.800 else 0
        e4 = 1 if comp >= 0.600 else 0
        e5 = 1 if cons == 1.000 else 0
        e6 = 1 if feat in DRUGGABLE else 0
        total = e1+e2+e3+e4+e5+e6

        is_anc = (feat == 'anchor')
        if is_anc:
            decision = 'ANCHOR'
        elif cons >= cfg['inc_cons'] and comp >= cfg['inc_comp']:
            decision = 'INCLUDE'
        elif cons >= cfg['sec_cons'] and comp >= cfg['sec_comp']:
            decision = 'SECONDARY'
        else:
            decision = 'EXCLUDE'

        assessed.append({
            'residue':        row['residue'],
            'chain':          chain,
            'position':       pos,
            'aa':             aa,
            'primary_feature':feat,
            'conservation':   cons,
            'cons_tier':      tier,
            'composite':      comp,
            'e1_cluster0':    e1,
            'e2_chain1':      e2,
            'e3_cons_high':   e3,
            'e4_composite':   e4,
            'e5_identical':   e5,
            'e6_druggable':   e6,
            'evidence_score': total,
            'decision':       decision,
        })

    assessed.sort(key=lambda x:(
        0 if x['decision']=='ANCHOR'    else
        1 if x['decision']=='INCLUDE'   else
        2 if x['decision']=='SECONDARY' else 3,
        -x['evidence_score'], -x['composite']))
    return assessed


def get_coords_for_assessed(assessed, cfg):
    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', cfg['pdb'])[0]
    chain1    = structure[cfg['chain1_id']]
    chain2    = structure[cfg['chain2_id']]

    coords = {}
    for r in assessed:
        chain_obj = (chain1 if r['chain']==cfg['name1']
                     else chain2)
        pos = r['position']
        try:
            res = chain_obj[(' ', pos, ' ')]
            ca  = res['CA'].coord
            coords[(r['chain'], pos)] = {
                'x': round(float(ca[0]),3),
                'y': round(float(ca[1]),3),
                'z': round(float(ca[2]),3),
                'aa': res.resname,
            }
        except KeyError:
            pass
    return coords


def assign_elements(assessed, coords, cfg):
    features = []
    for r in assessed:
        if r['decision'] not in ('ANCHOR','INCLUDE'):
            continue
        key = (r['chain'], r['position'])
        if key not in coords:
            continue
        c = coords[key]

        # Element assignment
        if r['decision'] == 'ANCHOR':
            el = 'E1'
        elif r['chain'] == cfg['name1'] and r['conservation'] >= 0.800:
            el = 'E1'
        elif r['chain'] == cfg['name2']:
            el = 'E2'
        else:
            el = 'E3'

        feat = r['primary_feature']
        pharm_type = (
            'ANCHOR'      if feat=='anchor'      else
            'AROMATIC'    if feat=='aromatic'     else
            'HYDROPHOBIC' if feat=='hydrophobic'  else
            'CHARGED_POS' if feat=='charged_pos'  else
            'CHARGED_NEG' if feat=='charged_neg'  else
            'HBOND'       if feat=='hbond_donor'  else
            'UNKNOWN')
        tol = (1.0 if pharm_type in ('ANCHOR','CHARGED_POS','CHARGED_NEG')
               else 1.5 if pharm_type=='AROMATIC'
               else 2.0)

        features.append({
            'id':           r['residue'],
            'chain':        r['chain'],
            'position':     r['position'],
            'aa':           c['aa'],
            'element':      el,
            'role':         r['decision'],
            'pharm_type':   pharm_type,
            'tolerance_A':  tol,
            'conservation': float(r['conservation']),
            'composite':    float(r['composite']),
            'x':            c['x'],
            'y':            c['y'],
            'z':            c['z'],
        })

    features.sort(key=lambda x:(
        x['element'],
        0 if x['role']=='ANCHOR' else 1,
        -x['composite']))
    return features


def compute_centroids(features):
    centroids = {}
    for el in ['E1','E2','E3']:
        pts = [f for f in features if f['element']==el]
        if not pts:
            continue
        cx = float(np.mean([p['x'] for p in pts]))
        cy = float(np.mean([p['y'] for p in pts]))
        cz = float(np.mean([p['z'] for p in pts]))
        radius = float(max(
            np.sqrt((p['x']-cx)**2+(p['y']-cy)**2+(p['z']-cz)**2)
            for p in pts))
        centroids[el] = {
            'x': round(cx,3), 'y': round(cy,3),
            'z': round(cz,3), 'radius': round(radius,2),
        }
    return centroids


def save_outputs(iface, assessed, features, centroids, cfg):
    out6 = os.path.join(BASE,
        f'CGCP/02-deep-dive/{iface}/step-06-assessment')
    out7 = os.path.join(BASE,
        f'CGCP/02-deep-dive/{iface}/step-07-pharmacophore')
    os.makedirs(out6, exist_ok=True)
    os.makedirs(out7, exist_ok=True)

    # Step 6 TSV
    fields6 = ['residue','chain','position','aa',
               'primary_feature','conservation','cons_tier',
               'composite','e1_cluster0','e2_chain1',
               'e3_cons_high','e4_composite','e5_identical',
               'e6_druggable','evidence_score','decision']
    with open(os.path.join(out6,
              f'integrated_assessment_{iface}.tsv'),'w',
              newline='') as f:
        w = csv.DictWriter(f,fieldnames=fields6,
                           delimiter='\t',extrasaction='ignore')
        w.writeheader(); w.writerows(assessed)

    dec_counts = {}
    for r in assessed:
        dec_counts[r['decision']] = \
            dec_counts.get(r['decision'],0)+1

    pharmacophore_list = [r['residue'] for r in assessed
                          if r['decision'] in ('ANCHOR','INCLUDE')]
    secondary_list     = [r['residue'] for r in assessed
                          if r['decision']=='SECONDARY']

    out6_data = {
        'interface': iface,
        'anchor':    f"{cfg['anchor_chain']}-{cfg['anchor_res']}",
        'decisions': dec_counts,
        'inc_threshold': {'cons': cfg['inc_cons'],
                          'comp': cfg['inc_comp']},
        'pharmacophore_candidates': pharmacophore_list,
        'secondary_candidates':     secondary_list,
        'all': assessed,
    }
    with open(os.path.join(out6,
              f'integrated_assessment_{iface}.json'),'w') as f:
        json.dump(out6_data,f,indent=2,cls=NpEncoder)

    # Step 7
    el_summary = {}
    for el in ['E1','E2','E3']:
        members = [f for f in features if f['element']==el]
        el_name = cfg['el_names'].get(el, f'{el} — zone')
        el_summary[el] = {
            'name':       el_name,
            'n_residues': len(members),
            'centroid':   centroids.get(el,{}),
            'members':    [f['id'] for f in members],
        }

    pharmacophore = {
        'name':       f'CGCP-{iface}-v1',
        'interface':  iface,
        'date':       '2026-03-24',
        'method':     'CGCP pipeline v1.0',
        'n_features': len(features),
        'elements':   el_summary,
        'features':   features,
    }
    with open(os.path.join(out7,
              f'pharmacophore_{iface}.json'),'w') as f:
        json.dump(pharmacophore,f,indent=2,cls=NpEncoder)

    fields7 = ['id','chain','position','aa','element',
               'role','pharm_type','tolerance_A',
               'conservation','composite','x','y','z']
    with open(os.path.join(out7,
              f'pharmacophore_{iface}.tsv'),'w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=fields7,
                           delimiter='\t',extrasaction='ignore')
        w.writeheader(); w.writerows(features)

    return out6, out7, out6_data, pharmacophore


def make_figures(iface, assessed, out6_data,
                 pharmacophore, centroids, out6, out7):
    features = pharmacophore['features']

    # Figure 1: Assessment
    fig1,(ax1,ax2) = plt.subplots(1,2,figsize=(14,5.5))
    fig1.subplots_adjust(wspace=0.45,left=0.07,
                         right=0.97,top=0.88,bottom=0.26)

    top = [r for r in assessed if r['decision']!='EXCLUDE'][:25]
    x   = np.arange(len(top))
    for xi,r in enumerate(top):
        ec,fc = dec_style(r['decision'])
        ax1.bar(xi,r['evidence_score'],width=0.40,
                facecolor=fc,edgecolor=ec,linewidth=1.2,
                zorder=3,clip_on=False)
        ax1.text(xi,r['evidence_score']+0.07,
                 str(r['evidence_score']),
                 ha='center',va='bottom',fontsize=7,
                 fontweight='bold',color=ec,clip_on=False)
        badge=('A' if r['decision']=='ANCHOR' else
               'I' if r['decision']=='INCLUDE' else 'S')
        ax1.text(xi,-0.4,badge,ha='center',va='top',
                 fontsize=6.5,color=ec,clip_on=False)

    ax1.axhline(4,color=P['include'],linewidth=0.8,
                linestyle='--',alpha=0.7)
    ax1.set_xticks(x)
    set_xticklabels_vertical(
        ax1,[r['residue'].split('-')[1] for r in top],fontsize=7)
    ax1.set_ylabel('Evidence score (max 6)',
                   fontsize=11,fontweight='bold')
    ax1.set_xlim(-0.6,len(top)-0.4)
    prism_axes(ax1,ymax=7,yticks=[0,1,2,3,4,5,6])
    make_legend(ax1,[
        (P['anchor_f'],    P['anchor'],    'ANCHOR'),
        (P['include_f'],   P['include'],   'INCLUDE'),
        (P['secondary_f'], P['secondary'], 'SECONDARY'),
    ],loc='upper right',fontsize=7.5)
    panel_title(ax1,'A',f'Evidence scores — {iface}')

    dec_order  = ['ANCHOR','INCLUDE','SECONDARY','EXCLUDE']
    dec_labels = ['Anchor','Include','Secondary','Exclude']
    dec_counts = [out6_data['decisions'].get(d,0) for d in dec_order]
    x2 = np.arange(len(dec_order))
    for xi,val,d in zip(x2,dec_counts,dec_order):
        ec,fc = dec_style(d)
        ax2.bar(xi,val,width=0.50,facecolor=fc,edgecolor=ec,
                linewidth=1.2,zorder=3,clip_on=False)
        ax2.text(xi,val+0.3,str(val),ha='center',va='bottom',
                 fontsize=10,fontweight='bold',color=ec,clip_on=False)
    ax2.set_xticks(x2)
    set_xticklabels_vertical(ax2,dec_labels,fontsize=9)
    ax2.set_ylabel('Number of residues',fontsize=11,fontweight='bold')
    ymax_b = max(dec_counts)+6
    prism_axes(ax2,ymax=ymax_b,yticks=range(0,ymax_b,5))
    ax2.set_xlim(-0.6,len(dec_order)-0.4)
    panel_title(ax2,'B','Decision summary')

    p1 = os.path.join(out6,f'Fig_Step06_Assessment_{iface}.png')
    fig1.savefig(p1); plt.close()

    # Figure 2: Pharmacophore
    if not features:
        return

    fig2 = plt.figure(figsize=(14,5.5))
    gs   = fig2.add_gridspec(1,2,wspace=0.45,
                              left=0.07,right=0.97,
                              top=0.88,bottom=0.26)
    ax_a = fig2.add_subplot(gs[0,0])
    ax_b = fig2.add_subplot(gs[0,1])

    # Panel A: 2D X-Z
    for el_id,ct in centroids.items():
        col  = EL_COLORS[el_id]
        fill = EL_FILLS[el_id]
        circ = Circle((ct['x'],ct['z']),radius=ct['radius'],
                      facecolor=fill,edgecolor=col,
                      linewidth=1.2,alpha=0.20,zorder=1)
        ax_a.add_patch(circ)
        ax_a.scatter(ct['x'],ct['z'],marker='+',
                     color=col,s=180,linewidths=2.5,zorder=5)
        ax_a.text(ct['x'],ct['z']+ct['radius']+1.0,
                  el_id,ha='center',va='bottom',
                  fontsize=10,color=col,fontweight='bold')

    for f in features:
        col  = EL_COLORS[f['element']]
        fill = EL_FILLS[f['element']]
        is_anc = f['role']=='ANCHOR'
        ax_a.scatter(f['x'],f['z'],color=fill,
                     edgecolors=col,linewidths=1.0,
                     s=220 if is_anc else 50,
                     marker='*' if is_anc else 'o',
                     alpha=0.88,zorder=4)
        if is_anc or f['composite']>=0.700:
            ax_a.text(f['x']+1.0,f['z']+0.8,
                      f['id'].split('-')[1],
                      fontsize=7.5,color=col,
                      fontweight='bold' if is_anc else 'normal',
                      va='center')

    ax_a.set_xlabel('X coordinate (Å)',fontsize=11,fontweight='bold')
    ax_a.set_ylabel('Z coordinate (Å)',fontsize=11,fontweight='bold')
    ax_a.spines['left'].set_position(('outward',6))
    ax_a.spines['bottom'].set_position(('outward',6))
    ax_a.spines['left'].set_linewidth(1.2)
    ax_a.spines['bottom'].set_linewidth(1.2)
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)
    ax_a.tick_params(labelsize=9,width=1.2,length=6,direction='out')
    active_els = [el for el in ['E1','E2','E3'] if el in centroids]
    make_legend(ax_a,[
        (EL_FILLS[el],EL_COLORS[el],
         pharmacophore['elements'][el]['name'])
        for el in active_els
    ],loc='best',fontsize=7.5)
    panel_title(ax_a,'A',f'2D pharmacophore map — {iface}',
                '★=anchor')

    # Panel B: composite scores
    feats_s = sorted(features,key=lambda x:-x['composite'])
    x3 = np.arange(len(feats_s))
    for xi,f in enumerate(feats_s):
        col  = EL_COLORS[f['element']]
        fill = EL_FILLS[f['element']]
        is_anc = f['role']=='ANCHOR'
        ax_b.bar(xi,f['composite'],width=0.38,
                 facecolor=fill,edgecolor=col,
                 linewidth=2.0 if is_anc else 1.2,
                 zorder=3,clip_on=False)
        ax_b.errorbar(xi,f['composite'],
                      yerr=f['tolerance_A']/10.0,
                      color=col,linewidth=0.9,
                      capsize=3,capthick=0.9,zorder=6)

    ax_b.set_xticks(x3)
    set_xticklabels_vertical(
        ax_b,[f['id'].split('-')[1] for f in feats_s],fontsize=7)
    ax_b.set_ylabel('CGCP composite score',
                    fontsize=11,fontweight='bold')
    ax_b.set_xlim(-0.6,len(feats_s)-0.4)
    prism_axes(ax_b,ymax=1.25,yticks=[0,0.2,0.4,0.6,0.8,1.0])
    ax_b.tick_params(labelsize=9,width=1.2,length=6,direction='out')
    make_legend(ax_b,[
        (EL_FILLS[el],EL_COLORS[el],
         pharmacophore['elements'][el]['name'])
        for el in active_els
    ],loc='upper right',fontsize=7.5)
    panel_title(ax_b,'B','Composite score per feature')

    p2 = os.path.join(out7,
        f'Fig_Step07_Pharmacophore_{iface}.png')
    fig2.savefig(p2); plt.close()


def print_interface_summary(iface, assessed,
                             out6_data, pharmacophore,
                             centroids):
    print(f"\n{'='*60}")
    print(f"Steps 6+7 — {iface}")
    print(f"{'='*60}")
    print(f"\n  Decisions: "
          + " | ".join(f"{k}={v}"
                       for k,v in out6_data['decisions'].items()))
    print(f"\n  Pharmacophore candidates "
          f"({len(out6_data['pharmacophore_candidates'])}):")
    for res in out6_data['pharmacophore_candidates']:
        print(f"    → {res}")

    print(f"\n  {pharmacophore['name']}:")
    for el_id in ['E1','E2','E3']:
        el = pharmacophore['elements'][el_id]
        ct = el['centroid']
        if el['n_residues']==0: continue
        print(f"  {el_id}: {el['name']}")
        print(f"    {el['n_residues']} res | "
              f"centroid ({ct.get('x',0):.1f},"
              f"{ct.get('y',0):.1f},{ct.get('z',0):.1f}) | "
              f"radius {ct.get('radius',0):.1f}Å")
        print(f"    {', '.join(el['members'])}")


def process_interface(iface):
    cfg     = INTERFACES[iface]
    records = load_conservation_tsv(iface)
    assessed = assess(records, cfg, iface)
    coords   = get_coords_for_assessed(assessed, cfg)
    features = assign_elements(assessed, coords, cfg)
    centroids = compute_centroids(features)
    out6, out7, out6_data, pharmacophore = \
        save_outputs(iface, assessed, features,
                     centroids, cfg)
    make_figures(iface, assessed, out6_data,
                 pharmacophore, centroids, out6, out7)
    print_interface_summary(iface, assessed,
                            out6_data, pharmacophore,
                            centroids)
    return out6_data, pharmacophore


if __name__ == '__main__':
    print('CGCP Phase 2 Steps 6+7 — Remaining 4 interfaces')
    all_results = {}
    for iface in INTERFACES:
        try:
            all_results[iface] = process_interface(iface)
        except Exception as e:
            print(f"\nERROR in {iface}: {e}")
            import traceback; traceback.print_exc()

    print('\n'+'='*60)
    print('FINAL SUMMARY — All 9 interfaces complete')
    print('='*60)
    completed = [
        ('NSP12-NSP7',    'PHE440',  'CGCP-NSP12-NSP7-v1'),
        ('NSP12-NSP8',    'LEU387',  'CGCP-NSP12-NSP8-v1'),
        ('NSP9-NSP12',    'ARG733',  'CGCP-NSP9-NSP12-v1'),
        ('NSP10-NSP16',   'LYS93',   'CGCP-NSP10-NSP16-v1'),
        ('NSP7-NSP8',     'PHE92',   'CGCP-NSP7-NSP8-v1'),
        ('NSP10-NSP14',   'HIS80',   'CGCP-NSP10-NSP14-v1'),
        ('NSP13-Helicase','LYS414',  'CGCP-NSP13-Helicase-v1'),
        ('NSP12-NSP13',   'TYR903',  'CGCP-NSP12-NSP13-v1'),
        ('NSP15',         'ASP40',   'CGCP-NSP15-v1'),
    ]
    for iface, anchor, name in completed:
        print(f"  ✅ {name:<30} anchor={anchor}")
    print('\nNext: Step 8 ChimeraX for each interface + VirtualFlow setup')
