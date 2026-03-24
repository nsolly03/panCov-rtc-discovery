#!/usr/bin/env python3
"""
CGCP Phase 2 Steps 4+5 — Four remaining interfaces
DBSCAN clustering + conservation overlay for:
  NSP10-NSP14, NSP13-Helicase, NSP12-NSP13, NSP15

Same method as all previous interfaces.

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step04_05_remaining_interfaces.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
from sklearn.cluster import DBSCAN
from prism_style import (apply_prism, prism_axes,
                          set_xticklabels_vertical,
                          make_legend, panel_title,
                          COLORS, CLUSTER_COLORS, CLUSTER_FILLS)

apply_prism()

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):  return int(obj)
        if isinstance(obj, np.floating): return float(obj)
        if isinstance(obj, np.ndarray):  return obj.tolist()
        if isinstance(obj, np.bool_):    return bool(obj)
        return super().default(obj)

BASE = os.path.expanduser('~/projects/rtc-pan-coronavirus')

EPS         = 8.0
MIN_SAMPLES = 2
NOISE_COL   = COLORS['gray']
NOISE_FILL  = COLORS['gray_fill']

TIER_COLORS = {
    'identical': '#D7191C',
    'high':      '#F46D43',
    'moderate':  '#FEE090',
    'variable':  '#AAAAAA',
}
TIER_FILLS = {
    'identical': '#FDAE61',
    'high':      '#FDD0A2',
    'moderate':  '#FFFFBF',
    'variable':  '#DDDDDD',
}

INTERFACES = [
    'NSP10-NSP14',
    'NSP13-Helicase',
    'NSP12-NSP13',
    'NSP15',
]


def get_tier(cons):
    if cons == 1.000:  return 'identical'
    if cons >= 0.800:  return 'high'
    if cons >= 0.600:  return 'moderate'
    return 'variable'


def load_features(iface):
    path = os.path.join(BASE,
        f'CGCP/02-deep-dive/{iface}/step-03-features/'
        f'feature_classification_{iface}.tsv')
    records = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['ca_x'] and row['ca_y'] and row['ca_z']:
                records.append({
                    'residue':         row['residue'],
                    'chain':           row['chain'],
                    'position':        int(row['position']),
                    'aa':              row['aa'],
                    'primary_feature': row['primary_feature'],
                    'contact_score':   int(row['contact_score']),
                    'conservation':    float(row['conservation']),
                    'composite':       float(row['composite']),
                    'is_anchor':       int(row['is_anchor']),
                    'x': float(row['ca_x']),
                    'y': float(row['ca_y']),
                    'z': float(row['ca_z']),
                })
    return records


def run_dbscan(records):
    coords = np.array([[r['x'],r['y'],r['z']] for r in records])
    db = DBSCAN(eps=EPS, min_samples=MIN_SAMPLES,
                metric='euclidean').fit(coords)
    return [int(l) for l in db.labels_]


def analyze_clusters(records, labels):
    cluster_info = {}
    for label in sorted(set(labels)):
        members = [records[i] for i,l
                   in enumerate(labels) if l==label]
        if not members: continue
        coords    = np.array([[r['x'],r['y'],r['z']]
                               for r in members])
        centroid  = coords.mean(axis=0)
        radius    = float(max(
            np.linalg.norm(c-centroid) for c in coords))
        mean_cons = float(np.mean([r['conservation']
                                   for r in members]))
        mean_comp = float(np.mean([r['composite']
                                   for r in members]))
        has_anc   = bool(any(r['primary_feature']=='anchor'
                             for r in members))

        if   has_anc and mean_cons >= 0.700:
            decision = 'PRIMARY PHARMACOPHORE'
        elif mean_cons >= 0.800 and mean_comp >= 0.500:
            decision = 'SECONDARY PHARMACOPHORE'
        elif mean_cons >= 0.600:
            decision = 'SUPPORTING'
        else:
            decision = 'DEPRIORITIZE'

        cluster_info[label] = {
            'label':      int(label),
            'n_residues': int(len(members)),
            'centroid_x': round(float(centroid[0]),3),
            'centroid_y': round(float(centroid[1]),3),
            'centroid_z': round(float(centroid[2]),3),
            'radius':     round(radius,2),
            'mean_cons':  round(mean_cons,3),
            'mean_comp':  round(mean_comp,3),
            'has_anchor': has_anc,
            'decision':   decision,
            'members':    [r['residue'] for r in members],
        }
    return cluster_info


def process_interface(iface):
    print(f"\n{'='*55}")
    print(f"Steps 4+5: {iface}")
    print(f"{'='*55}")

    out4 = os.path.join(BASE,
        f'CGCP/02-deep-dive/{iface}/step-04-clusters')
    out5 = os.path.join(BASE,
        f'CGCP/02-deep-dive/{iface}/step-05-conservation')
    os.makedirs(out4, exist_ok=True)
    os.makedirs(out5, exist_ok=True)

    records      = load_features(iface)
    labels       = run_dbscan(records)
    cluster_info = analyze_clusters(records, labels)

    # Save Step 4
    rows = [{
        'residue':         r['residue'],
        'chain':           r['chain'],
        'position':        int(r['position']),
        'aa':              r['aa'],
        'primary_feature': r['primary_feature'],
        'conservation':    float(r['conservation']),
        'composite':       float(r['composite']),
        'cluster':         int(labels[i]),
        'x':               float(r['x']),
        'y':               float(r['y']),
        'z':               float(r['z']),
    } for i,r in enumerate(records)]

    tsv4 = os.path.join(out4, f'clusters_{iface}.tsv')
    with open(tsv4,'w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=list(rows[0].keys()),
                           delimiter='\t')
        w.writeheader(); w.writerows(rows)

    json4 = {
        'interface':   iface,
        'dbscan_eps':  float(EPS),
        'min_samples': int(MIN_SAMPLES),
        'n_clusters':  int(len([l for l in cluster_info if l>=0])),
        'n_noise':     int(sum(1 for l in labels if l==-1)),
        'clusters':    {str(k):v for k,v in cluster_info.items()},
        'residues':    rows,
    }
    with open(os.path.join(out4,f'clusters_{iface}.json'),'w') as f:
        json.dump(json4,f,indent=2,cls=NpEncoder)
    print(f"  Step4 TSV: {tsv4}")

    # Save Step 5
    result = [{**r,'cons_tier':get_tier(r['conservation'])}
              for r in records]
    result.sort(key=lambda x:(
        0 if x['primary_feature']=='anchor' else 1,
        -x['conservation'],-x['composite']))

    rows5 = [{**r,'cluster':int(labels[i])}
             for i,r in enumerate(records)]
    rows5.sort(key=lambda x:(
        0 if x['primary_feature']=='anchor' else 1,
        -x['conservation'],-x['composite']))

    fields5 = ['residue','chain','position','aa',
               'primary_feature','cluster',
               'conservation','cons_tier','composite',
               'x','y','z']
    tsv5 = os.path.join(out5,
        f'conservation_overlay_{iface}.tsv')
    with open(tsv5,'w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=fields5,
                           delimiter='\t',extrasaction='ignore')
        w.writeheader(); w.writerows(rows5)

    tier_stats = {}
    for tier in ['identical','high','moderate','variable']:
        members = [r for r in result if r['cons_tier']==tier]
        if members:
            tier_stats[tier] = {
                'n':         len(members),
                'mean_comp': round(float(np.mean(
                    [r['composite'] for r in members])),3),
                'members':   [r['residue'] for r in members],
            }

    with open(os.path.join(out5,
              f'conservation_overlay_{iface}.json'),'w') as f:
        json.dump({'interface':iface,'tier_stats':tier_stats,
                   'records':rows5},f,indent=2,cls=NpEncoder)
    print(f"  Step5 TSV: {tsv5}")

    # ── Figures ──────────────────────────────────────────────
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(13,5.5))
    fig.subplots_adjust(wspace=0.45,left=0.07,
                        right=0.97,top=0.88,bottom=0.22)

    # Panel A: DBSCAN X-Z
    plotted = set()
    for i,r in enumerate(records):
        lbl = labels[i]
        if lbl==-1:
            col,fill = NOISE_COL,NOISE_FILL
            lbl_label = ('Noise' if 'noise' not in plotted
                         else None)
            plotted.add('noise')
        else:
            idx  = lbl%len(CLUSTER_COLORS)
            col  = CLUSTER_COLORS[idx]
            fill = CLUSTER_FILLS[idx]
            ck   = f'C{lbl}'
            dec  = cluster_info[lbl]['decision']
            lbl_label = (f"C{lbl}: {dec[:16]}"
                         if ck not in plotted else None)
            plotted.add(ck)

        is_anc = r['primary_feature']=='anchor'
        ax1.scatter(r['x'],r['z'],color=fill,
                    edgecolors=col,linewidths=1.0,
                    s=220 if is_anc else 45,
                    marker='*' if is_anc else 'o',
                    alpha=0.85,zorder=4,label=lbl_label)
        if is_anc or r['composite']>=0.800:
            ax1.annotate(r['residue'].split('-')[1],
                         (r['x'],r['z']),fontsize=7.5,
                         color=col,
                         fontweight='bold' if is_anc else 'normal',
                         xytext=(4,3),textcoords='offset points')

    for lbl,info in cluster_info.items():
        if lbl<0: continue
        idx = lbl%len(CLUSTER_COLORS)
        ax1.scatter(info['centroid_x'],info['centroid_z'],
                    marker='+',color=CLUSTER_COLORS[idx],
                    s=180,linewidths=2.5,zorder=6)

    ax1.set_xlabel('X coordinate (Å)',fontsize=11,fontweight='bold')
    ax1.set_ylabel('Z coordinate (Å)',fontsize=11,fontweight='bold')
    ax1.spines['left'].set_position(('outward',6))
    ax1.spines['bottom'].set_position(('outward',6))
    ax1.spines['left'].set_linewidth(1.2)
    ax1.spines['bottom'].set_linewidth(1.2)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.tick_params(labelsize=9,width=1.2,length=6,direction='out')
    ax1.legend(fontsize=7,frameon=True,facecolor='white',
               edgecolor='#CCCCCC',loc='best',framealpha=0.9)
    panel_title(ax1,'A',f'DBSCAN clusters — {iface}',
                f'eps={EPS}Å | ★=anchor')

    # Panel B: conservation tiers (bar chart)
    tier_order  = ['identical','high','moderate','variable']
    tier_labels = ['Identical','High','Moderate','Variable']
    tier_counts = [tier_stats.get(t,{}).get('n',0)
                   for t in tier_order]
    x = np.arange(len(tier_order))
    for xi,val,t in zip(x,tier_counts,tier_order):
        ec = TIER_COLORS[t]; fc = TIER_FILLS[t]
        ax2.bar(xi,val,width=0.50,facecolor=fc,edgecolor=ec,
                linewidth=1.2,zorder=3,clip_on=False)
        ax2.text(xi,val+0.3,str(val),
                 ha='center',va='bottom',fontsize=10,
                 fontweight='bold',color=ec,clip_on=False)

    # Overlay mean_comp as line
    ax2b = ax2.twinx()
    comp_vals = [tier_stats.get(t,{}).get('mean_comp',0)
                 for t in tier_order]
    ax2b.plot(x,comp_vals,'o--',color=COLORS['blue'],
              linewidth=1.5,markersize=6,zorder=5)
    ax2b.set_ylabel('Mean composite',fontsize=9,
                    color=COLORS['blue'])
    ax2b.tick_params(axis='y',labelsize=8,
                     colors=COLORS['blue'])
    ax2b.set_ylim(0,1.2)
    ax2b.spines['right'].set_linewidth(1.0)

    ax2.set_xticks(x)
    set_xticklabels_vertical(ax2,tier_labels,fontsize=9)
    ax2.set_ylabel('Number of residues',
                   fontsize=11,fontweight='bold')
    ymax = max(tier_counts)+5 if tier_counts else 10
    prism_axes(ax2,ymax=ymax,yticks=range(0,ymax+1,5))
    ax2.set_xlim(-0.6,len(tier_order)-0.4)
    ax2.tick_params(labelsize=9,width=1.2,length=6,direction='out')
    panel_title(ax2,'B','Conservation tiers + mean composite')

    path = os.path.join(out4,
        f'Fig_Steps04-05_{iface}.png')
    fig.savefig(path); plt.close()
    print(f"  Figure: {path}")

    # Print summary
    n_clusters = len([l for l in cluster_info if l>=0])
    n_noise    = sum(1 for l in labels if l==-1)
    print(f"  DBSCAN: {n_clusters} clusters | {n_noise} noise")
    for lbl in sorted(cluster_info.keys()):
        info = cluster_info[lbl]
        name = f"C{lbl}" if lbl>=0 else "Noise"
        print(f"  {name}: {info['n_residues']} res | "
              f"cons={info['mean_cons']:.3f} | "
              f"{info['decision']}")

    print(f"  Conservation tiers:")
    for tier in ['identical','high','moderate','variable']:
        if tier in tier_stats:
            ts = tier_stats[tier]
            print(f"    {tier:<12} n={ts['n']:>3} "
                  f"mean_comp={ts['mean_comp']:.3f}")

    # Pharmacophore preview
    candidates = [r for r in result
                  if r['conservation']>=0.800
                  and r['composite']>=0.500]
    candidates.sort(key=lambda x:-x['composite'])
    print(f"  Pharmacophore preview (cons>=0.800, comp>=0.500): "
          f"{len(candidates)}")
    for r in candidates[:8]:
        sym = '★' if r['primary_feature']=='anchor' else ' '
        print(f"    {sym} {r['residue']:<24} "
              f"cons={r['conservation']:.3f} "
              f"comp={r['composite']:.3f} "
              f"{r['cons_tier']}")

    return {
        'n_clusters':  n_clusters,
        'n_noise':     n_noise,
        'tier_stats':  tier_stats,
        'n_candidates': len(candidates),
        'top_candidates': [r['residue'] for r in candidates[:5]],
        'cluster_info': cluster_info,
    }


if __name__ == '__main__':
    print('CGCP Phase 2 Steps 4+5 — Remaining 4 interfaces')
    results = {}
    for iface in INTERFACES:
        try:
            results[iface] = process_interface(iface)
        except Exception as e:
            print(f"  ERROR in {iface}: {e}")
            import traceback; traceback.print_exc()

    print('\n'+'='*55)
    print('SUMMARY — Steps 4+5 all interfaces')
    print('='*55)
    for iface, info in results.items():
        print(f"\n{iface}:")
        print(f"  Clusters={info['n_clusters']} | "
              f"Noise={info['n_noise']}")
        print(f"  Candidates: {info['n_candidates']}")
        print(f"  Top: {info['top_candidates'][:3]}")
    print('\nNext: Steps 6+7 for each interface')
