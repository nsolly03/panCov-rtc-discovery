#!/usr/bin/env python3
"""
CGCP Phase 2 Step 7 - Pharmacophore Hypothesis: NSP12-NSP7
Formalizes the 8 pharmacophore candidates into named features
with 3D coordinates, tolerances, and interaction types.

Pharmacophore elements:
  Element 1 - Aromatic hydrophobic core (6 residues)
    PHE440 (anchor), PHE442, PHE415, TYR420, PRO412, GLY413
    Feature type: AROMATIC + HYDROPHOBIC
    Interaction: pi-stacking, van der Waals

  Element 2 - Electrostatic anchor (1 residue)
    GLU431
    Feature type: NEGATIVE CHARGE
    Interaction: salt bridge with LYS1(NSP7)

  Element 3 - Distal aromatic (1 residue)
    PHE843
    Feature type: AROMATIC
    Interaction: hydrophobic contact

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step07_pharmacophore_hypothesis_NSP12-NSP7.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
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
    'axes.titlepad':      8,
    'xtick.direction':    'out',
    'ytick.direction':    'out',
    'xtick.major.size':   4,
    'ytick.major.size':   4,
    'xtick.major.width':  0.75,
    'ytick.major.width':  0.75,
    'xtick.major.pad':    4,
    'ytick.major.pad':    4,
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

# Pharmacophore feature colors
FEAT_COLORS = {
    'AROMATIC':      '#E6550D',
    'HYDROPHOBIC':   '#636363',
    'NEG_CHARGE':    '#4DAC26',
    'ANCHOR':        '#D7191C',
}
FEAT_FILLS = {k: v+'55' for k, v in FEAT_COLORS.items()}

# ── Pharmacophore definition ─────────────────────────────────
# Based on Step 6 decisions — all cons=1.000
PHARMACOPHORE_DEF = {
    'PHE440': {
        'element':      'E1',
        'role':         'ANCHOR',
        'pharm_type':   'AROMATIC',
        'tolerance_A':  1.5,
        'interaction':  'pi-stacking + vdW core',
        'cons':         1.000,
        'composite':    1.000,
    },
    'PHE442': {
        'element':      'E1',
        'role':         'INCLUDE',
        'pharm_type':   'AROMATIC',
        'tolerance_A':  1.5,
        'interaction':  'pi-stacking',
        'cons':         1.000,
        'composite':    0.774,
    },
    'TYR420': {
        'element':      'E1',
        'role':         'INCLUDE',
        'pharm_type':   'AROMATIC',
        'tolerance_A':  1.5,
        'interaction':  'pi-stacking + H-bond donor',
        'cons':         1.000,
        'composite':    0.616,
    },
    'PHE415': {
        'element':      'E1',
        'role':         'INCLUDE',
        'pharm_type':   'AROMATIC',
        'tolerance_A':  1.5,
        'interaction':  'pi-stacking',
        'cons':         1.000,
        'composite':    0.527,
    },
    'PRO412': {
        'element':      'E1',
        'role':         'INCLUDE',
        'pharm_type':   'HYDROPHOBIC',
        'tolerance_A':  2.0,
        'interaction':  'vdW hydrophobic',
        'cons':         1.000,
        'composite':    0.895,
    },
    'GLY413': {
        'element':      'E1',
        'role':         'INCLUDE',
        'pharm_type':   'HYDROPHOBIC',
        'tolerance_A':  2.0,
        'interaction':  'vdW backbone',
        'cons':         1.000,
        'composite':    0.454,
    },
    'GLU431': {
        'element':      'E2',
        'role':         'INCLUDE',
        'pharm_type':   'NEG_CHARGE',
        'tolerance_A':  1.0,
        'interaction':  'salt bridge with LYS1(NSP7)',
        'cons':         1.000,
        'composite':    0.748,
    },
    'PHE843': {
        'element':      'E3',
        'role':         'INCLUDE',
        'pharm_type':   'AROMATIC',
        'tolerance_A':  2.0,
        'interaction':  'vdW hydrophobic distal',
        'cons':         1.000,
        'composite':    0.488,
    },
}

BASE     = os.path.expanduser('~/projects/rtc-pan-coronavirus')
PDB_FILE = os.path.join(BASE,
    '03-virtual-screening/NSP12-NSP7_3/receptor_NSP12-NSP7_3.pdb')
S6_DIR   = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP7/step-06-assessment')
OUT_DIR  = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP7/step-07-pharmacophore')
os.makedirs(OUT_DIR, exist_ok=True)


# ── Get Ca coordinates ───────────────────────────────────────
def get_coords():
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', PDB_FILE)[0]
    chain_a = structure['A']

    coords = {}
    res_num = {
        'PHE440': 440, 'PHE442': 442, 'TYR420': 420,
        'PHE415': 415, 'PRO412': 412, 'GLY413': 413,
        'GLU431': 431, 'PHE843': 843,
    }

    for name, pos in res_num.items():
        try:
            res  = chain_a[(' ', pos, ' ')]
            ca   = res['CA'].coord
            coords[name] = {
                'x': round(float(ca[0]), 3),
                'y': round(float(ca[1]), 3),
                'z': round(float(ca[2]), 3),
            }
        except KeyError:
            print(f"  Warning: {name} not found")

    return coords


# ── Compute element centroids ────────────────────────────────
def compute_centroids(coords):
    elements = {}
    for name, defn in PHARMACOPHORE_DEF.items():
        el = defn['element']
        if el not in elements:
            elements[el] = []
        if name in coords:
            elements[el].append(coords[name])

    centroids = {}
    for el, pts in elements.items():
        cx = np.mean([p['x'] for p in pts])
        cy = np.mean([p['y'] for p in pts])
        cz = np.mean([p['z'] for p in pts])
        # Radius = max distance from centroid to any member
        radius = max(
            np.sqrt((p['x']-cx)**2 +
                    (p['y']-cy)**2 +
                    (p['z']-cz)**2)
            for p in pts
        )
        centroids[el] = {
            'x':      round(float(cx), 3),
            'y':      round(float(cy), 3),
            'z':      round(float(cz), 3),
            'radius': round(float(radius), 2),
        }

    return centroids


# ── Build pharmacophore JSON ─────────────────────────────────
def build_pharmacophore(coords, centroids):
    features = []
    for name, defn in PHARMACOPHORE_DEF.items():
        if name not in coords:
            continue
        feat = {
            'id':           name,
            'residue':      f'NSP12-{name}',
            'element':      defn['element'],
            'role':         defn['role'],
            'pharm_type':   defn['pharm_type'],
            'tolerance_A':  defn['tolerance_A'],
            'interaction':  defn['interaction'],
            'conservation': defn['cons'],
            'composite':    defn['composite'],
            'x':            coords[name]['x'],
            'y':            coords[name]['y'],
            'z':            coords[name]['z'],
        }
        features.append(feat)

    # Sort by element then composite
    features.sort(key=lambda x: (x['element'],
                                  -x['composite']))

    pharmacophore = {
        'name':       'CGCP-NSP12-NSP7-v1',
        'interface':  'NSP12-NSP7',
        'date':       '2026-03-18',
        'method':     'CGCP pipeline v1.0',
        'reference_pdb': '7BV2',
        'chain':      'A (NSP12)',
        'n_features': len(features),
        'elements': {
            'E1': {
                'name':        'Aromatic hydrophobic core',
                'description': 'pi-stacking cluster around PHE440 anchor',
                'n_residues':  len([f for f in features
                                    if f['element']=='E1']),
                'centroid':    centroids.get('E1', {}),
            },
            'E2': {
                'name':        'Electrostatic anchor',
                'description': 'GLU431 salt bridge with LYS1(NSP7)',
                'n_residues':  len([f for f in features
                                    if f['element']=='E2']),
                'centroid':    centroids.get('E2', {}),
            },
            'E3': {
                'name':        'Distal aromatic patch',
                'description': 'PHE843 hydrophobic contact',
                'n_residues':  len([f for f in features
                                    if f['element']=='E3']),
                'centroid':    centroids.get('E3', {}),
            },
        },
        'features': features,
        'screening_notes': (
            'Compounds targeting this pharmacophore should: '
            '(1) present aromatic/hydrophobic groups to engage '
            'E1 core (PHE440/PHE442/TYR420/PHE415); '
            '(2) optionally carry a negative charge or H-bond '
            'acceptor to mimic GLU431 interaction with NSP7; '
            '(3) fit within 8-12A radius of E1 centroid.'
        ),
    }

    return pharmacophore


# ── Save outputs ─────────────────────────────────────────────
def save_json(pharmacophore):
    path = os.path.join(OUT_DIR,
        'pharmacophore_NSP12-NSP7.json')
    with open(path, 'w') as f:
        json.dump(pharmacophore, f, indent=2)
    print(f"  JSON: {path}")


def save_tsv(pharmacophore):
    path = os.path.join(OUT_DIR,
        'pharmacophore_NSP12-NSP7.tsv')
    fields = ['id','residue','element','role',
              'pharm_type','tolerance_A','interaction',
              'conservation','composite','x','y','z']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields,
                           delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(pharmacophore['features'])
    print(f"  TSV: {path}")


# ── Prism figure — 2x2 ──────────────────────────────────────
def prism_axes(ax, ymax=None, yticks=None):
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    if ymax is not None:
        ax.spines['left'].set_bounds(0, ymax)
        ax.set_ylim(-ymax*0.04, ymax*1.1)
    if yticks is not None:
        ax.set_yticks(yticks)


def make_figure(pharmacophore, centroids):
    features = pharmacophore['features']

    fig = plt.figure(figsize=(13.0, 11.0))
    gs  = fig.add_gridspec(2, 2,
                           hspace=0.62, wspace=0.45,
                           left=0.09, right=0.97,
                           top=0.93, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, 0], projection='3d')
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # Element colors
    EL_COLOR = {
        'E1': '#D7191C',
        'E2': '#4DAC26',
        'E3': '#2166AC',
    }
    EL_FILL = {k: v+'55' for k,v in EL_COLOR.items()}
    EL_LABEL = {
        'E1': 'E1: Aromatic core',
        'E2': 'E2: Electrostatic',
        'E3': 'E3: Distal aromatic',
    }

    # ── Panel A: 3D pharmacophore map ─────────────────────────
    ax_a.set_facecolor('white')
    for pane in [ax_a.xaxis.pane, ax_a.yaxis.pane,
                 ax_a.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor('#DDDDDD')
    ax_a.grid(False)

    plotted_el = set()
    for f in features:
        el  = f['element']
        col = EL_COLOR[el]
        is_anchor = f['role'] == 'ANCHOR'
        size   = 200 if is_anchor else 80
        marker = '*' if is_anchor else 'o'
        label  = EL_LABEL[el] \
                 if el not in plotted_el else None

        ax_a.scatter(f['x'], f['y'], f['z'],
                     color=col, s=size,
                     marker=marker,
                     edgecolors='white',
                     linewidths=0.6,
                     zorder=5, label=label,
                     alpha=0.9)

        ax_a.text(f['x']+0.8, f['y']+0.8,
                  f['z']+1.5,
                  f['id'], fontsize=7,
                  color=col,
                  fontweight='bold'
                  if is_anchor else 'normal')
        plotted_el.add(el)

    # Draw element centroids as transparent spheres
    for el, ct in centroids.items():
        col = EL_COLOR[el]
        ax_a.scatter(ct['x'], ct['y'], ct['z'],
                     color=col, s=400,
                     alpha=0.15, zorder=2,
                     edgecolors=col,
                     linewidths=0.5)

    ax_a.set_xlabel('X (A)', fontsize=7, labelpad=2)
    ax_a.set_ylabel('Y (A)', fontsize=7, labelpad=2)
    ax_a.set_zlabel('Z (A)', fontsize=7, labelpad=2)
    ax_a.tick_params(labelsize=6)
    ax_a.view_init(elev=25, azim=-55)
    ax_a.legend(fontsize=7, loc='upper left',
                frameon=True, facecolor='white',
                edgecolor='#CCCCCC', framealpha=0.9)
    ax_a.set_title(
        'A   3D pharmacophore map\n'
        '(star=PHE440 anchor; sphere=element centroid)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    # ── Panel B: Feature composite bars ───────────────────────
    labels_b = [f['id'] for f in features]
    vals_b   = [f['composite'] for f in features]
    cols_b   = [EL_COLOR[f['element']] for f in features]
    fills_b  = [EL_FILL[f['element']] for f in features]
    x        = np.arange(len(features))
    w        = 0.42

    for xi, val, fc, ec, f in zip(
            x, vals_b, fills_b, cols_b, features):
        ax_b.bar(xi, val, width=w,
                 facecolor=fc, edgecolor=ec,
                 linewidth=0.9, zorder=3,
                 clip_on=False)
        ax_b.scatter(xi, val, color=ec,
                     s=18, zorder=5, clip_on=False)
        ax_b.text(xi, val+0.012, f'{val:.3f}',
                  ha='center', va='bottom',
                  fontsize=6.5, color=ec,
                  fontweight='bold', clip_on=False)
        # Tolerance band
        tol = f['tolerance_A'] / 10.0
        ax_b.errorbar(xi, val, yerr=tol,
                      color=ec, linewidth=0.8,
                      capsize=3, capthick=0.8,
                      zorder=6)

    ax_b.set_xticks(x)
    ax_b.set_xticklabels(labels_b, fontsize=7.5,
                         rotation=90, ha='center',
                         va='top')
    ax_b.tick_params(axis='x', pad=2,
                     length=3, width=0.75)
    ax_b.set_xlim(-0.6, len(features)-0.4)
    ax_b.set_ylabel('CGCP composite score', fontsize=9)
    prism_axes(ax_b, ymax=1.25,
               yticks=[0,0.2,0.4,0.6,0.8,1.0])
    ax_b.set_title(
        'B   Pharmacophore feature scores\n'
        '(error bars = spatial tolerance in Angstrom/10)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    # Element color legend
    handles_b = [
        mpatches.Patch(facecolor=EL_FILL[el],
                       edgecolor=EL_COLOR[el],
                       linewidth=0.9,
                       label=EL_LABEL[el])
        for el in ['E1','E2','E3']
    ]
    ax_b.legend(handles=handles_b, loc='upper right',
                fontsize=7.5, frameon=True,
                facecolor='white', edgecolor='#CCCCCC')

    # ── Panel C: Element summary table ───────────────────────
    ax_c.axis('off')
    ax_c.set_title(
        'C   Pharmacophore element summary\n'
        '(CGCP-NSP12-NSP7-v1)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    table_data = []
    for el_id, el_def in pharmacophore['elements'].items():
        ct = el_def['centroid']
        members = [f['id'] for f in features
                   if f['element'] == el_id]
        table_data.append([
            el_id,
            el_def['name'],
            ', '.join(members),
            f"({ct.get('x',0):.1f}, "
            f"{ct.get('y',0):.1f}, "
            f"{ct.get('z',0):.1f})",
            f"{ct.get('radius',0):.1f} A",
        ])

    col_labels = ['Element','Name','Residues',
                  'Centroid','Radius']
    col_widths = [0.06, 0.22, 0.38, 0.22, 0.10]

    y0 = 0.85
    dy = 0.12

    # Header
    xpos = 0.02
    for lbl, w2 in zip(col_labels, col_widths):
        ax_c.text(xpos, y0, lbl,
                  transform=ax_c.transAxes,
                  fontsize=8, fontweight='bold',
                  color=P_['gray'], va='top')
        xpos += w2

    ax_c.plot([0.01, 0.99], [y0-0.04, y0-0.04],
              color=P_['gray'], linewidth=0.5,
              transform=ax_c.transAxes)

    for i, row in enumerate(table_data):
        yi   = y0 - 0.08 - i*dy
        xpos = 0.02
        el   = row[0]
        for val, w2 in zip(row, col_widths):
            ax_c.text(xpos, yi, val,
                      transform=ax_c.transAxes,
                      fontsize=7.5, color=EL_COLOR.get(
                          el, P_['text']),
                      va='top')
            xpos += w2

    # ── Panel D: Tolerance radius per feature ────────────────
    tols  = [f['tolerance_A'] for f in features]
    types = [f['pharm_type'] for f in features]
    tcols = [FEAT_COLORS.get(t, '#636363')
             for t in types]
    tfils = [FEAT_FILLS.get(t, '#63636355')
             for t in types]

    for xi, val, fc, ec in zip(
            x, tols, tfils, tcols):
        ax_d.bar(xi, val, width=w,
                 facecolor=fc, edgecolor=ec,
                 linewidth=0.9, zorder=3,
                 clip_on=False)
        ax_d.text(xi, val+0.04, f'{val} A',
                  ha='center', va='bottom',
                  fontsize=7.5, color=ec,
                  fontweight='bold', clip_on=False)

    ax_d.set_xticks(x)
    ax_d.set_xticklabels(labels_b, fontsize=7.5,
                         rotation=90, ha='center',
                         va='top')
    ax_d.tick_params(axis='x', pad=2,
                     length=3, width=0.75)
    ax_d.set_xlim(-0.6, len(features)-0.4)
    ax_d.set_ylabel('Spatial tolerance (A)', fontsize=9)
    prism_axes(ax_d, ymax=3.0,
               yticks=[0, 0.5, 1.0, 1.5, 2.0, 2.5])
    ax_d.set_title(
        'D   Spatial tolerance per feature\n'
        '(tolerance = acceptable displacement for docking)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    # Feature type legend
    handles_d = [
        mpatches.Patch(facecolor=FEAT_FILLS[t],
                       edgecolor=FEAT_COLORS[t],
                       linewidth=0.9, label=t)
        for t in ['AROMATIC','HYDROPHOBIC','NEG_CHARGE']
    ]
    ax_d.legend(handles=handles_d, loc='upper right',
                fontsize=7.5, frameon=True,
                facecolor='white', edgecolor='#CCCCCC')

    path = os.path.join(OUT_DIR,
        'Fig_Step07_PharmacophorHypothesis.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


# ── Print summary ────────────────────────────────────────────
def print_summary(pharmacophore, centroids):
    print('\n' + '='*65)
    print('STEP 7 PHARMACOPHORE HYPOTHESIS — NSP12-NSP7')
    print(f"Name: {pharmacophore['name']}")
    print('='*65)

    for el_id, el_def in pharmacophore['elements'].items():
        ct = el_def['centroid']
        print(f"\n  {el_id}: {el_def['name']}")
        print(f"    Description: {el_def['description']}")
        print(f"    Residues: {el_def['n_residues']}")
        print(f"    Centroid: ({ct.get('x',0):.2f}, "
              f"{ct.get('y',0):.2f}, "
              f"{ct.get('z',0):.2f})")
        print(f"    Radius: {ct.get('radius',0):.2f} A")
        print(f"    Members:")
        for f in pharmacophore['features']:
            if f['element'] == el_id:
                print(f"      {f['id']:<10} "
                      f"{f['pharm_type']:<15} "
                      f"tol={f['tolerance_A']} A  "
                      f"cons={f['conservation']:.3f}  "
                      f"comp={f['composite']:.3f}")

    print(f"\n  Screening notes:")
    print(f"  {pharmacophore['screening_notes']}")
    print('='*65)


# ── Patch for panel C color ref ──────────────────────────────
P_ = {'gray': '#636363', 'text': '#1A1A1A'}


# ── Main ─────────────────────────────────────────────────────
if __name__ == '__main__':
    print('CGCP Phase 2 Step 7 - Pharmacophore Hypothesis')

    coords    = get_coords()
    centroids = compute_centroids(coords)
    pharma    = build_pharmacophore(coords, centroids)

    save_json(pharma)
    save_tsv(pharma)
    make_figure(pharma, centroids)
    print_summary(pharma, centroids)

    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Phase 2 Step 8 - controls + ChimeraX figures")
