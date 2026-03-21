import json, os, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio.PDB import PDBParser

matplotlib.rcParams.update({
    'font.family':'DejaVu Sans','font.size':11,
    'axes.titlesize':12,'axes.titleweight':'bold',
    'axes.labelsize':11,'axes.labelweight':'bold',
    'axes.linewidth':1.2,'xtick.labelsize':10,
    'ytick.labelsize':10,'figure.dpi':150,
})

RES_DIR = "results"
HOTSPOT_PDB = [
    10,11,12,13,14,15,26,28,29,33,36,38,39,40,41,42,43,
    47,48,49,59,62,64,91,95,96,97,163,164,166,168,169,
    170,171,172,173,203,241,242,243,244,263,265,267,269,
    270,271,272,280,282,283,284,285,286,287,288,289,291,292
]
STRUCTS = {
    "6VWW":{"res":1.90,"label":"6VWW\n(apo)","ligand":None,"primary":True},
    "6W01":{"res":2.20,"label":"6W01\n(apo)","ligand":None,"primary":False},
    "6WLC":{"res":1.70,"label":"6WLC\n(+UMP)","ligand":"UMP\n(active site)","primary":False},
    "6WXC":{"res":1.80,"label":"6WXC\n(+tipiracil)","ligand":"tipiracil","primary":False},
    "2H85":{"res":2.00,"label":"2H85\n(SARS-1)","ligand":None,"primary":False},
}

# B-factors
bf = {}
for sid in STRUCTS:
    path = f"00-reference/known_interfaces/NSP15/{sid}_NSP15-true-dimer.pdb"
    if not os.path.exists(path):
        path = f"00-reference/pdb_structures/NSP15/{sid}.pdb"
    if not os.path.exists(path):
        bf[sid] = 0.0; continue
    s = PDBParser(QUIET=True).get_structure(sid, path)
    vals = [r['CA'].get_bfactor()
            for c in s[0].get_chains()
            for r in c.get_residues()
            if r.get_id()[0]==' ' and r.get_id()[1] in HOTSPOT_PDB and 'CA' in r]
    bf[sid] = float(np.mean(vals)) if vals else 0.0

sids   = list(STRUCTS.keys())
labels = [STRUCTS[s]["label"] for s in sids]
res    = [STRUCTS[s]["res"]   for s in sids]
bfv    = [bf[s]               for s in sids]
cols   = ['#B22222' if STRUCTS[s]["primary"] else '#5B9BD5' for s in sids]
x      = np.arange(len(sids))
w      = 0.55

# ── Figure: 3 rows — title | charts | caption ──────────────────
fig = plt.figure(figsize=(13, 7.5))
fig.patch.set_facecolor('white')

# Row heights: charts get most space, caption gets fixed bottom strip
gs = fig.add_gridspec(2, 2,
    height_ratios=[10, 1.2],
    hspace=0.05,
    left=0.07, right=0.97,
    top=0.88, bottom=0.08,
    wspace=0.36)

ax1  = fig.add_subplot(gs[0, 0])
ax2  = fig.add_subplot(gs[0, 1])
axc  = fig.add_subplot(gs[1, :])   # caption row

# ── Main title ─────────────────────────────────────────────────
fig.suptitle("NSP15 — Tier 0: Structure Selection Evidence",
             fontsize=14, fontweight='bold', y=0.97)

# ── Panel 1: Resolution ────────────────────────────────────────
bars1 = ax1.bar(x, res, width=w, color=cols,
                edgecolor='white', linewidth=1.0, zorder=3)
ax1.set_ylim(0, 2.8)
ax1.invert_yaxis()
ax1.set_xticks(x)
ax1.set_xticklabels(labels, fontsize=10, fontweight='bold')
ax1.set_ylabel("Resolution (Å)", fontweight='bold')
ax1.set_title("Crystal Structure Resolution\n(lower = better)", pad=6)
ax1.yaxis.grid(True, alpha=0.3, linestyle='--', zorder=0)
ax1.set_axisbelow(True)
for sp in ['top','right']:
    ax1.spines[sp].set_visible(False)

for bar, val, sid in zip(bars1, res, sids):
    ax1.text(bar.get_x() + bar.get_width()/2,
             bar.get_height() + 0.05,
             f"{val:.1f} Å", ha='center', va='top',
             fontsize=10, fontweight='bold',
             color='#8B0000' if STRUCTS[sid]["primary"] else '#1A5276')
    if STRUCTS[sid]["ligand"]:
        ax1.text(bar.get_x() + bar.get_width()/2,
                 bar.get_height()/2,
                 STRUCTS[sid]["ligand"],
                 ha='center', va='center', fontsize=8,
                 color='white', fontweight='bold', rotation=90)

ax1.annotate('Best resolution\nbut ligand-bound\n→ not selected',
             xy=(2, 1.70), xytext=(3.5, 1.20),
             arrowprops=dict(arrowstyle='->', color='#555', lw=1.2),
             fontsize=8.5, color='#555555', ha='center',
             bbox=dict(boxstyle='round,pad=0.3', fc='#FFFDE7',
                       ec='#F9A825', lw=0.9))

# ── Panel 2: B-factor ──────────────────────────────────────────
bars2 = ax2.bar(x, bfv, width=w, color=cols,
                edgecolor='white', linewidth=1.0, zorder=3)
ax2.set_ylim(0, 50)
ax2.set_xticks(x)
ax2.set_xticklabels(labels, fontsize=10, fontweight='bold')
ax2.set_ylabel("Mean B-factor at hotspots (Å²)", fontweight='bold')
ax2.set_title("Interface Hotspot Flexibility\n(lower = more ordered)", pad=6)
ax2.yaxis.grid(True, alpha=0.3, linestyle='--', zorder=0)
ax2.set_axisbelow(True)
for sp in ['top','right']:
    ax2.spines[sp].set_visible(False)

for bar, val, sid in zip(bars2, bfv, sids):
    ax2.text(bar.get_x() + bar.get_width()/2,
             bar.get_height() + 0.6,
             f"{val:.1f}", ha='center', va='bottom',
             fontsize=10, fontweight='bold',
             color='#8B0000' if STRUCTS[sid]["primary"] else '#1A5276')

ax2.annotate('6VWW selected:\napo form + verified\nbiological assembly\n(75 contacts)',
             xy=(0, bfv[0]), xytext=(1.6, 45),
             arrowprops=dict(arrowstyle='->', color='#B22222', lw=1.4),
             fontsize=8.5, color='#B22222', ha='center',
             bbox=dict(boxstyle='round,pad=0.3', fc='#FFF0F0',
                       ec='#B22222', lw=0.9))

# ── Shared legend between charts and caption ───────────────────
legend_patches = [
    mpatches.Patch(color='#B22222',
        label='Selected receptor — 6VWW (apo, 1.90 Å, biological assembly A+A-2 verified)'),
    mpatches.Patch(color='#5B9BD5',
        label='Other structures — used for cross-validation or active-site ground truth'),
]
fig.legend(handles=legend_patches, loc='lower center',
           ncol=1, frameon=True, fontsize=9,
           bbox_to_anchor=(0.5, 0.095),
           fancybox=False, edgecolor='#CCCCCC',
           handlelength=1.5, handletextpad=0.6)

# ── Caption axes — no frame, just text ────────────────────────
axc.axis('off')
caption = (
    "Fig. A1 | Structure selection evidence for NSP15 NendoU homodimer. "
    "6VWW (apo, 1.90 Å) was selected as the primary receptor: it has no active-site ligand distortion, "
    "and its biological assembly (chains A+A-2) shows 75 inter-chain contacts versus only 7 in the asymmetric unit — "
    "confirming the true dimer interface. 6WLC has the best resolution (1.70 Å) but UMP binding "
    "distorts the active site; it is retained as ground truth for druggability validation (Script 16_9)."
)
axc.text(0.5, 0.8, caption,
         ha='center', va='top', fontsize=8.5,
         color='#333333', style='italic',
         transform=axc.transAxes,
         wrap=True)

plt.savefig(f"{RES_DIR}/FigA1_NSP15_structure_justification_9.png",
            dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved FigA1")
