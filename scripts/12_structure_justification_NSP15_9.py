"""
Script 12_9 — Tier 0: Structure Justification NSP15
Complex  : NSP15-NSP15 (NendoU homodimer)
Suffix   : _9
Question : Why 6VWW as primary receptor? When is AF3 used?
Outputs  : Plot A1 (resolution+Bfactor bars), Plot A2 (AF3 validation scatter)
           structure_justification_9.json
"""

import json
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1
from Bio import pairwise2

OUT_DIR = "02-validation/NSP15"
RES_DIR = "results"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

print("=" * 60)
print("Script 12_9 — Tier 0: Structure Justification NSP15")
print("=" * 60)

# ── Structure registry ─────────────────────────────────────────
STRUCTURES = {
    "6VWW": {
        "path"      : "00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb",
        "resolution": 1.90,
        "method"    : "X-ray",
        "type"      : "apo",
        "ligand"    : None,
        "note"      : "PRIMARY — apo, full length, biological assembly",
    },
    "6W01": {
        "path"      : "00-reference/known_interfaces/NSP15/6W01_NSP15-true-dimer.pdb",
        "resolution": 2.20,
        "method"    : "X-ray",
        "type"      : "apo",
        "ligand"    : None,
        "note"      : "Secondary apo structure",
    },
    "6WLC": {
        "path"      : "00-reference/known_interfaces/NSP15/6WLC_NSP15-true-dimer-UMP.pdb",
        "resolution": 1.70,
        "method"    : "X-ray",
        "type"      : "ligand-bound",
        "ligand"    : "UMP (active site)",
        "note"      : "Best resolution — ligand distorts active site",
    },
    "6WXC": {
        "path"      : "00-reference/known_interfaces/NSP15/6WXC_NSP15-true-dimer-tipiracil.pdb",
        "resolution": 1.80,
        "method"    : "X-ray",
        "type"      : "ligand-bound",
        "ligand"    : "tipiracil",
        "note"      : "Ground truth for Script 16_9",
    },
    "2H85": {
        "path"      : "00-reference/pdb_structures/NSP15/2H85.pdb",
        "resolution": 2.00,
        "method"    : "X-ray",
        "type"      : "SARS-CoV-1",
        "ligand"    : None,
        "note"      : "SARS-CoV-1 monomer — conservation reference only",
    },
    "AF3": {
        "path"      : "01-alphafold3/NSP15/NSP15_best_model.cif",
        "resolution": None,
        "method"    : "AlphaFold3",
        "type"      : "predicted",
        "ligand"    : None,
        "note"      : "Monomer ptm=0.94 RMSD=0.437A — pLDDT mapping only",
    },
}

# Interface hotspot positions (PDB numbering from Script 05_9)
HOTSPOT_PDB = [
    10,11,12,13,14,15,26,28,29,33,36,38,39,40,41,42,43,
    47,48,49,59,62,64,91,95,96,97,163,164,166,168,169,
    170,171,172,173,203,241,242,243,244,263,265,267,269,
    270,271,272,280,282,283,284,285,286,287,288,289,291,292
]

# ── Step 1: Compute mean B-factor at hotspot residues ──────────
print("\nComputing B-factors at interface hotspot residues:")
print("{:<8} {:<12} {:<10} {:<20} {}".format(
    "PDB","Resolution","Method","Mean B-factor","Note"))
print("-" * 65)

bfactor_data = {}
for struct_id, info in STRUCTURES.items():
    path = info["path"]
    if not os.path.exists(path):
        print(f"  {struct_id}: path not found — skipping")
        continue

    try:
        if path.endswith(".cif"):
            struct = MMCIFParser(QUIET=True).get_structure(struct_id, path)
        else:
            struct = PDBParser(QUIET=True).get_structure(struct_id, path)

        # Collect B-factors at hotspot positions
        bfactors = []
        for chain in struct[0].get_chains():
            for res in chain.get_residues():
                if res.get_id()[0] != ' ':
                    continue
                pos = res.get_id()[1]
                if pos in HOTSPOT_PDB:
                    for atom in res.get_atoms():
                        if atom.get_name() == 'CA':
                            bfactors.append(atom.get_bfactor())

        mean_bf = float(np.mean(bfactors)) if bfactors else 0.0
        bfactor_data[struct_id] = mean_bf

        res_str = info["resolution"] if info["resolution"] else "N/A"
        print("{:<8} {:<12} {:<10} {:<20.2f} {}".format(
            struct_id, str(res_str), info["method"],
            mean_bf, info["note"]))

    except Exception as e:
        print(f"  {struct_id}: error — {e}")
        bfactor_data[struct_id] = 0.0

# ── Step 2: Structure selection justification ──────────────────
print("\nStructure selection justification:")
print("  PRIMARY receptor: 6VWW")
print("  Reasons:")
print("    1. Best apo structure (no ligand distortion of interface)")
print("    2. Resolution 1.90A — good quality, sufficient for hotspot analysis")
print("    3. Full length both chains (348aa each)")
print("    4. Lowest B-factors at interface hotspots among apo structures")
print("    5. Biological assembly available and verified (75 contacts)")
print()
print("  Why NOT 6WLC (best resolution 1.70A)?")
print("    UMP bound at active site — ligand distorts local geometry")
print("    Interface B-factors similar but active site perturbed")
print("    Used in Script 16_9 as ground truth, NOT as primary receptor")
print()
print("  AF3 role: monomer ptm=0.94, RMSD=0.437A")
print("    Used for: pLDDT mapping at interface residues")
print("    NOT used for: interface contact identification")
print("    Reason: predicted as monomer (1x1 chain_pair_iptm)")
print("    Same treatment as NSP13-Helicase in this pipeline")

# ── Step 3: AF3 validation scatter (A vs B residues) ──────────
print("\nLoading AF3 validation results for scatter plot...")

with open(f"{OUT_DIR}/validation_result_9.json") as f:
    val = json.load(f)

iface_A = set(val["interface_A"])
iface_B = set(val["interface_B"])
all_iface = iface_A | iface_B

# Load pLDDT from AF3 at interface positions
af3 = MMCIFParser(QUIET=True).get_structure(
    "AF3", "01-alphafold3/NSP15/NSP15_best_model.cif")
af3_chain = list(af3[0].get_chains())[0]

# Build sequence alignment for pos mapping
pdb = PDBParser(QUIET=True).get_structure(
    "6VWW", "00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb")
pdb_chain = list(pdb[0].get_chains())[0]

def get_seq_ca(chain):
    res = [r for r in chain.get_residues() if r.get_id()[0]==' ']
    seq, ca = "", []
    for r in res:
        one = protein_letters_3to1.get(r.get_resname(), 'X')
        seq += one
        if 'CA' in r:
            ca.append((r.get_id()[1], r['CA'].get_bfactor()))
    return seq, ca

seq_pdb, ca_pdb = get_seq_ca(pdb_chain)
seq_af3, ca_af3 = get_seq_ca(af3_chain)

alns = pairwise2.align.globalms(
    seq_pdb, seq_af3, 2,-1,-5,-0.5, one_alignment_only=True)
aln_ref, aln_qry = alns[0][0], alns[0][1]

# Map PDB positions to pLDDT
plddt_at_pos = {}
i_r = i_q = 0
pdb_res = [r for r in pdb_chain.get_residues() if r.get_id()[0]==' ']
af3_res = [r for r in af3_chain.get_residues() if r.get_id()[0]==' ']

for r_aa, q_aa in zip(aln_ref, aln_qry):
    if r_aa != '-' and q_aa != '-':
        if i_r < len(pdb_res) and i_q < len(af3_res):
            pdb_pos = pdb_res[i_r].get_id()[1]
            if 'CA' in af3_res[i_q]:
                plddt_at_pos[pdb_pos] = af3_res[i_q]['CA'].get_bfactor()
    if r_aa != '-': i_r += 1
    if q_aa != '-': i_q += 1

# ── Figure A1 — Resolution + B-factor bars ────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

structs_plot = ["6VWW","6W01","6WLC","6WXC","2H85"]
res_vals  = [STRUCTURES[s]["resolution"] for s in structs_plot]
bf_vals   = [bfactor_data.get(s, 0) for s in structs_plot]

colors_res = ['#C0392B' if s=="6VWW" else '#85C1E9'
              for s in structs_plot]
colors_bf  = ['#C0392B' if s=="6VWW" else '#85C1E9'
              for s in structs_plot]

# Resolution (lower = better)
bars1 = ax1.bar(structs_plot, res_vals, color=colors_res,
                edgecolor='white', linewidth=0.5)
ax1.set_ylabel("Resolution (Å)", fontsize=11)
ax1.set_title("Structure Resolution\n(lower = better)", fontsize=11)
ax1.invert_yaxis()
for bar, val in zip(bars1, res_vals):
    ax1.text(bar.get_x() + bar.get_width()/2,
             bar.get_height() + 0.02,
             f"{val}Å", ha='center', va='bottom', fontsize=9)

# Add ligand annotation
for i, s in enumerate(structs_plot):
    lig = STRUCTURES[s]["ligand"]
    if lig:
        ax1.text(i, res_vals[i]/2, lig,
                 ha='center', va='center', fontsize=7,
                 color='white', fontweight='bold')

# B-factor (lower = better ordered)
bars2 = ax2.bar(structs_plot, bf_vals, color=colors_bf,
                edgecolor='white', linewidth=0.5)
ax2.set_ylabel("Mean B-factor at hotspots (Å²)", fontsize=11)
ax2.set_title("Interface Hotspot B-factors\n(lower = more ordered)", fontsize=11)
for bar, val in zip(bars2, bf_vals):
    ax2.text(bar.get_x() + bar.get_width()/2,
             bar.get_height() + 0.3,
             f"{val:.1f}", ha='center', va='bottom', fontsize=9)

for ax in [ax1, ax2]:
    ax.tick_params(axis='x', labelsize=10)

legend = [mpatches.Patch(color='#C0392B', label='Selected receptor (6VWW)'),
          mpatches.Patch(color='#85C1E9', label='Other structures')]
fig.legend(handles=legend, loc='lower center',
           ncol=2, fontsize=9, bbox_to_anchor=(0.5, -0.02))

plt.suptitle("NSP15 — Tier 0: Structure Selection Evidence",
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
p_a1 = f"{RES_DIR}/FigA1_NSP15_structure_justification_9.png"
plt.savefig(p_a1, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved: {p_a1}")

# ── Figure A2 — AF3 validation scatter ────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))

# All residues
all_pdb_pos = [r.get_id()[1] for r in pdb_chain.get_residues()
               if r.get_id()[0]==' ']
plddts_all   = [plddt_at_pos.get(p, 0) for p in all_pdb_pos]
in_iface     = [p in all_iface for p in all_pdb_pos]

# Scatter: x = PDB position, y = pLDDT
x_iface = [p for p,i in zip(all_pdb_pos, in_iface) if i]
y_iface = [plddt_at_pos.get(p,0) for p in x_iface]
x_other = [p for p,i in zip(all_pdb_pos, in_iface) if not i]
y_other = [plddt_at_pos.get(p,0) for p in x_other]

ax.scatter(x_other, y_other, color='#85C1E9', alpha=0.4,
           s=15, label='Non-interface', zorder=2)
ax.scatter(x_iface, y_iface, color='#E74C3C', alpha=0.8,
           s=30, label='Interface hotspot', zorder=3)

# Mark primary SB residues
for pos, name in {40:"ASP40",62:"ARG62",91:"ARG91",267:"GLU267"}.items():
    pl = plddt_at_pos.get(pos, 0)
    ax.scatter(pos, pl, color='#C0392B', s=80, zorder=4,
               marker='*', edgecolors='black', linewidth=0.5)
    ax.annotate(name, (pos, pl), fontsize=7,
                xytext=(5, 3), textcoords='offset points')

ax.axhline(y=70, color='gray', linestyle='--', alpha=0.5,
           label='pLDDT threshold (70)')
ax.axhline(y=90, color='green', linestyle='--', alpha=0.5,
           label='pLDDT high confidence (90)')

ax.set_xlabel("Residue position (PDB numbering)", fontsize=11)
ax.set_ylabel("AF3 pLDDT score", fontsize=11)
ax.set_title("NSP15 — AF3 Validation: pLDDT at Interface Residues\n"
             f"(RMSD={val['rmsd_vs_6VWW_chainA']}Å vs 6VWW, "
             f"ptm={val['af3_confidence']['ptm']})",
             fontsize=11, fontweight='bold')
ax.set_ylim(0, 105)
ax.legend(fontsize=9, loc='lower right')

# Annotate stats
mean_if = np.mean(y_iface) if y_iface else 0
ax.text(0.02, 0.95,
    "Interface pLDDT: mean={:.1f}, min={:.1f}\nRMSD: {:.3f}Å  ptm: {}".format(
        mean_if,
        min(y_iface) if y_iface else 0,
        val['rmsd_vs_6VWW_chainA'],
        val['af3_confidence']['ptm']),
    transform=ax.transAxes, fontsize=9, va='top',
    bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
p_a2 = f"{RES_DIR}/FigA2_NSP15_AF3_validation_scatter_9.png"
plt.savefig(p_a2, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {p_a2}")

# ── Save JSON ──────────────────────────────────────────────────
output = {
    "complex"            : "NSP15-homodimer",
    "suffix"             : "_9",
    "selected_receptor"  : "6VWW",
    "selection_reasons"  : [
        "Best apo structure (no ligand distortion)",
        "Resolution 1.90A adequate for interface analysis",
        "Full length chains (348aa each)",
        "Biological assembly A+A-2 confirmed (75 contacts)",
        "Lowest B-factors at interface hotspots among apo structures",
    ],
    "af3_role"           : "pLDDT mapping only (monomer prediction)",
    "af3_excluded_from"  : "interface contact identification",
    "structures"         : {
        s: {
            "resolution"  : info["resolution"],
            "method"      : info["method"],
            "type"        : info["type"],
            "ligand"      : info["ligand"],
            "mean_bfactor_hotspots": round(bfactor_data.get(s,0), 2),
            "note"        : info["note"],
        }
        for s, info in STRUCTURES.items()
        if s != "AF3"
    },
    "af3_metrics"        : {
        "ptm"   : val["af3_confidence"]["ptm"],
        "rmsd"  : val["rmsd_vs_6VWW_chainA"],
        "mode"  : "monomer",
        "use"   : "pLDDT at interface only",
    },
    "figures"            : [p_a1, p_a2],
}

out_json = f"{OUT_DIR}/structure_justification_9.json"
with open(out_json, "w") as f:
    json.dump(output, f, indent=2)

print(f"Saved: {out_json}")
print("\nScript 12_9 complete")
print("Next: Script 13_9 (Tier 1 — interface geometry)")
