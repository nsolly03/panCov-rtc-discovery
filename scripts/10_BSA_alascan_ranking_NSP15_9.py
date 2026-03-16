"""
Script 10_9 — BSA + Alanine Scanning + Composite Ranking NSP15
Complex  : NSP15-NSP15 (NendoU homodimer)
Suffix   : _9
Receptor : 6VWW_NSP15-true-dimer.pdb
Controls :
  Positive: interface residues validated by hexamerization-abolishing
            mutations (Deng and Baker 2018, Virology 517:157-163)
  Negative: His235, His250, Lys290 (catalytic triad, active site,
            ZERO dimer interface contacts)
"""

import json
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

OUT_DIR  = "02-validation/NSP15"
RES_DIR  = "results"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

PDB_PATH = "00-reference/known_interfaces/NSP15/6VWW_NSP15-true-dimer.pdb"

print("=" * 60)
print("Script 10_9 — BSA + AlaScan + Ranking NSP15")
print("=" * 60)

# ── Load data ──────────────────────────────────────────────────
with open(f"{OUT_DIR}/interface_analysis_9.json") as f:
    iface = json.load(f)

with open(f"{OUT_DIR}/conservation_summary_9.json") as f:
    cons_data = json.load(f)

# Build conservation map: UniProt pos -> score
# PDB pos = UniProt pos + 1 (offset confirmed in Script 06_9)
OFFSET = 1
cons_map = {}
for h in cons_data["hotspot_details"]:
    uni_pos = h["position"]
    pdb_pos = uni_pos + OFFSET
    cons_map[pdb_pos] = h["score"]

# Consensus hotspots from Script 05_9 (PDB numbering)
consensus_A = set(iface["consensus_A"])
consensus_B = set(iface["consensus_B"])
all_hotspots_pdb = sorted(consensus_A | consensus_B)

# Controls (PDB numbering)
POS_CONTROLS = {40: "ASP40", 62: "ARG62", 91: "ARG91", 267: "GLU267"}
NEG_CONTROLS = {235: "HIS235", 250: "HIS250", 290: "LYS290"}

print(f"\nHotspots loaded: {len(all_hotspots_pdb)} (PDB numbering)")
print(f"Positive controls: {list(POS_CONTROLS.values())}")
print(f"Negative controls: {list(NEG_CONTROLS.values())}")

# ── Load structure ─────────────────────────────────────────────
pdb = PDBParser(QUIET=True).get_structure("6VWW", PDB_PATH)
model = pdb[0]
chains = list(model.get_chains())
chain_A = chains[0]
chain_B = chains[1]

# ── BSA calculation ────────────────────────────────────────────
print("\nCalculating BSA (ShrakeRupley)...")

sr = ShrakeRupley()

# Unbound SASA (each chain alone)
from Bio.PDB import Structure as PDBStruct, Model as PDBModel

def get_chain_sasa(chain):
    tmp = PDBStruct.Structure("tmp")
    m   = PDBModel.Model(0)
    tmp.add(m)
    m.add(chain.copy())
    sr.compute(tmp, level="R")
    sasa = {}
    for res in m.get_residues():
        if res.get_id()[0] == ' ':
            sasa[res.get_id()[1]] = res.sasa
    return sasa

sasa_A_unbound = get_chain_sasa(chain_A)
sasa_B_unbound = get_chain_sasa(chain_B)

# Complex SASA
sr.compute(pdb, level="R")
sasa_A_complex = {}
sasa_B_complex = {}
for chain in model.get_chains():
    for res in chain.get_residues():
        if res.get_id()[0] == ' ':
            pos = res.get_id()[1]
            if chain.id == chain_A.id:
                sasa_A_complex[pos] = res.sasa
            else:
                sasa_B_complex[pos] = res.sasa

# BSA per residue
bsa = {}
for pos in all_hotspots_pdb:
    bsa_a = (sasa_A_unbound.get(pos,0) - sasa_A_complex.get(pos,0))
    bsa_b = (sasa_B_unbound.get(pos,0) - sasa_B_complex.get(pos,0))
    bsa[pos] = max(0, max(bsa_a, bsa_b))

# Also compute for negative controls
for pos in NEG_CONTROLS:
    bsa_a = (sasa_A_unbound.get(pos,0) - sasa_A_complex.get(pos,0))
    bsa[pos] = max(0, bsa_a)

print("BSA top 10:")
for pos, val in sorted(bsa.items(), key=lambda x: -x[1])[:10]:
    if pos in all_hotspots_pdb:
        aa = next((r.get_resname() for r in chain_A.get_residues()
                   if r.get_id()[1]==pos and r.get_id()[0]==' '), '?')
        print(f"  {aa}{pos}: {val:.1f} A2")

# ── Computational alanine scanning ────────────────────────────
print("\nComputational alanine scanning...")

# Contact types per residue from Script 05_9
contact_loss = {}
sb_contacts  = {}
hb_contacts  = {}
hy_contacts  = {}

# Load raw contacts from interface analysis
vww_data = iface["per_structure"].get("6VWW", {})
contacts_raw = iface.get("per_structure",{}).get("6VWW",{})

# Recompute from structure directly for accuracy
CUTOFF = 5.0
POS_CH = {'ARG','LYS','HIS'}
NEG_CH = {'ASP','GLU'}
HYDRO  = {'ALA','VAL','LEU','ILE','MET','PHE','TRP','PRO','TYR'}

res_A = [r for r in chain_A.get_residues() if r.get_id()[0]==' ']
res_B = [r for r in chain_B.get_residues() if r.get_id()[0]==' ']

for rA in res_A:
    posA = rA.get_id()[1]
    if posA not in all_hotspots_pdb:
        continue
    sb = hb = hy = 0
    for rB in res_B:
        in_contact = False
        for aA in rA.get_atoms():
            for aB in rB.get_atoms():
                if aA - aB < CUTOFF:
                    in_contact = True
                    break
            if in_contact: break
        if in_contact:
            nA, nB = rA.get_resname(), rB.get_resname()
            if ((nA in POS_CH and nB in NEG_CH) or
                (nA in NEG_CH and nB in POS_CH)):
                sb += 1
            elif any(aA.element in ('N','O')
                     for aA in rA.get_atoms()
                     if any(aB.element in ('N','O')
                            for aB in rB.get_atoms()
                            if aA - aB < 3.5)):
                hb += 1
            else:
                hy += 1
    total = sb + hb + hy
    contact_loss[posA] = total
    sb_contacts[posA]  = sb
    hb_contacts[posA]  = hb
    hy_contacts[posA]  = hy

# Negative controls should show near-zero interface contacts
print("Negative control contact check:")
for pos, name in NEG_CONTROLS.items():
    cl = contact_loss.get(pos, 0)
    print(f"  {name}: interface contacts = {cl}  "
          f"(expected ~0 — active site, not dimer interface)")

# ── Composite ranking ──────────────────────────────────────────
print("\nComputing composite ranking...")

max_bsa  = max(bsa[p] for p in all_hotspots_pdb if p in bsa) or 1
max_loss = max(contact_loss.get(p,0) for p in all_hotspots_pdb) or 1

scores = {}
for pos in all_hotspots_pdb:
    aa   = next((r.get_resname() for r in chain_A.get_residues()
                 if r.get_id()[1]==pos and r.get_id()[0]==' '), 'UNK')
    cons = cons_map.get(pos, 0.2)
    bur  = bsa.get(pos, 0) / max_bsa
    enrg = contact_loss.get(pos, 0) / max_loss
    # Interface score (present in how many structures)
    n_structs = sum(1 for s in ["6VWW","6W01","6WLC","6WXC"]
                    if pos in iface["per_structure"].get(s,{}).get("iface_A",[]) +
                               iface["per_structure"].get(s,{}).get("iface_B",[]))
    iface_score = n_structs / 4.0

    composite = iface_score * cons * (0.5*bur + 0.5*enrg)
    scores[pos] = {
        "aa"           : aa,
        "position_pdb" : pos,
        "position_uni" : pos - OFFSET,
        "conservation" : round(cons, 3),
        "bsa"          : round(bsa.get(pos,0), 1),
        "contact_loss" : contact_loss.get(pos, 0),
        "sb_loss"      : sb_contacts.get(pos, 0),
        "hb_loss"      : hb_contacts.get(pos, 0),
        "hy_loss"      : hy_contacts.get(pos, 0),
        "iface_score"  : round(iface_score, 3),
        "composite"    : round(composite, 4),
        "is_pos_ctrl"  : pos in POS_CONTROLS,
        "is_neg_ctrl"  : pos in NEG_CONTROLS,
        "residue_aa"   : f"{aa}{pos}",
    }

ranked = sorted(scores.values(), key=lambda x: -x["composite"])

print(f"\nComposite ranking top 10:")
print("{:<6} {:<10} {:<8} {:<8} {:<10} {:<12} {}".format(
    "Rank","Residue","Cons","BSA","Loss","Composite","Control"))
print("-" * 65)
for i, r in enumerate(ranked[:10], 1):
    ctrl = "★ POS" if r["is_pos_ctrl"] else (
           "✗ NEG" if r["is_neg_ctrl"] else "")
    print("{:<6} {:<10} {:<8.3f} {:<8.1f} {:<10} {:<12.4f} {}".format(
        i, r["residue_aa"], r["conservation"],
        r["bsa"], r["contact_loss"], r["composite"], ctrl))

# Pipeline validation check
print("\nControl validation:")
pos_ranks = {r["position_pdb"]: i+1
             for i, r in enumerate(ranked)
             if r["position_pdb"] in POS_CONTROLS}
for pos, name in POS_CONTROLS.items():
    rank = pos_ranks.get(pos, "not ranked")
    status = "✓ PASS" if isinstance(rank,int) and rank <= 10 else "✗ CHECK"
    print(f"  {name}: rank={rank}  {status}")

# ── Save CSV ───────────────────────────────────────────────────
csv_path = f"{OUT_DIR}/composite_ranking_NSP15_9.csv"
with open(csv_path, "w") as f:
    f.write("rank,residue_aa,position_pdb,position_uni,aa,"
            "conservation,bsa,contact_loss,sb_loss,hb_loss,hy_loss,"
            "iface_score,composite,is_pos_ctrl,is_neg_ctrl\n")
    for i, r in enumerate(ranked, 1):
        f.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
            i, r["residue_aa"], r["position_pdb"], r["position_uni"],
            r["aa"], r["conservation"], r["bsa"], r["contact_loss"],
            r["sb_loss"], r["hb_loss"], r["hy_loss"],
            r["iface_score"], r["composite"],
            r["is_pos_ctrl"], r["is_neg_ctrl"]))

print(f"\nSaved: {csv_path}")

# ── Figure 4 — BSA bar chart ───────────────────────────────────
top20 = ranked[:20]
fig, ax = plt.subplots(figsize=(10, 6))
colors = ['#C0392B' if r["is_pos_ctrl"] else
          '#7F8C8D' if r["is_neg_ctrl"] else
          '#2980B9' for r in top20]
ax.barh([r["residue_aa"] for r in top20],
        [r["bsa"] for r in top20],
        color=colors, edgecolor='white')
ax.set_xlabel("BSA (Å²)", fontsize=11)
ax.set_title("NSP15 Homodimer — BSA per Hotspot Residue (Top 20)",
             fontsize=12, fontweight='bold')
ax.invert_yaxis()
legend_elements = [
    mpatches.Patch(color='#C0392B', label='Primary SB anchor ★'),
    mpatches.Patch(color='#2980B9', label='Interface hotspot'),
]
ax.legend(handles=legend_elements, fontsize=9)
plt.tight_layout()
p4 = f"{RES_DIR}/Fig4_NSP15_BSA_9.png"
plt.savefig(p4, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {p4}")

# ── Figure 5 — Alanine scan bar chart ─────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
top20_loss = sorted(scores.values(),
                    key=lambda x: -x["contact_loss"])[:20]
colors5 = ['#C0392B' if r["is_pos_ctrl"] else '#2ECC71'
           for r in top20_loss]
bottom = np.zeros(len(top20_loss))
for ctype, col, key in [("SB","#E74C3C","sb_loss"),
                         ("HB","#3498DB","hb_loss"),
                         ("HY","#2ECC71","hy_loss")]:
    vals = [r[key] for r in top20_loss]
    ax.barh([r["residue_aa"] for r in top20_loss],
            vals, left=bottom, label=ctype,
            color=col, edgecolor='white')
    bottom += np.array(vals)

ax.set_xlabel("Contacts lost on Ala mutation", fontsize=11)
ax.set_title("NSP15 Homodimer — Computational Alanine Scanning (Top 20)",
             fontsize=12, fontweight='bold')
ax.invert_yaxis()
ax.legend(fontsize=9)
plt.tight_layout()
p5 = f"{RES_DIR}/Fig5_NSP15_AlaScan_9.png"
plt.savefig(p5, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {p5}")

# ── Figure 6 — Composite ranking ──────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))
top25 = ranked[:25]
colors6 = []
for r in top25:
    if r["is_pos_ctrl"]:
        colors6.append('#C0392B')
    elif r["is_neg_ctrl"]:
        colors6.append('#7F8C8D')
    elif r["conservation"] >= 0.8:
        colors6.append('#2980B9')
    else:
        colors6.append('#85C1E9')

ax.barh([r["residue_aa"] for r in top25],
        [r["composite"] for r in top25],
        color=colors6, edgecolor='white')
ax.set_xlabel("Composite score", fontsize=11)
ax.set_title("NSP15 Homodimer — Composite Hotspot Ranking (Top 25)",
             fontsize=12, fontweight='bold')
ax.invert_yaxis()

legend_elements = [
    mpatches.Patch(color='#C0392B', label='Primary SB anchor (pos. ctrl) ★'),
    mpatches.Patch(color='#2980B9', label='Conserved hotspot (≥0.8)'),
    mpatches.Patch(color='#85C1E9', label='Variable hotspot (<0.8)'),
    mpatches.Patch(color='#7F8C8D', label='Negative control (active site)'),
]
ax.legend(handles=legend_elements, fontsize=8, loc='lower right')
plt.tight_layout()
p6 = f"{RES_DIR}/Fig6_NSP15_composite_ranking_9.png"
plt.savefig(p6, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {p6}")

# ── Save JSON ──────────────────────────────────────────────────
json_path = f"{OUT_DIR}/bsa_alascan_NSP15_9.json"
with open(json_path, "w") as f:
    json.dump({
        "complex"         : "NSP15-homodimer",
        "suffix"          : "_9",
        "top_ranked"      : ranked[:10],
        "positive_controls": {str(k):v for k,v in POS_CONTROLS.items()},
        "negative_controls": {str(k):v for k,v in NEG_CONTROLS.items()},
        "control_validation": {
            name: pos_ranks.get(pos,"not ranked")
            for pos, name in POS_CONTROLS.items()
        },
    }, f, indent=2)

print(f"Saved: {json_path}")
print("\nScript 10_9 complete")
print("Next: Script 11_9 (3D visualization)")
