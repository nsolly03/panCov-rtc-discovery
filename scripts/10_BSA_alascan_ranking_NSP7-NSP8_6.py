"""
Script 10_6: BSA + Alanine Scanning + Composite Ranking — NSP7-NSP8
====================================================================
Dual mode:
  Mode A: 7BV2 crystal (Chain B=NSP8, Chain C=NSP7)
  Mode B: AF3 primary  (Chain A=NSP8, Chain B=NSP7)
Ranked separately then combined — Mode B gets 1.20 bonus
"""
import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import PDB
from Bio.PDB.SASA import ShrakeRupley
from pathlib import Path
import tempfile, os

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP7-NSP8"
RES_DIR  = PROJECT / "02-validation" / "NSP7-NSP8"

plt.rcParams.update({
    "font.family":"DejaVu Sans","font.size":11,
    "axes.spines.top":False,"axes.spines.right":False,
    "figure.dpi":150,"savefig.dpi":300,"savefig.bbox":"tight",
})
BLUE="#2C5F8A"; RED="#C0392B"; ORANGE="#E67E22"
GOLD="#F39C12"; GREEN="#27AE60"; GREY="#95A5A6"; PURPLE="#7D3C98"

HYDROPHOBIC   = {"ALA","VAL","ILE","LEU","MET","PHE","TRP","TYR","PRO"}
CHARGED_POS   = {"LYS","ARG","HIS"}
CHARGED_NEG   = {"ASP","GLU"}
HBOND_ATOMS   = {"N","O","NZ","NH1","NH2","NE","OG","OG1",
                  "OE1","OE2","ND1","ND2","NE2","OH","OD1","OD2"}

# Mode A — crystal
MODE_A = {
    "pdb":    PDB_DIR/"7BV2_NSP7-NSP8.pdb",
    "ch8":    "B", "ch7":    "C",
    "hot8":   {163,178,179,180},
    "hot7":   {24,26,27},
    "label":  "A",
    "bonus":  1.0,
}
# Mode B — AF3 primary
MODE_B = {
    "pdb":    AF3_DIR/"NSP7_NSP8_best_model.pdb",
    "ch8":    "A", "ch7":    "B",
    "hot8":   {84,87,89,90,91,92,94,95,96,98,102,
                103,106,107,110,111,116,119,120,150,190},
    "hot7":   {2,6,12,13,16,19,28,35,49,50,52,53,
                54,56,57,58,59,60,66,68,69,71,74,75,76},
    "label":  "B",
    "bonus":  1.20,   # physiological interface bonus
}
PRIMARY_SB = {("B",190), ("B",50)}   # ARG190 NSP8, GLU50 NSP7

def compute_bsa(structure, ch1, ch2):
    sr = ShrakeRupley()
    sr.compute(structure, level="R")
    sasa_c = {(c.id, r.id[1]): r.sasa
              for m in structure for c in m
              if c.id in {ch1,ch2}
              for r in c if r.id[0]==" "}
    parser = PDB.PDBParser(QUIET=True)
    sasa_m = {}
    for cid in [ch1, ch2]:
        class CS(PDB.Select):
            def accept_chain(self, c): return c.id==cid
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb",delete=False)
        io  = PDB.PDBIO(); io.set_structure(structure)
        io.save(tmp.name, CS()); tmp.close()
        ms = parser.get_structure("m", tmp.name)
        sr.compute(ms, level="R")
        for m2 in ms:
            for c2 in m2:
                for r2 in c2:
                    if r2.id[0]==" ":
                        sasa_m[(cid,r2.id[1])] = r2.sasa
        os.unlink(tmp.name)
    return {k: max(0.0, float(sasa_m.get(k,0)-v))
            for k,v in sasa_c.items()}

def count_contacts(res_a, res_b_dict):
    rn = res_a.resname.strip()
    lost = {"sb":0,"hb":0,"hy":0}
    seen_sb=set(); seen_hb=set(); seen_hy=set()
    chpos = {"LYS":["NZ"],"ARG":["NH1","NH2","NE"],"HIS":["ND1","NE2"]}
    chneg = {"ASP":["OD1","OD2"],"GLU":["OE1","OE2"]}
    for pos2,r2 in res_b_dict.items():
        rn2 = r2.resname.strip()
        # SB — residue pair level
        pa = chpos.get(rn,[]); na = chneg.get(rn,[])
        pa2= chpos.get(rn2,[]); na2= chneg.get(rn2,[])
        pairs = []
        if pa and na2: pairs += [(a,b) for a in pa for b in na2]
        if na and pa2: pairs += [(a,b) for a in na for b in pa2]
        for an1,an2 in pairs:
            try:
                d = res_a[an1] - r2[an2]
                if d<=5.0 and (pos2 not in seen_sb):
                    seen_sb.add(pos2); lost["sb"]+=1
            except: pass
        for a1 in res_a.get_atoms():
            for a2 in r2.get_atoms():
                try: d = a1-a2
                except: continue
                if (d<=3.5 and a1.name in HBOND_ATOMS
                        and a2.name in HBOND_ATOMS
                        and pos2 not in seen_hb):
                    seen_hb.add(pos2); lost["hb"]+=1
                if (d<=4.5 and rn in HYDROPHOBIC
                        and rn2 in HYDROPHOBIC
                        and a1.element=="C"
                        and a2.element=="C"
                        and pos2 not in seen_hy):
                    seen_hy.add(pos2); lost["hy"]+=1
    return lost

def process_mode(mode_cfg, cons_8, cons_7):
    parser = PDB.PDBParser(QUIET=True)
    s = parser.get_structure("x", mode_cfg["pdb"])
    ch8=mode_cfg["ch8"]; ch7=mode_cfg["ch7"]
    hot8=mode_cfg["hot8"]; hot7=mode_cfg["hot7"]
    mode_lbl=mode_cfg["label"]; bonus=mode_cfg["bonus"]

    res8 = {r.id[1]:r for m in s for c in m if c.id==ch8
            for r in c if r.id[0]==" "}
    res7 = {r.id[1]:r for m in s for c in m if c.id==ch7
            for r in c if r.id[0]==" "}

    print(f"\n    Computing BSA Mode {mode_lbl}...")
    bsa = compute_bsa(s, ch8, ch7)

    rows = []
    for chain_lbl, ch_id, hotspots, res_self, res_other, cons_df in [
        ("NSP8", ch8, hot8, res8, res7, cons_8),
        ("NSP7", ch7, hot7, res7, res8, cons_7),
    ]:
        cons_map = dict(zip(cons_df["position"],
                             cons_df["conservation"]))
        for pos in hotspots:
            if pos not in res_self: continue
            bsa_val = float(bsa.get((ch_id,pos),0))
            lost    = count_contacts(res_self[pos], res_other)
            cons    = float(cons_map.get(pos,0))
            total   = lost["sb"]*3+lost["hb"]*2+lost["hy"]
            is_sb   = ("B",pos) in PRIMARY_SB if chain_lbl=="NSP8" \
                      else ("B",pos) in {("B",50)}
            # Correct primary_sb check
            is_psb  = ((chain_lbl=="NSP8" and pos==190) or
                       (chain_lbl=="NSP7" and pos==50))
            aa = res_self[pos].resname.strip()
            rows.append({
                "residue":     f"{chain_lbl}-{aa}{pos}",
                "chain":       chain_lbl,
                "position":    pos,
                "aa":          aa,
                "mode":        mode_lbl,
                "bsa":         round(bsa_val,1),
                "sb_loss":     lost["sb"],
                "hb_loss":     lost["hb"],
                "hy_loss":     lost["hy"],
                "total_loss":  total,
                "conservation":cons,
                "primary_sb":  is_psb,
                "mode_bonus":  bonus,
            })
    return rows

def composite(df):
    mb  = df["bsa"].max(); mc = df["total_loss"].max()
    df["bsa_norm"] = df["bsa"]/mb if mb>0 else 0
    df["ct_norm"]  = df["total_loss"]/mc if mc>0 else 0
    df["composite"]= (0.35*df["bsa_norm"] +
                       0.35*df["ct_norm"]  +
                       0.20*df["conservation"] +
                       0.10*df["primary_sb"].astype(float)
                      ) * df["mode_bonus"]
    mx = df["composite"].max()
    if mx>0: df["composite"] /= mx
    return df.sort_values("composite",ascending=False)

def get_color(row):
    if row.get("primary_sb",False): return RED
    if row.get("mode","B")=="B" and row["conservation"]>=0.8:
        return BLUE if "NSP8" in row["chain"] else GREEN
    if row.get("mode","B")=="A": return PURPLE
    return GREY

def make_ylabels(df):
    lbls = []
    for _,row in df.iterrows():
        lbl = f"[{row['mode']}] {row['residue']}"
        if row.get("primary_sb",False): lbl = f"★  {lbl}"
        lbls.append(lbl)
    return lbls

patches_legend = [
    mpatches.Patch(color=RED,    label="★ Primary SB"),
    mpatches.Patch(color=BLUE,   label="NSP8 Mode B cons≥0.8"),
    mpatches.Patch(color=GREEN,  label="NSP7 Mode B cons≥0.8"),
    mpatches.Patch(color=PURPLE, label="Mode A crystal"),
    mpatches.Patch(color=GREY,   label="Variable"),
]

# ── Main ────────────────────────────────────────────────
print("\n"+"="*58)
print("  Script 10_6: BSA+AlaScan+Ranking — NSP7-NSP8")
print("  Dual mode: A=crystal, B=AF3(x1.20 bonus)")
print("="*58)

cons8 = pd.read_csv(RES_DIR/"conservation_NSP8.csv")
cons7 = pd.read_csv(RES_DIR/"conservation_NSP7.csv")

all_rows = []
for mode in [MODE_A, MODE_B]:
    all_rows.extend(process_mode(mode, cons8, cons7))

df = composite(pd.DataFrame(all_rows))
df["residue_aa"] = df["residue"]   # already has AA

print(f"\n  {'Rank':<5} {'Residue':<28} {'BSA':>7} "
      f"{'Loss':>5} {'Cons':>6} {'Score':>7} {'Mode'}")
print(f"  {'-'*68}")
for rank,(_,row) in enumerate(df.head(12).iterrows(),1):
    flags = "★" if row.get("primary_sb",False) else " "
    print(f"  {rank:<5} {row['residue']:<28} "
          f"{row['bsa']:>7.1f} "
          f"{int(row['total_loss']):>5} "
          f"{row['conservation']:>6.3f} "
          f"{row['composite']:>7.4f} "
          f"[{row['mode']}]{flags}")

# Fig4 BSA
top = df.nlargest(15,"bsa")
colors = [get_color(r) for _,r in top.iterrows()]
fig,ax = plt.subplots(figsize=(12,6))
y = np.arange(len(top))
ax.barh(y,top["bsa"],color=colors,edgecolor="white",
        linewidth=0.6,height=0.72)
for i,(_,row) in enumerate(top.iterrows()):
    ax.text(row["bsa"]+0.3,i,f"{row['bsa']:.1f} Å²",
            va="center",fontsize=8.5)
ax.set_yticks(y); ax.set_yticklabels(make_ylabels(top),fontsize=9)
ax.set_xlabel("Buried Surface Area (Å²)")
ax.set_title("NSP7-NSP8: Top 15 Residues by BSA",fontweight="bold")
ax.grid(axis="x",alpha=0.2)
ax.legend(handles=patches_legend,fontsize=8.5,loc="lower right")
plt.tight_layout()
plt.savefig(RES_DIR/"Fig4_NSP7-NSP8_BSA_6.png"); plt.close()
print(f"\n  Saved: Fig4_NSP7-NSP8_BSA_6.png")

# Fig5 AlaScan
top5 = df[df["total_loss"]>0].nlargest(15,"total_loss")
if top5.empty: top5 = df.nlargest(15,"bsa")
x = np.arange(len(top5))
fig,ax = plt.subplots(figsize=(14,5))
ax.bar(x-0.25,top5["sb_loss"]*3,0.25,
       label="Salt bridge (x3)",color=RED,alpha=0.85)
ax.bar(x,     top5["hb_loss"]*2,0.25,
       label="H-bond (x2)",color=BLUE,alpha=0.85)
ax.bar(x+0.25,top5["hy_loss"],  0.25,
       label="Hydrophobic (x1)",color=ORANGE,alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(
    [f"[{r['mode']}]{r['residue']}" for _,r in top5.iterrows()],
    rotation=45,ha="right",fontsize=8.5)
ax.set_ylabel("Estimated energy loss (AU)")
ax.set_title("NSP7-NSP8: Alanine Scanning (Dual Mode)",fontweight="bold")
ax.legend(fontsize=9); ax.grid(axis="y",alpha=0.25)
plt.tight_layout()
plt.savefig(RES_DIR/"Fig5_NSP7-NSP8_AlaScan_6.png"); plt.close()
print(f"  Saved: Fig5_NSP7-NSP8_AlaScan_6.png")

# Fig6 Ranking
top6 = df.head(15)
colors = [get_color(r) for _,r in top6.iterrows()]
fig,ax = plt.subplots(figsize=(12,6))
y = np.arange(len(top6))
ax.barh(y,top6["composite"],color=colors,edgecolor="white",
        linewidth=0.6,height=0.72)
for i,(_,row) in enumerate(top6.iterrows()):
    ax.text(row["composite"]+0.005,i,
            f"{row['composite']:.4f}",va="center",fontsize=8.5)
ax.set_yticks(y); ax.set_yticklabels(make_ylabels(top6),fontsize=9)
ax.set_xlabel("Composite Drug Target Score")
ax.set_title(
    "NSP7-NSP8: Composite Hotspot Ranking\n"
    "BSA(35%)+Contact(35%)+Cons(20%)+PrimarySB(10%) | Mode B x1.20",
    fontweight="bold")
ax.set_xlim(0,1.12); ax.grid(axis="x",alpha=0.2)
ax.legend(handles=patches_legend,fontsize=8.5,loc="lower right")
plt.tight_layout()
plt.savefig(RES_DIR/"Fig6_NSP7-NSP8_composite_ranking_6.png"); plt.close()
print(f"  Saved: Fig6_NSP7-NSP8_composite_ranking_6.png")

# Save CSV + JSON
df.to_csv(RES_DIR/"composite_ranking_NSP7-NSP8_6.csv",index=False)
print(f"  Saved: composite_ranking_NSP7-NSP8_6.csv")

top_pharm = df.iloc[0]
print(f"\n  Primary pharmacophore: {top_pharm['residue']} "
      f"[Mode {top_pharm['mode']}] "
      f"score={top_pharm['composite']:.4f}")
print("="*58+"\n")
