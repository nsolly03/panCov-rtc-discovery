"""
Script 08_6: Docking Preparation — NSP7-NSP8
=============================================
Two receptor files:
  Mode A: 7BV2 (Crystal — for completeness)
  Mode B: AF3  (Primary druggable target)
"""
import json
from pathlib import Path
from Bio import PDB

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP7-NSP8"
RES_DIR  = PROJECT / "02-validation" / "NSP7-NSP8"
DOCK_DIR = PROJECT / "03-virtual-screening" / "NSP7-NSP8_6"
DOCK_DIR.mkdir(parents=True, exist_ok=True)

CHAIN_NSP8_PDB = "B"; CHAIN_NSP7_PDB = "C"
CHAIN_NSP8_AF3 = "A"; CHAIN_NSP7_AF3 = "B"

MODE_A_NSP8 = {163,178,179,180}; MODE_A_NSP7 = {24,26,27}
MODE_B_NSP8 = {87,91,92,94,98,110,111,116,120,150,190}
MODE_B_NSP7 = {49,50,56}

HOTSPOTS_NSP8_B = [84,87,89,90,91,92,94,95,96,98,
                    102,103,106,107,110,111,116,119,
                    120,150,190]
HOTSPOTS_NSP7_B = [2,6,12,13,16,19,28,35,49,50,
                    52,53,54,56,57,58,59,60,66,68,
                    69,71,74,75,76]

class Select2Chain(PDB.Select):
    def __init__(self, chains):
        self.chains = chains
    def accept_chain(self, c):
        return c.id in self.chains
    def accept_residue(self, r):
        return r.id[0]==" " and r.resname.strip()!="HOH"
    def accept_atom(self, a):
        return a.altloc in (" ","A")

def write_receptor(structure, chains, out_path):
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(str(out_path), Select2Chain(chains))
    # Count
    p  = PDB.PDBParser(QUIET=True)
    s2 = p.get_structure("r", out_path)
    res = sum(1 for m in s2 for c in m
              for r in c if r.id[0]==" ")
    ats = sum(1 for m in s2 for c in m
              for r in c if r.id[0]==" "
              for _ in r.get_atoms())
    return res, ats

def write_vina(rec_path, box, out_path):
    with open(out_path,"w") as f:
        f.write(f"receptor = {rec_path}\n\n")
        f.write(f"center_x = {box['center_x']}\n")
        f.write(f"center_y = {box['center_y']}\n")
        f.write(f"center_z = {box['center_z']}\n\n")
        f.write(f"size_x = {box['size_x']}\n")
        f.write(f"size_y = {box['size_y']}\n")
        f.write(f"size_z = {box['size_z']}\n\n")
        f.write("exhaustiveness = 16\n")
        f.write("num_modes = 9\n")
        f.write("energy_range = 3\n")

print("\n"+"="*58)
print("  Script 08_6: Docking Prep — NSP7-NSP8")
print("  Mode A: 7BV2 crystal | Mode B: AF3 primary")
print("="*58+"\n")

parser = PDB.PDBParser(QUIET=True)
pocket = json.load(open(RES_DIR/"pocket_analysis_6.json"))
box_A  = pocket["docking_box_A"]
box_B  = pocket["docking_box_B"]

# Verify hotspots
def verify(structure, chain, hotspots, label):
    res_map = {r.id[1]: r.resname.strip()
               for m in structure for c in m
               if c.id==chain
               for r in c if r.id[0]==" "}
    missing = [p for p in hotspots if p not in res_map]
    print(f"    {label}: {len(hotspots)-len(missing)}/"
          f"{len(hotspots)} hotspots present "
          f"{'✅' if not missing else '⚠️ missing:'+str(missing)}")
    return res_map

# ── Mode A: 7BV2 ──────────────────────────────────────
print("  Mode A: 7BV2 crystal receptor...")
s7bv2  = parser.get_structure("7BV2",
    PDB_DIR/"7BV2_NSP7-NSP8.pdb")
rec_A  = DOCK_DIR/"receptor_NSP7-NSP8_ModeA_6.pdb"
res,ats = write_receptor(s7bv2,
    {CHAIN_NSP8_PDB, CHAIN_NSP7_PDB}, rec_A)
print(f"    Receptor: {res} residues, {ats} atoms")
verify(s7bv2, CHAIN_NSP8_PDB, sorted(MODE_A_NSP8), "NSP8 Mode A")
verify(s7bv2, CHAIN_NSP7_PDB, sorted(MODE_A_NSP7), "NSP7 Mode A")
vina_A = DOCK_DIR/"vina_config_NSP7-NSP8_ModeA_6.txt"
write_vina(rec_A, box_A, vina_A)
print(f"    Vina config: {vina_A.name}")

# ── Mode B: AF3 ───────────────────────────────────────
print("\n  Mode B: AF3 receptor (PRIMARY)...")
saf3   = parser.get_structure("AF3",
    AF3_DIR/"NSP7_NSP8_best_model.pdb")
rec_B  = DOCK_DIR/"receptor_NSP7-NSP8_ModeB_AF3_6.pdb"
res,ats = write_receptor(saf3,
    {CHAIN_NSP8_AF3, CHAIN_NSP7_AF3}, rec_B)
print(f"    Receptor: {res} residues, {ats} atoms")
verify(saf3, CHAIN_NSP8_AF3, HOTSPOTS_NSP8_B, "NSP8 Mode B")
verify(saf3, CHAIN_NSP7_AF3, HOTSPOTS_NSP7_B, "NSP7 Mode B")
vina_B = DOCK_DIR/"vina_config_NSP7-NSP8_ModeB_AF3_6.txt"
write_vina(rec_B, box_B, vina_B)
print(f"    Vina config: {vina_B.name}")

# ── VirtualFlow config (Mode B primary) ──────────────
vf_cfg = DOCK_DIR/"virtualflow_config_NSP7-NSP8_6.json"
json.dump({
    "target":           "NSP7-NSP8",
    "primary_mode":     "B_AF3",
    "receptor_ModeA":   str(rec_A),
    "receptor_ModeB":   str(rec_B),
    "center_ModeB":     [box_B["center_x"],
                          box_B["center_y"],
                          box_B["center_z"]],
    "size_ModeB":       [box_B["size_x"],
                          box_B["size_y"],
                          box_B["size_z"]],
    "exhaustiveness":    16,
    "fpocket_gate_pass": True,
    "fpocket_ModeB":     0.531,
    "dual_mode":         True,
    "note":              "Mode A crystal undruggable "
                         "(score=0.276). Mode B AF3 "
                         "primary target (score=0.531). "
                         "PHE92/LEU91/ARG190 pan-cov core.",
    "pan_cov_NSP8":     [87,91,92,94,98,110,111,
                          116,120,150,190],
    "primary_sb":        "ARG190-GLU50",
}, open(vf_cfg,"w"), indent=2)
print(f"\n  VF config: {vf_cfg.name}")
print(f"\n  Output: 03-virtual-screening/NSP7-NSP8_6/")
print("="*58+"\n")
