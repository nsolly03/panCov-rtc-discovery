import json, subprocess, shutil, numpy as np
from pathlib import Path
from Bio import PDB

PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
AF3_DIR  = PROJECT / "01-alphafold3" / "NSP7-NSP8"
RES_DIR  = PROJECT / "02-validation" / "NSP7-NSP8"
TMP_DIR  = PROJECT / "tmp"
TMP_DIR.mkdir(exist_ok=True)

PDB_FILES = {
    "7BV2": PDB_DIR / "7BV2_NSP7-NSP8.pdb",
    "6NUR": PDB_DIR / "6NUR_NSP7-NSP8.pdb",
    "AF3":  AF3_DIR / "NSP7_NSP8_best_model.pdb",
}
CHAIN_NSP8_PDB = "B"; CHAIN_NSP7_PDB = "C"
CHAIN_NSP8_AF3 = "A"; CHAIN_NSP7_AF3 = "B"
FPOCKET_GATE   = 0.50
MODE_A_NSP8 = {163,178,179,180}; MODE_A_NSP7 = {24,26,27}
MODE_B_NSP8 = {87,91,92,94,98,110,111,116,120,150,190}
MODE_B_NSP7 = {49,50,56}

def run_fpocket(pdb_file):
    subprocess.run(["fpocket","-f",str(pdb_file)],
        capture_output=True, text=True, cwd=str(TMP_DIR))

def parse_fpocket(pdb_file):
    stem = pdb_file.stem
    for d in [TMP_DIR/f"{stem}_out", pdb_file.parent/f"{stem}_out"]:
        info = d / f"{stem}_info.txt"
        if info.exists():
            pockets, cur = [], {}
            for line in open(info):
                line = line.strip()
                if line.startswith("Pocket"):
                    if cur: pockets.append(cur)
                    cur = {"id": int(line.split()[1])}
                elif ":" in line:
                    k,_,v = line.partition(":")
                    try: cur[k.strip()] = float(v.strip())
                    except: cur[k.strip()] = v.strip()
            if cur: pockets.append(cur)
            return pockets
    return []

def define_box(structure, c8, c7, h8, h7, pad=5.0):
    coords = []
    for m in structure:
        for c in m:
            hs = h8 if c.id==c8 else h7 if c.id==c7 else None
            if hs:
                for r in c:
                    if r.id[0]==" " and r.id[1] in hs:
                        coords.extend([a.get_vector().get_array()
                                        for a in r.get_atoms()])
    if not coords: return None
    coords = np.array(coords)
    center = coords.mean(axis=0)
    size   = coords.max(axis=0)-coords.min(axis=0)+2*pad
    return {k:round(float(v),3) for k,v in zip(
        ["center_x","center_y","center_z",
         "size_x","size_y","size_z"],
        list(center)+list(size))}

print("\n"+"="*58)
print("  Script 07_6: Pocket Detection — NSP7-NSP8")
print(f"  Extra gate: fpocket >= {FPOCKET_GATE}")
print("="*58+"\n")

parser = PDB.PDBParser(QUIET=True)
print("  1. fpocket analysis...")
all_pockets = {}; best_scores = {}
for sid, pdb_f in PDB_FILES.items():
    if not pdb_f.exists():
        print(f"     {sid}: not found — skip"); continue
    run_fpocket(pdb_f)
    pockets = parse_fpocket(pdb_f)
    all_pockets[sid] = pockets
    if pockets:
        best  = max(pockets, key=lambda p:
            p.get("Druggability Score", p.get("Drug Score",0)))
        score = best.get("Druggability Score",
                          best.get("Drug Score",0))
        vol   = best.get("Volume",0)
        best_scores[sid] = score
        gate  = "✅" if score >= FPOCKET_GATE else "❌"
        print(f"     {sid}: {len(pockets)} pockets | "
              f"best={score:.3f} vol={vol:.1f} A3 {gate}")
    else:
        best_scores[sid]=0
        print(f"     {sid}: no pockets parsed")

p7b = all_pockets.get("7BV2",[])
if p7b:
    print(f"\n     7BV2 top 5 pockets:")
    print(f"     {'Pocket':<8} {'DrugScore':>10} {'Volume':>10}")
    print(f"     {'-'*30}")
    for p in sorted(p7b, key=lambda x:x.get("Druggability Score",
            x.get("Drug Score",0)), reverse=True)[:5]:
        print(f"     {p.get('id','?'):<8} "
              f"{p.get('Druggability Score',p.get('Drug Score',0)):>10.3f} "
              f"{p.get('Volume',0):>10.1f}")

max_score = max(best_scores.values()) if best_scores else 0
gate_pass = max_score >= FPOCKET_GATE
print(f"\n  fpocket gate: max={max_score:.3f} "
      f"{'✅ PASS' if gate_pass else '❌ FAIL'}")

print("\n  2. Mode A docking box (crystal interface)...")
s7bv2 = parser.get_structure("7BV2", PDB_FILES["7BV2"])
box_A = define_box(s7bv2, CHAIN_NSP8_PDB, CHAIN_NSP7_PDB,
                    MODE_A_NSP8, MODE_A_NSP7)
if box_A:
    vol_A = box_A["size_x"]*box_A["size_y"]*box_A["size_z"]
    print(f"     Center: ({box_A['center_x']}, "
          f"{box_A['center_y']}, {box_A['center_z']})")
    print(f"     Size:   {box_A['size_x']} x "
          f"{box_A['size_y']} x {box_A['size_z']} A")
    print(f"     Volume: {vol_A:.0f} A3")

print("\n  3. Mode B docking box (AF3 interface)...")
saf3  = parser.get_structure("AF3", PDB_FILES["AF3"])
box_B = define_box(saf3, CHAIN_NSP8_AF3, CHAIN_NSP7_AF3,
                    MODE_B_NSP8, MODE_B_NSP7)
if box_B:
    vol_B = box_B["size_x"]*box_B["size_y"]*box_B["size_z"]
    print(f"     Center: ({box_B['center_x']}, "
          f"{box_B['center_y']}, {box_B['center_z']})")
    print(f"     Size:   {box_B['size_x']} x "
          f"{box_B['size_y']} x {box_B['size_z']} A")
    print(f"     Volume: {vol_B:.0f} A3")

out = RES_DIR / "pocket_analysis_6.json"
json.dump({
    "complex":"NSP7-NSP8","fpocket_gate":FPOCKET_GATE,
    "fpocket_pass":gate_pass,"best_druggability":max_score,
    "fpocket_scores":best_scores,
    "docking_box_A":box_A,"docking_box_B":box_B,
    "mode_A_hotspots":{"NSP8":sorted(MODE_A_NSP8),
                        "NSP7":sorted(MODE_A_NSP7)},
    "mode_B_hotspots":{"NSP8":sorted(MODE_B_NSP8),
                        "NSP7":sorted(MODE_B_NSP7)},
}, open(out,"w"), indent=2)
print(f"\n  Saved: pocket_analysis_6.json")
print("="*58+"\n")
