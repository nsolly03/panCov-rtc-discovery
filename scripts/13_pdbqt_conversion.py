#!/usr/bin/env python3
"""
Script 13 - PDBQT Conversion for AutoDock Vina
Pan-coronavirus RTC Inhibitor Discovery Pipeline
Converts filtered ZINC20 SMILES -> 3D PDBQT files using Meeko
"""

import sys, json, subprocess
from pathlib import Path
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

BASE     = Path(__file__).parent.parent
DATA     = BASE / "data" / "zinc20"
PDBQT    = BASE / "data" / "pdbqt"
PDBQT.mkdir(parents=True, exist_ok=True)
OUT_SUM  = DATA / "pdbqt_conversion_summary.json"

IN_SMI   = DATA / "zinc20_filtered.smi"

def smiles_to_pdbqt(smiles, zinc_id, out_dir):
    """Convert SMILES to 3D PDBQT via RDKit + Meeko."""
    out_path = out_dir / f"{zinc_id}.pdbqt"
    if out_path.exists():
        return "cached"
    try:
        # Generate 3D conformer
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "invalid_smiles"
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result == -1:
            return "embed_failed"
        AllChem.MMFFOptimizeMolecule(mol)

        # Meeko preparation
        preparator = MoleculePreparation()
        mol_setups  = preparator.prepare(mol)
        if not mol_setups:
            return "meeko_failed"

        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
        if not is_ok:
            return f"pdbqt_error: {error_msg}"

        out_path.write_text(pdbqt_string)
        return "ok"

    except Exception as e:
        return f"exception: {str(e)[:60]}"

def main():
    print("="*65)
    print("Script 13 - PDBQT Conversion (SMILES -> AutoDock Vina)")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*65)

    USE_TEST = "--test" in sys.argv

    if not IN_SMI.exists():
        print(f"ERROR: {IN_SMI} not found. Run Script 12 first.")
        sys.exit(1)

    # Read filtered SMILES
    compounds = []
    with open(IN_SMI) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                compounds.append((parts[0], parts[1]))

    if USE_TEST:
        compounds = compounds[:200]
        print(f"TEST mode: converting first 200 compounds")
    else:
        print(f"FULL mode: converting {len(compounds):,} compounds")

    stats = {
        "total_input": len(compounds),
        "ok": 0, "cached": 0, "failed": 0,
        "fail_reasons": {},
        "started": datetime.now().isoformat(),
    }

    for i, (smi, zid) in enumerate(compounds):
        result = smiles_to_pdbqt(smi, zid, PDBQT)
        if result == "ok":
            stats["ok"] += 1
        elif result == "cached":
            stats["cached"] += 1
        else:
            stats["failed"] += 1
            stats["fail_reasons"][result] = stats["fail_reasons"].get(result, 0) + 1

        if (i+1) % 50 == 0 or (i+1) == len(compounds):
            done = stats["ok"] + stats["cached"]
            print(f"  [{i+1}/{len(compounds)}] OK:{done} Failed:{stats['failed']}")

    stats["completed"]  = datetime.now().isoformat()
    stats["output_dir"] = str(PDBQT)
    stats["total_pdbqt_files"] = len(list(PDBQT.glob("*.pdbqt")))

    with open(OUT_SUM, 'w') as f:
        json.dump(stats, f, indent=2)

    print("\n"+"="*65)
    print("SUMMARY")
    print("="*65)
    print(f"Total input:  {stats['total_input']:>6,}")
    print(f"Converted:    {stats['ok']:>6,}")
    print(f"Cached:       {stats['cached']:>6,}")
    print(f"Failed:       {stats['failed']:>6,}")
    if stats["fail_reasons"]:
        print("\nFail reasons:")
        for r,c in sorted(stats["fail_reasons"].items(), key=lambda x:-x[1]):
            print(f"  {r:<35} {c:>5,}")
    print(f"\nPDBQT files: {stats['total_pdbqt_files']:,} in {PDBQT}")
    print("Status: COMPLETE")

if __name__ == "__main__":
    main()
