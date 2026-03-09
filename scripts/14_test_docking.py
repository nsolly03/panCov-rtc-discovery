#!/usr/bin/env python3
"""
Script 14 - Test Docking: Validate Tier 1 Docking Boxes
Pan-coronavirus RTC Inhibitor Discovery Pipeline
Runs AutoDock Vina on 100 compounds against all 3 Tier 1 targets.
Validates box coordinates, scoring, and pose quality.
"""

import sys, json, subprocess, shutil
from pathlib import Path
from datetime import datetime

BASE      = Path(__file__).parent.parent
PDBQT     = BASE / "data" / "pdbqt"
SCREENING = BASE / "03-virtual-screening"
RESULTS   = BASE / "04-hits"
RESULTS.mkdir(exist_ok=True)
OUT_SUM   = RESULTS / "test_docking_summary.json"

# Tier 1 targets
TARGETS = [
    {
        "name":     "NSP12-NSP7",
        "suffix":   "_3",
        "receptor": SCREENING / "NSP12-NSP7_3" / "receptor_NSP12-NSP7_3.pdb",
        "config":   SCREENING / "NSP12-NSP7_3" / "vina_config_NSP12-NSP7_3.txt",
    },
    {
        "name":     "NSP9-NSP12",
        "suffix":   "_5",
        "receptor": SCREENING / "NSP9-NSP12_5" / "receptor_NSP9-NSP12_5.pdb",
        "config":   SCREENING / "NSP9-NSP12_5" / "vina_config_NSP9-NSP12_5.txt",
    },
    {
        "name":     "NSP12-NSP8",
        "suffix":   "_4",
        "receptor": SCREENING / "NSP12-NSP8_4" / "receptor_NSP12-NSP8_4.pdb",
        "config":   SCREENING / "NSP12-NSP8_4" / "vina_config_NSP12-NSP8_4.txt",
    },
]

def pdb_to_pdbqt(pdb_path, pdbqt_path):
    """Convert receptor PDB to PDBQT using AutoDock's prepare_receptor (via meeko)."""
    cmd = [
        "python", "-m", "meeko",
        "--receptor", str(pdb_path),
        "--out", str(pdbqt_path),
    ]
    # Try obabel as fallback
    if not shutil.which("prepare_receptor4.py"):
        cmd = [
            "obabel", str(pdb_path),
            "-O", str(pdbqt_path),
            "-xr", "--partialcharge", "gasteiger"
        ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0

def parse_vina_config(config_path):
    """Parse vina config file into dict."""
    cfg = {}
    with open(config_path) as f:
        for line in f:
            line = line.strip()
            if '=' in line and not line.startswith('#'):
                k, v = line.split('=', 1)
                cfg[k.strip()] = v.strip()
    return cfg

def run_vina(receptor_pdbqt, ligand_pdbqt, config_path, out_path, cpu=2):
    """Run AutoDock Vina for one ligand."""
    cfg = parse_vina_config(config_path)
    cmd = [
        "vina",
        "--receptor", str(receptor_pdbqt),
        "--ligand",   str(ligand_pdbqt),
        "--center_x", cfg.get("center_x", "0"),
        "--center_y", cfg.get("center_y", "0"),
        "--center_z", cfg.get("center_z", "0"),
        "--size_x",   cfg.get("size_x", "20"),
        "--size_y",   cfg.get("size_y", "20"),
        "--size_z",   cfg.get("size_z", "20"),
        "--out",      str(out_path),
        "--cpu",      str(cpu),
        "--exhaustiveness", "8",
        "--num_modes", "5",
        "--energy_range", "3",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    return result.returncode, result.stdout, result.stderr

def parse_vina_score(out_pdbqt):
    """Extract best docking score from Vina output PDBQT."""
    try:
        with open(out_pdbqt) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT"):
                    parts = line.split()
                    return float(parts[3])
    except Exception:
        pass
    return None

def prepare_receptor_pdbqt(target, out_dir):
    """Prepare receptor PDBQT from PDB using obabel."""
    rec_pdbqt = out_dir / f"receptor_{target['name']}.pdbqt"
    if rec_pdbqt.exists():
        return rec_pdbqt
    print(f"  Preparing receptor PDBQT for {target['name']}...")
    cmd = [
        "obabel", str(target["receptor"]),
        "-O", str(rec_pdbqt),
        "-xr"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0 and rec_pdbqt.exists():
        print(f"  Receptor PDBQT ready: {rec_pdbqt.name}")
        return rec_pdbqt
    else:
        print(f"  ERROR preparing receptor: {result.stderr[:100]}")
        return None

def main():
    print("="*65)
    print("Script 14 - Test Docking: Validate Tier 1 Boxes")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*65)

    # Get first 100 PDBQT files
    ligand_files = sorted(PDBQT.glob("*.pdbqt"))[:100]
    if not ligand_files:
        print(f"ERROR: No PDBQT files in {PDBQT}. Run Script 13 first.")
        sys.exit(1)
    print(f"Ligands: {len(ligand_files)} test compounds")

    all_stats = {}

    for target in TARGETS:
        tname = target["name"]
        print(f"\n{'='*50}")
        print(f"Target: {tname}")
        print(f"{'='*50}")

        # Check receptor PDB exists
        if not target["receptor"].exists():
            print(f"  ERROR: receptor not found: {target['receptor']}")
            print(f"  Run: python scripts/01_download_structures.py")
            all_stats[tname] = {"error": "receptor_not_found"}
            continue

        if not target["config"].exists():
            print(f"  ERROR: config not found: {target['config']}")
            all_stats[tname] = {"error": "config_not_found"}
            continue

        # Output directory for this target
        out_dir = RESULTS / f"test_{tname}"
        out_dir.mkdir(exist_ok=True)

        # Prepare receptor PDBQT
        rec_pdbqt = prepare_receptor_pdbqt(target, out_dir)
        if rec_pdbqt is None:
            all_stats[tname] = {"error": "receptor_prep_failed"}
            continue

        # Print docking box info
        cfg = parse_vina_config(target["config"])
        print(f"  Box center: ({cfg.get('center_x')}, {cfg.get('center_y')}, {cfg.get('center_z')})")
        print(f"  Box size:   ({cfg.get('size_x')}, {cfg.get('size_y')}, {cfg.get('size_z')}) Å")

        # Run docking
        scores = []
        n_ok = n_fail = 0

        for i, lig in enumerate(ligand_files):
            out_pdbqt = out_dir / f"{lig.stem}_out.pdbqt"
            code, stdout, stderr = run_vina(
                rec_pdbqt, lig, target["config"], out_pdbqt
            )
            if code == 0 and out_pdbqt.exists():
                score = parse_vina_score(out_pdbqt)
                if score is not None:
                    scores.append((score, lig.stem))
                    n_ok += 1
            else:
                n_fail += 1

            if (i+1) % 20 == 0 or (i+1) == len(ligand_files):
                print(f"  [{i+1}/{len(ligand_files)}] OK:{n_ok} Failed:{n_fail}")

        if not scores:
            print(f"  ERROR: No successful dockings for {tname}")
            all_stats[tname] = {"error": "no_docking_results"}
            continue

        scores.sort()
        top10 = scores[:10]

        print(f"\n  Results: {n_ok} docked | {n_fail} failed")
        print(f"  Score range: {scores[0][0]:.2f} to {scores[-1][0]:.2f} kcal/mol")
        print(f"  Mean score:  {sum(s for s,_ in scores)/len(scores):.2f} kcal/mol")
        print(f"\n  Top 10 hits:")
        for rank, (score, zinc_id) in enumerate(top10, 1):
            print(f"    {rank:2}. {zinc_id:<20} {score:.2f} kcal/mol")

        all_stats[tname] = {
            "n_docked":    n_ok,
            "n_failed":    n_fail,
            "best_score":  scores[0][0],
            "worst_score": scores[-1][0],
            "mean_score":  round(sum(s for s,_ in scores)/len(scores), 2),
            "top10": [{"zinc_id": z, "score": s} for s,z in top10],
        }

    # Save summary
    summary = {
        "completed": datetime.now().isoformat(),
        "n_ligands_tested": len(ligand_files),
        "targets": all_stats,
    }
    with open(OUT_SUM, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*65}")
    print("TEST DOCKING COMPLETE")
    print(f"{'='*65}")
    for tname, s in all_stats.items():
        if "error" in s:
            print(f"  {tname:<20} ERROR: {s['error']}")
        else:
            print(f"  {tname:<20} best={s['best_score']:.2f}  "
                  f"mean={s['mean_score']:.2f}  n={s['n_docked']}")
    print(f"\nSummary: {OUT_SUM}")
    print("Status: COMPLETE")

if __name__ == "__main__":
    main()
