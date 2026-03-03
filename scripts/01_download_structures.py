"""
Script 01: Download all PDB structures for 8 binary RTC targets
================================================================
Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

What this script does:
  1. Downloads all PDB files from RCSB Protein Data Bank
  2. Saves them to 00-reference/pdb_structures/
  3. Prints a quality summary (resolution, chains, atom count)

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/01_download_structures.py
"""

import time
import urllib.request
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────
PROJECT  = Path.home() / "projects" / "rtc-pan-coronavirus"
PDB_DIR  = PROJECT / "00-reference" / "pdb_structures"
PDB_DIR.mkdir(parents=True, exist_ok=True)

# ── All PDB IDs needed for 8 targets ──────────────────────
# Target             PDB IDs
# NSP10-NSP16     -> 6W4H, 6WVN, 7DFG
# NSP12-NSP7      -> 7BV2, 6NUR, 7C2K  (extract chains A+C)
# NSP12-NSP8      -> 7BV2, 6NUR, 7C2K  (extract chains A+B)
# NSP7-NSP8       -> 7BV2, 6NUR        (extract chains B+C)
# NSP9-NSP12      -> 9FW2, 8SQK        (chain IDs to confirm)
# NSP10-NSP14     -> 7DIY, 5C8T
# NSP13-Helicase  -> 6XEZ, 7NIO
# NSP12-NSP13     -> 7CXM, 7RDY

ALL_PDBS = [
    "6W4H", "6WVN", "7DFG",
    "7BV2", "6NUR", "7C2K",
    "9FW2", "8SQK",
    "7DIY", "5C8T",
    "6XEZ", "7NIO",
    "7CXM", "7RDY",
]


def download_pdb(pdb_id):
    out = PDB_DIR / f"{pdb_id}.pdb"
    if out.exists():
        print(f"  SKIP  {pdb_id}  (already downloaded)")
        return True
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, out)
        kb = out.stat().st_size / 1024
        print(f"  OK    {pdb_id}  ({kb:.0f} KB)")
        time.sleep(0.5)
        return True
    except Exception as e:
        print(f"  FAIL  {pdb_id}  ERROR: {e}")
        return False


def check_pdb(pdb_id):
    path = PDB_DIR / f"{pdb_id}.pdb"
    if not path.exists():
        return None
    resolution = "N/A"
    chains     = set()
    n_atoms    = 0
    with open(path) as f:
        for line in f:
            if "RESOLUTION." in line and "REMARK   2" in line:
                parts = line.split()
                for i, p in enumerate(parts):
                    if p == "RESOLUTION.":
                        try:
                            resolution = parts[i+1] + " A"
                        except:
                            pass
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chains.add(line[21])
                n_atoms += 1
    return {
        "pdb":        pdb_id,
        "resolution": resolution,
        "chains":     "".join(sorted(chains)),
        "n_atoms":    n_atoms,
    }


def main():
    print("\n" + "="*55)
    print("  Downloading PDB structures — 8 RTC targets")
    print("  Saving to: 00-reference/pdb_structures/")
    print("="*55 + "\n")

    ok = fail = 0
    for pdb_id in ALL_PDBS:
        if download_pdb(pdb_id):
            ok += 1
        else:
            fail += 1

    print(f"\n  Downloaded: {ok}  |  Failed: {fail}")

    print("\n" + "="*55)
    print("  Quality Summary")
    print("="*55)
    print(f"  {'PDB':<8} {'Resolution':<12} {'Chains':<10} {'Atoms':>7}")
    print(f"  {'-'*8} {'-'*12} {'-'*10} {'-'*7}")
    for pdb_id in ALL_PDBS:
        q = check_pdb(pdb_id)
        if q:
            print(f"  {q['pdb']:<8} "
                  f"{q['resolution']:<12} "
                  f"{q['chains']:<10} "
                  f"{q['n_atoms']:>7}")
        else:
            print(f"  {pdb_id:<8} MISSING")

    print("\n  IMPORTANT — check NSP9-NSP12:")
    print("  9FW2 and 8SQK chains need to be confirmed.")
    print("  Look at the Chains column above for these two.")
    print("="*55 + "\n")


if __name__ == "__main__":
    main()
