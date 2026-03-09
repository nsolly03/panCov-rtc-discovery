#!/usr/bin/env python3
"""
Script 12 - ZINC20 Drug-Like Subset Download + ADMET Filtering
Pan-coronavirus RTC Inhibitor Discovery Pipeline
"""

import sys, json, time
import urllib.request, urllib.error
from pathlib import Path
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import Descriptors, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

BASE    = Path(__file__).parent.parent
DATA    = BASE / "data" / "zinc20"
DATA.mkdir(parents=True, exist_ok=True)
OUT_SMI = DATA / "zinc20_filtered.smi"
OUT_SUM = DATA / "admet_filter_summary.json"

ZINC20_BASE = "https://files.docking.org/2D"
HEADERS     = {"User-Agent": "Mozilla/5.0 (academic screening pipeline)"}

# Sub-tranche suffix letters
LETTERS = "ABCDEFGHIJ"

# TEST: 2 tranches x 3 sub-tranches = ~8,000 compounds
TEST_TRANCHES = ["BA", "CA"]
TEST_SUBS     = ["AA", "AB", "AC"]

# FULL: all drug-like tranches
ALL_TRANCHES = [
    "BA","BB","BC","BD","BE",
    "CA","CB","CC","CD","CE",
    "DA","DB","DC","DD","DE",
    "EA","EB","EC","ED","EE",
    "FA","FB","FC","FD","FE",
    "GA","GB","GC","GD","GE",
    "HA","HB","HC","HD","HE",
]
ALL_SUBS = [a+b for a in LETTERS for b in LETTERS]  # AA..JJ = 100 sub-tranches

def setup_pains():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    return FilterCatalog.FilterCatalog(params)

def passes_admet(mol, pains):
    if mol is None:
        return False, "invalid_mol"
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = Descriptors.NumHDonors(mol)
    hba  = Descriptors.NumHAcceptors(mol)
    rb   = Descriptors.NumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)
    if not (150 <= mw  <= 500): return False, "MW"
    if logp > 5:                 return False, "LogP"
    if hbd > 5:                 return False, "HBD"
    if hba > 10:                return False, "HBA"
    if rb  > 10:                return False, "RotBond"
    if tpsa > 140:              return False, "TPSA"
    if pains.HasMatch(mol):     return False, "PAINS"
    return True, "PASS"

def fetch_subtranche(tranche, sub, retries=3):
    """Download one sub-tranche file e.g. BAAA.smi"""
    code  = tranche + sub          # e.g. BAAA
    url   = f"{ZINC20_BASE}/{tranche}/{code}.smi"
    local = DATA / f"{code}.smi"

    if local.exists():
        pass  # use cache silently
    else:
        req = urllib.request.Request(url, headers=HEADERS)
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(req, timeout=30) as resp:
                    local.write_bytes(resp.read())
                break
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    return []   # sub-tranche doesn't exist — normal
                print(f"    [WARN] {code} attempt {attempt+1}: HTTP {e.code}")
                time.sleep(2)
                if attempt == retries-1:
                    return []
            except Exception as e:
                print(f"    [WARN] {code} attempt {attempt+1}: {e}")
                time.sleep(2)
                if attempt == retries-1:
                    return []

    compounds = []
    try:
        with open(local) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('smiles') or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    compounds.append((parts[0], parts[1]))
                elif len(parts) == 1:
                    compounds.append((parts[0], "UNKNOWN"))
    except Exception as e:
        print(f"    [ERROR] parsing {code}: {e}")
    return compounds

def main():
    print("="*65)
    print("Script 12 - ZINC20 Download + ADMET Filtering")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*65)

    USE_TEST = "--test" in sys.argv
    tranches = TEST_TRANCHES if USE_TEST else ALL_TRANCHES
    subs     = TEST_SUBS     if USE_TEST else ALL_SUBS
    mode     = "TEST" if USE_TEST else "FULL"
    print(f"Mode: {mode} | Tranches: {len(tranches)} | Sub-tranches each: {len(subs)}")

    pains = setup_pains()
    stats = {
        "mode": mode,
        "total_input": 0, "total_passed": 0, "total_failed": 0,
        "fail_reasons": {},
        "started": datetime.now().isoformat(),
    }
    passed = []

    for i, tranche in enumerate(tranches):
        print(f"\n[{i+1}/{len(tranches)}] Tranche {tranche}...")
        tp = tf = sub_count = 0
        for sub in subs:
            compounds = fetch_subtranche(tranche, sub)
            if not compounds:
                continue
            sub_count += 1
            for smi, zid in compounds:
                stats["total_input"] += 1
                mol = Chem.MolFromSmiles(smi)
                ok, reason = passes_admet(mol, pains)
                if ok:
                    passed.append((smi, zid))
                    stats["total_passed"] += 1
                    tp += 1
                else:
                    stats["total_failed"] += 1
                    stats["fail_reasons"][reason] = stats["fail_reasons"].get(reason,0)+1
                    tf += 1
        rate = tp/(tp+tf)*100 if (tp+tf)>0 else 0
        print(f"  Sub-tranches: {sub_count} | Input:{tp+tf:,} | "
              f"Passed:{tp:,} ({rate:.1f}%)")

    print(f"\nWriting {len(passed):,} compounds to {OUT_SMI}...")
    with open(OUT_SMI, 'w') as f:
        f.write("# ZINC20 filtered | Lipinski+PAINS | panCov RTC pipeline\n")
        for smi, zid in passed:
            f.write(f"{smi} {zid}\n")

    stats["completed"]     = datetime.now().isoformat()
    stats["output_file"]   = str(OUT_SMI)
    stats["total_compounds_written"] = len(passed)
    stats["pass_rate_pct"] = round(
        stats["total_passed"]/stats["total_input"]*100, 2
    ) if stats["total_input"] > 0 else 0

    with open(OUT_SUM, 'w') as f:
        json.dump(stats, f, indent=2)

    print("\n"+"="*65)
    print("SUMMARY")
    print("="*65)
    print(f"Total input:   {stats['total_input']:>8,}")
    print(f"Total passed:  {stats['total_passed']:>8,}  ({stats['pass_rate_pct']}%)")
    print(f"Total failed:  {stats['total_failed']:>8,}")
    print("\nFail reasons:")
    for r,c in sorted(stats["fail_reasons"].items(), key=lambda x:-x[1]):
        print(f"  {r:<20} {c:>8,}")
    print(f"\nOutput: {OUT_SMI}")
    print("Status: COMPLETE")

if __name__ == "__main__":
    main()
