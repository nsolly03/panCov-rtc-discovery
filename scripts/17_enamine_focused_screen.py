#!/usr/bin/env python3
"""
Script 17 — ZINC20 Similarity-Focused Expansion Screen
Pan-coronavirus RTC Inhibitor Discovery Pipeline
Olivier Nsekuye | GIGA-VIN Lab, University of Liège | 2026-03-10

Downloads additional ZINC20 tranches, filters by:
  1. Tanimoto similarity >= 0.3 to any of 11 triple-target hits
  2. Lipinski RO5 + ADMET (same as Script 12)
  3. PAINS filters
Then converts to PDBQT for NIC5 SLURM screen.
"""

import os, sys, json, time, requests
from pathlib import Path
from datetime import datetime

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
import meeko

BASE    = Path(__file__).parent.parent
ZINC20  = BASE / "data" / "zinc20"
PDBQT   = BASE / "data" / "pdbqt_enamine"
HITS    = BASE / "04-hits"
FIGS    = BASE / "figures" / "script17"
ZINC20.mkdir(parents=True, exist_ok=True)
PDBQT.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

# ── 11 triple-target hits SMILES ────────────────────────────────────────────
LEAD_SMILES = [
    ('351017',   'O=C1c2ccccc2[C@]2(O)Nc3[nH]c(=O)[nH]c(=O)c3[C@]12O'),
    ('13633807', 'O=C1NC(=O)C2=Nn3c(nc4ccccc4c3=O)C[C@]2(O)N1'),
    ('13633805', 'O=C1NC(=O)C2=Nn3c(nc4ccccc4c3=O)C[C@@]2(O)N1'),
    ('351016',   'O=C1c2ccccc2[C@@]2(O)Nc3[nH]c(=O)[nH]c(=O)c3[C@]12O'),
    ('5024943',  'CC1(C)C[C@@]2(CCO1)NNC(=O)[C@@H]1C(=O)NNC(=O)[C@@H]12'),
    ('351018',   'O=C1c2ccccc2[C@@]2(O)Nc3[nH]c(=O)[nH]c(=O)c3[C@@]12O'),
    ('3169307',  'NC(NC(=O)c1ccccc1)=C1C(=O)NC(=O)NC1=O'),
    ('5024944',  'CC1(C)C[C@]2(CCO1)NNC(=O)[C@@H]1C(=O)NNC(=O)[C@H]12'),
    ('351019',   'O=C1c2ccccc2[C@]2(O)Nc3[nH]c(=O)[nH]c(=O)c3[C@@]12O'),
    ('13662104', 'CC1(C)C[C@]2(CCO1)NNC(=O)[C@H]1C(=O)NNC(=O)[C@@H]12'),
    ('5024945',  'CC1(C)C[C@@]2(CCO1)NNC(=O)[C@@H]1C(=O)NNC(=O)[C@H]12'),
]

SIM_THRESHOLD = 0.30   # Tanimoto similarity cutoff
MAX_COMPOUNDS = 10000  # cap for NIC5 screen

# ── Tranches to download (skip already-done BA/CA sub AA-AC) ────────────────
HEADERS = {'User-Agent': 'Mozilla/5.0'}
ALL_TRANCHES = ["BA","BB","BC","BD","BE",
                "CA","CB","CC","CD","CE",
                "DA","DB","DC","DD","DE",
                "EA","EB","EC","ED","EE"]
LETTERS = "ABCDEFGHIJ"
ALL_SUBS = [a+b for a in LETTERS for b in LETTERS]

# Already downloaded
DONE = set(p.stem for p in ZINC20.glob("*.smi") if p.stem != "zinc20_filtered")

def setup_pains():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    return FilterCatalog.FilterCatalog(params)

def passes_admet(mol):
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = Descriptors.NumHDonors(mol)
    hba  = Descriptors.NumHAcceptors(mol)
    rb   = Descriptors.NumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)
    return (150 <= mw <= 500 and logp <= 5 and
            hbd <= 5 and hba <= 10 and rb <= 10 and tpsa <= 140)

def tanimoto_to_leads(fp, lead_fps):
    sims = DataStructs.BulkTanimotoSimilarity(fp, lead_fps)
    return max(sims)

def pdbqt_convert(smi, zinc_id, out_dir):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None: return False
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0: return False
        AllChem.MMFFOptimizeMolecule(mol)
        mk  = meeko.MoleculePreparation()
        mk.prepare(mol)
        pdbqt_str = meeko.PDBQTWriterLegacy.write_string(mk)[0]
        out = out_dir / f"{zinc_id}.pdbqt"
        out.write_text(pdbqt_str)
        return True
    except:
        return False

def main():
    print("="*65)
    print("Script 17 — ZINC20 Similarity-Focused Expansion Screen")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*65)

    # Build lead fingerprints
    print("\nBuilding lead fingerprints...")
    lead_fps = []
    for zid, smi in LEAD_SMILES:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            lead_fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048))
    print(f"  {len(lead_fps)} lead fingerprints ready")

    pains = setup_pains()

    # ── Step 1: Download + similarity filter ───────────────────────────────
    print(f"\nStep 1: Downloading ZINC20 tranches (sim >= {SIM_THRESHOLD})...")
    similar_compounds = []  # [(zinc_id, smiles, max_sim)]
    n_downloaded = 0
    n_tested     = 0
    n_passed_sim = 0

    for tranche in ALL_TRANCHES:
        if len(similar_compounds) >= MAX_COMPOUNDS:
            break
        for sub in ALL_SUBS:
            key = f"{tranche}{sub}"
            if key in DONE:
                continue
            if len(similar_compounds) >= MAX_COMPOUNDS:
                break

            url = f"https://files.docking.org/2D/{tranche}/{key}.smi"
            try:
                r = requests.get(url, headers=HEADERS, timeout=20)
                if r.status_code != 200:
                    continue
                content = r.text.strip()
                if not content:
                    continue

                # Save raw file
                (ZINC20 / f"{key}.smi").write_text(content)
                n_downloaded += 1

                # Parse and filter
                for line in content.splitlines():
                    parts = line.split()
                    if len(parts) < 2:
                        continue
                    smi, zinc_id = parts[0], parts[1]
                    n_tested += 1

                    mol = Chem.MolFromSmiles(smi)
                    if mol is None:
                        continue
                    if not passes_admet(mol):
                        continue
                    if pains.HasMatch(mol):
                        continue

                    fp      = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
                    max_sim = tanimoto_to_leads(fp, lead_fps)

                    if max_sim >= SIM_THRESHOLD:
                        similar_compounds.append((zinc_id, smi, round(max_sim, 3)))
                        n_passed_sim += 1

                if n_downloaded % 10 == 0:
                    print(f"  Tranches: {n_downloaded} | Tested: {n_tested:,} | "
                          f"Sim-hits: {n_passed_sim} ({len(similar_compounds)} total)")

                time.sleep(0.05)

            except Exception as e:
                continue

        # Progress per tranche
        print(f"  Tranche {tranche} done — sim-hits so far: {len(similar_compounds)}")

    print(f"\nDownload complete: {n_downloaded} tranches, {n_tested:,} compounds tested")
    print(f"Similarity-filtered: {len(similar_compounds)} compounds (>= {SIM_THRESHOLD})")

    # ── Step 2: Save filtered SMILES ───────────────────────────────────────
    out_smi = HITS / "enamine_focused_filtered.smi"
    with open(out_smi, 'w') as f:
        for zinc_id, smi, sim in similar_compounds:
            f.write(f"{smi} {zinc_id} sim={sim}\n")
    print(f"Saved: {out_smi}")

    # ── Step 3: PDBQT conversion ───────────────────────────────────────────
    print(f"\nStep 2: Converting {len(similar_compounds)} compounds to PDBQT...")
    n_ok = n_fail = 0
    for zinc_id, smi, sim in similar_compounds:
        if pdbqt_convert(smi, zinc_id, PDBQT):
            n_ok += 1
        else:
            n_fail += 1
        if (n_ok + n_fail) % 500 == 0:
            print(f"  Converted: {n_ok} | Failed: {n_fail}")

    print(f"PDBQT done: {n_ok} OK, {n_fail} failed")
    print(f"Output: {PDBQT}/")

    # ── Summary JSON ───────────────────────────────────────────────────────
    summary = {
        "date": datetime.now().isoformat(),
        "tranches_downloaded": n_downloaded,
        "compounds_tested": n_tested,
        "similarity_threshold": SIM_THRESHOLD,
        "similarity_hits": len(similar_compounds),
        "pdbqt_converted": n_ok,
        "pdbqt_failed": n_fail,
        "output_smi": str(out_smi),
        "output_pdbqt_dir": str(PDBQT),
    }
    with open(HITS / "enamine_focused_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    print("\n" + "="*65)
    print("STATUS: COMPLETE")
    print(f"Next: scp data/pdbqt_enamine/ to NIC5 and submit SLURM array")
    print("="*65)

if __name__ == "__main__":
    main()
