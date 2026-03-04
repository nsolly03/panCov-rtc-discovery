"""
Patch all Script 10_x: add AA identity to residue labels
Regenerate Figs 4-6 and ranking CSVs for all complexes
"""
import re
from pathlib import Path

PROJECT = Path.home() / "projects" / "rtc-pan-coronavirus"

# Map each script to its PDB file and chains
TARGETS = {
    "scripts/10_BSA_alascan_ranking_NSP10-NSP16_2.py": {
        "pdb":    "00-reference/pdb_structures/6W4H.pdb",
        "ch12":   "B",  # NSP10
        "ch16":   "A",  # NSP16
        "labels": {"B":"NSP10","A":"NSP16"},
    },
    "scripts/10_BSA_alascan_ranking_NSP10-NSP14_2.py": {
        "pdb":    "00-reference/pdb_structures/7DIY.pdb",
        "ch10":   "B",
        "ch14":   "A",
        "labels": {"B":"NSP10","A":"NSP14"},
    },
    "scripts/10_BSA_alascan_ranking_NSP12-NSP7_3.py": {
        "pdb":    "00-reference/pdb_structures/7BV2_NSP12-NSP7.pdb",
        "ch12":   "A",
        "ch7":    "C",
        "labels": {"A":"NSP12","C":"NSP7"},
    },
    "scripts/10_BSA_alascan_ranking_NSP12-NSP8_4.py": {
        "pdb":    "00-reference/pdb_structures/7BV2_NSP12-NSP8.pdb",
        "ch12":   "A",
        "ch8":    "B",
        "labels": {"A":"NSP12","B":"NSP8"},
    },
}

# Check which scripts exist
for script_rel in TARGETS:
    p = PROJECT / script_rel
    if p.exists():
        print(f"  Found: {script_rel}")
    else:
        print(f"  Missing: {script_rel}")
