#!/usr/bin/env python3
"""
Script 08: Docking Preparation for NSP10-NSP14
==============================================

Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

What this script does:
  Prepares the AF3 receptor file for virtual screening by:
  - Loading the consensus docking box from Script 07
  - Preparing receptor (cleaned PDB)
  - Generating AutoDock Vina / VirtualFlow configuration files
  - Creating a summary report for the screening campaign

Input: 02-validation/NSP10-NSP14/docking_config.json
Output: 03-virtual-screening/NSP10-NSP14/ directory with receptor and configs

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/08_docking_prep_NSP10-NSP14.py
"""

import json
import subprocess
import shutil
from pathlib import Path
from Bio import PDB
from Bio.PDB import PDBIO, Select

# Project paths
PROJECT_ROOT = Path.home() / "projects" / "rtc-pan-coronavirus"
VALIDATION_DIR = PROJECT_ROOT / "02-validation" / "NSP10-NSP14"
VS_DIR = PROJECT_ROOT / "03-virtual-screening" / "NSP10-NSP14"
AF3_DIR = PROJECT_ROOT / "01-alphafold3" / "NSP10-NSP14"

# Input/output files
DOCKING_CONFIG = VALIDATION_DIR / "docking_config.json"
AF3_PDB = AF3_DIR / "NSP10_NSP14_best_model.pdb"
RECEPTOR_PDB = VS_DIR / "NSP10_NSP14_receptor.pdb"
RECEPTOR_PDBQT = VS_DIR / "NSP10_NSP14_receptor.pdbqt"
VINA_CONFIG = VS_DIR / "vina_config.txt"
VF_CONFIG = VS_DIR / "virtualflow_config.json"
SUMMARY_FILE = VS_DIR / "docking_summary.json"


class ReceptorSelect(Select):
    """Select only protein residues (no water, no heteroatoms)."""
    
    def accept_residue(self, residue):
        return residue.id[0] == " "


def load_docking_config():
    """Load docking configuration from Script 07."""
    if not DOCKING_CONFIG.exists():
        print(f"[ERROR] Docking config not found: {DOCKING_CONFIG}")
        print(f"[ERROR] Run Script 07 first!")
        return None
    
    with open(DOCKING_CONFIG, 'r') as f:
        config = json.load(f)
    
    print(f"[Load] Config loaded: {config['complex']}")
    print(f"[Load] Box center: ({config['center_x']}, {config['center_y']}, {config['center_z']})")
    print(f"[Load] Box size: ({config['size_x']}, {config['size_y']}, {config['size_z']})")
    
    return config


def prepare_receptor(input_pdb, output_pdb):
    """Clean receptor structure (remove waters)."""
    print(f"\n[Receptor] Preparing from {input_pdb}")
    
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("receptor", input_pdb)
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_pdb), select=ReceptorSelect())
        
        atom_count = sum(1 for _ in structure.get_atoms())
        residue_count = sum(1 for _ in structure.get_residues())
        print(f"[Receptor] Atoms: {atom_count}, Residues: {residue_count}")
        print(f"[Receptor] Saved: {output_pdb}")
        
        return True
        
    except Exception as e:
        print(f"[ERROR] Failed: {e}")
        return False


def generate_vina_config(config, output_file):
    """Generate AutoDock Vina configuration."""
    print(f"\n[Vina] Generating config: {output_file}")
    
    vina_text = f"""# AutoDock Vina Configuration
# Target: {config['complex']}
receptor = {RECEPTOR_PDB.name}

# Search space
center_x = {config['center_x']}
center_y = {config['center_y']}
center_z = {config['center_z']}
size_x = {config['size_x']}
size_y = {config['size_y']}
size_z = {config['size_z']}

# Search parameters
exhaustiveness = 32
num_modes = 9
energy_range = 3
"""
    
    with open(output_file, 'w') as f:
        f.write(vina_text)
    
    print(f"[Vina] Config saved")
    return True


def generate_vf_config(config, output_file):
    """Generate VirtualFlow configuration."""
    print(f"\n[VirtualFlow] Generating config: {output_file}")
    
    vf_config = {
        "target": {
            "name": config['complex'],
            "receptor_file": str(RECEPTOR_PDB.name),
            "receptor_pdbqt": str(RECEPTOR_PDBQT.name)
        },
        "docking_box": {
            "center_x": config['center_x'],
            "center_y": config['center_y'],
            "center_z": config['center_z'],
            "size_x": config['size_x'],
            "size_y": config['size_y'],
            "size_z": config['size_z'],
            "source": "fpocket_consensus_7DIY"
        },
        "anchor_residues": config.get('anchor_residues', {"NSP10": 80, "NSP14": 126}),
        "docking_parameters": {
            "exhaustiveness": 32,
            "num_modes": 9,
            "energy_range": 3.0
        },
        "library": {
            "filter": {
                "mw_min": 200,
                "mw_max": 500,
                "logp_min": -2.0,
                "logp_max": 5.0
            }
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(vf_config, f, indent=2)
    
    print(f"[VirtualFlow] Config saved")
    return True


def generate_summary(config, output_file):
    """Generate summary report."""
    print(f"\n[Summary] Generating: {output_file}")
    
    summary = {
        "campaign": {
            "target": config['complex'],
            "date": subprocess.check_output(['date', '+%Y-%m-%d']).decode().strip(),
            "status": "ready_for_screening"
        },
        "receptor": {
            "source": "AlphaFold3",
            "file": str(RECEPTOR_PDB.name),
            "validation": {"iptm": 0.89, "ptm": 0.88, "f1_score": 0.952}
        },
        "binding_site": {
            "definition": "fpocket_consensus_7DIY",
            "anchor_residues": config.get('anchor_residues', {}),
            "conservation": {
                "NSP10_HIS80": 1.0,
                "note": "Fully conserved across 5 coronaviruses"
            }
        },
        "docking_box": {
            "center": {"x": config['center_x'], "y": config['center_y'], "z": config['center_z']},
            "size": {"x": config['size_x'], "y": config['size_y'], "z": config['size_z']}
        },
        "next_steps": [
            "Prepare receptor PDBQT with MGLTools/OpenBabel",
            "Prepare ligand library",
            "Run virtual screening",
            "Analyze top hits"
        ]
    }
    
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"[Summary] Saved")
    return True


def main():
    """Main execution."""
    
    print("="*70)
    print("SCRIPT 08: Docking Preparation")
    print("Target: NSP10-NSP14")
    print("="*70)
    
    # Create output directory
    VS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"\n[Setup] Output: {VS_DIR}")
    
    # Load config
    config = load_docking_config()
    if not config:
        return 1
    
    # Check AF3 model
    if not AF3_PDB.exists():
        print(f"[ERROR] AF3 model not found: {AF3_PDB}")
        return 1
    
    # Prepare receptor
    if not prepare_receptor(AF3_PDB, RECEPTOR_PDB):
        return 1
    
    # Create PDBQT placeholder (user must prepare properly)
    shutil.copy(RECEPTOR_PDB, RECEPTOR_PDBQT)
    print(f"\n[PDBQT] Placeholder: {RECEPTOR_PDBQT}")
    print(f"[PDBQT] ⚠️  Replace with properly prepared PDBQT before docking!")
    print(f"[PDBQT] Command: prepare_receptor4.py -r {RECEPTOR_PDB} -o {RECEPTOR_PDBQT}")
    
    # Generate configs
    generate_vina_config(config, VINA_CONFIG)
    generate_vf_config(config, VF_CONFIG)
    generate_summary(config, SUMMARY_FILE)
    
    # Final output
    print(f"\n{'='*70}")
    print("SCRIPT 08 COMPLETED")
    print(f"{'='*70}")
    
    print(f"\n[Output Files]")
    print(f"  Receptor (PDB):   {RECEPTOR_PDB}")
    print(f"  Receptor (PDBQT): {RECEPTOR_PDBQT} ⚠️ needs preparation")
    print(f"  Vina config:      {VINA_CONFIG}")
    print(f"  VirtualFlow:      {VF_CONFIG}")
    print(f"  Summary:          {SUMMARY_FILE}")
    
    print(f"\n[Docking Box]")
    print(f"  Center: ({config['center_x']}, {config['center_y']}, {config['center_z']})")
    print(f"  Size:   ({config['size_x']}, {config['size_y']}, {config['size_z']})")
    print(f"  Anchor: HIS80(NSP10)-ASP126(NSP14)")
    
    print(f"\n{'='*70}")
    print("Ready for virtual screening!")
    print(f"{'='*70}\n")
    
    return 0


if __name__ == "__main__":
    exit(main())
