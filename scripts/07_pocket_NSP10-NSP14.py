#!/usr/bin/env python3
"""
Script 07: Pocket Detection and Druggability Analysis - NSP10-NSP14
======================================================

Project : panCov-rtc-discovery
Author  : Olivier Nsekuye
Date    : 2025

What this script does:
  Runs fpocket on crystal structures (7DIY, 5C8T) and AF3 model,
  identifies pockets overlapping with conserved hotspot residues (HIS80, ASP126),
  scores druggability, and defines docking box coordinates for virtual screening.

Structures: 7DIY (SARS-CoV-2), 5C8T (SARS-CoV-1), AF3 model

How to run:
  conda activate rtc-discovery
  cd ~/projects/rtc-pan-coronavirus
  python scripts/07_pocket_NSP10-NSP14.py
"""

import json
import subprocess
import re
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np

# Project paths
PROJECT_ROOT = Path.home() / "projects" / "rtc-pan-coronavirus"
REFERENCE_DIR = PROJECT_ROOT / "00-reference" / "pdb_structures"
AF3_DIR = PROJECT_ROOT / "01-alphafold3" / "NSP10-NSP14"
VALIDATION_DIR = PROJECT_ROOT / "02-validation" / "NSP10-NSP14"
OUTPUT_DIR = VALIDATION_DIR

# Conserved hotspot residues from Script 06 (Entry 013)
CONSERVED_HOTSPOTS = {
    "NSP10": [5, 19, 21, 40, 42, 44, 45, 80, 93],
    "NSP14": [4, 7, 8, 9, 10, 20, 25, 27, 127, 201]
}

# Primary salt bridge residues
PRIMARY_TARGET = {
    "NSP10": 80,
    "NSP14": 126
}


def run_fpocket(pdb_path: Path, structure_name: str) -> Path:
    """
    Run fpocket on a PDB file.
    Returns path to output directory.
    """
    pdb_name = pdb_path.stem
    standard_output = pdb_path.parent / f"{pdb_name}_out"
    custom_output = pdb_path.parent / f"{structure_name}_out"
    
    print(f"\n[FPocket] Running fpocket on {pdb_path.name}")
    
    # If custom output exists, use it
    if custom_output.exists():
        print(f"[FPocket] Using existing output: {custom_output}")
        return custom_output
    
    # If standard output exists, rename it to custom name
    if standard_output.exists():
        print(f"[FPocket] Renaming existing output to: {custom_output}")
        standard_output.rename(custom_output)
        return custom_output
    
    # Run fpocket
    cmd = ["fpocket", "-f", str(pdb_path)]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=pdb_path.parent)
        
        if result.returncode != 0:
            print(f"[ERROR] fpocket failed: {result.stderr}")
            return None
        
        # fpocket creates {pdb_name}_out - rename to {structure_name}_out
        if standard_output.exists():
            standard_output.rename(custom_output)
            print(f"[FPocket] ✓ Completed: {custom_output}")
            return custom_output
        else:
            print(f"[ERROR] Expected output not created: {standard_output}")
            return None
        
    except FileNotFoundError:
        print(f"[ERROR] fpocket not found. Run: conda install -c bioconda fpocket")
        return None
    except Exception as e:
        print(f"[ERROR] fpocket execution failed: {e}")
        return None


def parse_pocket_residues(pocket_pdb: Path) -> List[Dict]:
    """Parse residues from a pocket PDB file."""
    residues = []
    seen = set()
    
    with open(pocket_pdb, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    chain_id = line[21].strip()
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    
                    key = (chain_id, res_num)
                    if key not in seen:
                        seen.add(key)
                        residues.append({
                            'chain': chain_id,
                            'residue_number': res_num,
                            'residue_name': res_name,
                            'x': x, 'y': y, 'z': z
                        })
                except (ValueError, IndexError):
                    continue
    
    return residues


def calculate_pocket_center(residues: List[Dict]) -> Tuple[float, float, float]:
    """Calculate geometric center of pocket residues."""
    if not residues:
        return (0.0, 0.0, 0.0)
    
    x_coords = [r['x'] for r in residues]
    y_coords = [r['y'] for r in residues]
    z_coords = [r['z'] for r in residues]
    
    return (float(np.mean(x_coords)), float(np.mean(y_coords)), float(np.mean(z_coords)))


def parse_fpocket_pockets(output_dir: Path, structure_name: str) -> List[Dict]:
    """
    Parse fpocket output to extract pocket information.
    """
    pockets = []
    
    # Look for {structure_name}_info.txt in output directory
    info_file = output_dir / f"{structure_name}_info.txt"
    
    print(f"[Parse] Looking for: {info_file}")
    
    if not info_file.exists():
        # Try to find any _info.txt file
        info_files = list(output_dir.glob("*_info.txt"))
        if info_files:
            info_file = info_files[0]
            print(f"[Parse] Found alternative: {info_file}")
        else:
            print(f"[ERROR] No info file found in {output_dir}")
            # List directory contents for debugging
            print(f"[DEBUG] Directory contents:")
            for item in output_dir.iterdir():
                print(f"  - {item.name}")
            return pockets
    
    with open(info_file, 'r') as f:
        content = f.read()
    
    # Parse pocket entries
    lines = content.splitlines()
    current_pocket = None
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('Pocket'):
            match = re.search(r'Pocket\s+(\d+)', line)
            if match:
                current_pocket = {
                    'id': int(match.group(1)),
                    'pocket_number': int(match.group(1)),
                    'residues': [],
                    'center': None,
                    'score': 0.0,
                    'druggability_score': 0.0,
                    'volume': 0.0
                }
                pockets.append(current_pocket)
        
        elif current_pocket and 'Score' in line and ':' in line and 'Druggability' not in line:
            parts = line.split(':')
            if len(parts) == 2:
                try:
                    current_pocket['score'] = float(parts[1].strip())
                except ValueError:
                    pass
        
        elif current_pocket and 'Druggability Score' in line and ':' in line:
            parts = line.split(':')
            if len(parts) == 2:
                try:
                    current_pocket['druggability_score'] = float(parts[1].strip())
                except ValueError:
                    pass
        
        elif current_pocket and 'Volume' in line and ':' in line:
            parts = line.split(':')
            if len(parts) == 2:
                try:
                    current_pocket['volume'] = float(parts[1].strip())
                except ValueError:
                    pass
    
    print(f"[Parse] Parsed {len(pockets)} pockets from info file")
    
    # Parse pocket PDB files
    pockets_dir = output_dir / "pockets"
    
    if not pockets_dir.exists():
        print(f"[WARNING] Pockets subdirectory not found: {pockets_dir}")
        return pockets
    
    for pocket in pockets:
        pocket_pdb = pockets_dir / f"pocket{pocket['pocket_number']}_atm.pdb"
        
        if pocket_pdb.exists():
            residues = parse_pocket_residues(pocket_pdb)
            pocket['residues'] = residues
            if residues:
                pocket['center'] = calculate_pocket_center(residues)
        else:
            print(f"[WARNING] Missing pocket file: {pocket_pdb}")
    
    return pockets


def check_hotspot_overlap(pocket: Dict, hotspots: Dict[str, List[int]], 
                          chain_map: Dict[str, str]) -> Dict:
    """Check if pocket overlaps with hotspot residues."""
    overlap = {
        'has_overlap': False,
        'overlapping_residues': [],
        'primary_target_included': False,
        'nsp10_hits': [],
        'nsp14_hits': []
    }
    
    pocket_residue_set = set()
    for res in pocket['residues']:
        key = (res['chain'], res['residue_number'])
        pocket_residue_set.add(key)
    
    # Check NSP10 hotspots
    nsp10_chain = chain_map.get('NSP10', 'A')
    for res_num in hotspots.get('NSP10', []):
        if (nsp10_chain, res_num) in pocket_residue_set:
            overlap['has_overlap'] = True
            overlap['nsp10_hits'].append(res_num)
            overlap['overlapping_residues'].append(f"NSP10:{res_num}")
            if res_num == PRIMARY_TARGET['NSP10']:
                overlap['primary_target_included'] = True
    
    # Check NSP14 hotspots
    nsp14_chain = chain_map.get('NSP14', 'B')
    for res_num in hotspots.get('NSP14', []):
        if (nsp14_chain, res_num) in pocket_residue_set:
            overlap['has_overlap'] = True
            overlap['nsp14_hits'].append(res_num)
            overlap['overlapping_residues'].append(f"NSP14:{res_num}")
            if res_num == PRIMARY_TARGET['NSP14']:
                overlap['primary_target_included'] = True
    
    return overlap


def calculate_docking_box(pocket: Dict, padding: float = 10.0) -> Dict:
    """Calculate docking box parameters from pocket residues."""
    if not pocket['residues']:
        return None
    
    x_coords = [r['x'] for r in pocket['residues']]
    y_coords = [r['y'] for r in pocket['residues']]
    z_coords = [r['z'] for r in pocket['residues']]
    
    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)
    z_min, z_max = min(z_coords), max(z_coords)
    
    center_x = (x_min + x_max) / 2
    center_y = (y_min + y_max) / 2
    center_z = (z_min + z_max) / 2
    
    size_x = (x_max - x_min) + 2 * padding
    size_y = (y_max - y_min) + 2 * padding
    size_z = (z_max - z_min) + 2 * padding
    
    return {
        'center_x': round(center_x, 3),
        'center_y': round(center_y, 3),
        'center_z': round(center_z, 3),
        'size_x': round(size_x, 3),
        'size_y': round(size_y, 3),
        'size_z': round(size_z, 3),
        'padding': padding
    }


def select_best_pocket(pockets: List[Dict], hotspots: Dict[str, List[int]], 
                       chain_map: Dict[str, str]) -> Optional[Dict]:
    """Select the best pocket for docking."""
    if not pockets:
        return None
    
    scored_pockets = []
    
    for pocket in pockets:
        overlap = check_hotspot_overlap(pocket, hotspots, chain_map)
        
        score = 0.0
        if overlap['primary_target_included']:
            score += 1000
        score += len(overlap['nsp10_hits']) * 10
        score += len(overlap['nsp14_hits']) * 10
        score += pocket.get('druggability_score', 0) * 5
        score += pocket.get('score', 0)
        
        scored_pockets.append({
            'pocket': pocket,
            'overlap': overlap,
            'score': score,
            'has_primary_target': overlap['primary_target_included']
        })
    
    scored_pockets.sort(key=lambda x: x['score'], reverse=True)
    
    if scored_pockets:
        best = scored_pockets[0]
        print(f"\n[SELECT] Pocket {best['pocket']['pocket_number']} selected (score: {best['score']:.1f})")
        print(f"[SELECT] Primary target: {best['has_primary_target']}")
        print(f"[SELECT] Hotspots: {best['overlap']['overlapping_residues']}")
        return best
    
    return None


def analyze_structure(pdb_path: Path, structure_name: str, 
                      chain_map: Dict[str, str]) -> Dict:
    """Complete pocket analysis pipeline for one structure."""
    print(f"\n{'='*60}")
    print(f"ANALYZING: {structure_name}")
    print(f"FILE: {pdb_path}")
    print(f"{'='*60}")
    
    output_dir = run_fpocket(pdb_path, structure_name)
    
    if not output_dir or not output_dir.exists():
        return {
            'structure': structure_name,
            'status': 'FAILED',
            'error': 'fpocket execution failed'
        }
    
    pockets = parse_fpocket_pockets(output_dir, structure_name)
    
    if not pockets:
        return {
            'structure': structure_name,
            'status': 'NO_POCKETS',
            'error': 'No pockets detected'
        }
    
    print(f"\n[Pockets Summary]")
    for p in pockets:
        print(f"  Pocket {p['pocket_number']}: Score={p.get('score', 0):.3f}, "
              f"Druggability={p.get('druggability_score', 0):.3f}, "
              f"Volume={p.get('volume', 0):.1f} Å³")
    
    best = select_best_pocket(pockets, CONSERVED_HOTSPOTS, chain_map)
    
    if not best:
        return {
            'structure': structure_name,
            'status': 'NO_SUITABLE_POCKET',
            'pockets_found': len(pockets),
            'error': 'No pocket overlaps with hotspots'
        }
    
    docking_box = calculate_docking_box(best['pocket'], padding=10.0)
    
    return {
        'structure': structure_name,
        'pdb_file': str(pdb_path),
        'status': 'SUCCESS',
        'total_pockets': len(pockets),
        'selected_pocket': {
            'pocket_number': best['pocket']['pocket_number'],
            'fpocket_score': float(best['pocket'].get('score', 0)),
            'druggability_score': float(best['pocket'].get('druggability_score', 0)),
            'volume': float(best['pocket'].get('volume', 0)),
            'num_residues': len(best['pocket']['residues']),
            'center': best['pocket']['center']
        },
        'hotspot_overlap': {
            'primary_target_included': best['overlap']['primary_target_included'],
            'nsp10_hits': best['overlap']['nsp10_hits'],
            'nsp14_hits': best['overlap']['nsp14_hits'],
            'all_overlapping': best['overlap']['overlapping_residues']
        },
        'docking_box': docking_box
    }


def main():
    """Main execution function."""
    
    print("="*70)
    print("SCRIPT 07: Pocket Detection and Druggability Analysis")
    print("Target: NSP10-NSP14")
    print("="*70)
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    structures = []
    
    pdb_7diy = REFERENCE_DIR / "7DIY.pdb"
    if pdb_7diy.exists():
        structures.append({
            'path': pdb_7diy,
            'name': '7DIY',
            'chain_map': {'NSP10': 'A', 'NSP14': 'B'}
        })
    
    pdb_5c8t = REFERENCE_DIR / "5C8T.pdb"
    if pdb_5c8t.exists():
        structures.append({
            'path': pdb_5c8t,
            'name': '5C8T',
            'chain_map': {'NSP10': 'A', 'NSP14': 'B'}
        })
    
    af3_pdb = AF3_DIR / "NSP10_NSP14_best_model.pdb"
    if af3_pdb.exists():
        structures.append({
            'path': af3_pdb,
            'name': 'AF3',
            'chain_map': {'NSP10': 'A', 'NSP14': 'B'}
        })
    
    if not structures:
        print("[ERROR] No structures found. Exiting.")
        return 1
    
    print(f"\n[Setup] Will analyze {len(structures)} structure(s)")
    
    all_results = []
    
    for struct in structures:
        result = analyze_structure(struct['path'], struct['name'], struct['chain_map'])
        all_results.append(result)
    
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    
    successful = [r for r in all_results if r['status'] == 'SUCCESS']
    failed = [r for r in all_results if r['status'] != 'SUCCESS']
    
    print(f"\nTotal: {len(all_results)}, Successful: {len(successful)}, Failed: {len(failed)}")
    
    if failed:
        print(f"\nFailed:")
        for r in failed:
            print(f"  - {r['structure']}: {r.get('error', 'Unknown')}")
    
    if successful:
        print(f"\n[Docking Boxes]")
        for r in successful:
            box = r['docking_box']
            print(f"\n{r['structure']}:")
            print(f"  Center: ({box['center_x']}, {box['center_y']}, {box['center_z']})")
            print(f"  Size:   ({box['size_x']}, {box['size_y']}, {box['size_z']})")
            print(f"  Druggability: {r['selected_pocket']['druggability_score']:.3f}")
    
    consensus = None
    for r in successful:
        if r['structure'] == '7DIY':
            consensus = r
            break
    if not consensus and successful:
        consensus = successful[0]
    
    if consensus:
        print(f"\n[Consensus] Based on {consensus['structure']}")
        box = consensus['docking_box']
        print(f"  Center: ({box['center_x']}, {box['center_y']}, {box['center_z']})")
        print(f"  Size:   ({box['size_x']}, {box['size_y']}, {box['size_z']})")
    
    output_file = OUTPUT_DIR / "pocket_analysis.json"
    
    try:
        date_str = subprocess.check_output(['date', '+%Y-%m-%d']).decode().strip()
    except:
        date_str = "2026-03-03"
    
    final_output = {
        'target_complex': 'NSP10-NSP14',
        'date': date_str,
        'primary_anchor_residue': 'HIS80(NSP10)-ASP126(NSP14)',
        'conserved_hotspots': CONSERVED_HOTSPOTS,
        'consensus_docking_box': consensus['docking_box'] if consensus else None,
        'structure_results': all_results
    }
    
    with open(output_file, 'w') as f:
        json.dump(final_output, f, indent=2)
    
    print(f"\n[Output] Saved: {output_file}")
    
    if consensus:
        docking_config = {
            'complex': 'NSP10-NSP14',
            'receptor_file': str(AF3_DIR / "NSP10_NSP14_best_model.pdb"),
            'center_x': consensus['docking_box']['center_x'],
            'center_y': consensus['docking_box']['center_y'],
            'center_z': consensus['docking_box']['center_z'],
            'size_x': consensus['docking_box']['size_x'],
            'size_y': consensus['docking_box']['size_y'],
            'size_z': consensus['docking_box']['size_z'],
            'anchor_residues': {'NSP10': 80, 'NSP14': 126}
        }
        
        docking_file = OUTPUT_DIR / "docking_config.json"
        with open(docking_file, 'w') as f:
            json.dump(docking_config, f, indent=2)
        print(f"[Output] Saved: {docking_file}")
    
    print(f"\n{'='*70}")
    print("Script 07 completed")
    print("Next: Script 08 - Docking preparation")
    print(f"{'='*70}\n")
    
    return 0


if __name__ == "__main__":
    exit(main())
