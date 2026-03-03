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
    "NSP10": [5, 19, 21, 40, 42, 44, 45, 80, 93],  # HIS80 is primary target
    "NSP14": [4, 7, 8, 9, 10, 20, 25, 27, 127, 201]  # ASP126 is primary target partner
}

# Primary salt bridge residues (highest priority for docking box)
PRIMARY_TARGET = {
    "NSP10": 80,   # HIS80
    "NSP14": 126   # ASP126
}


def run_fpocket(pdb_path: Path, output_suffix: str = "") -> Path:
    """
    Run fpocket on a PDB file.
    """
    pdb_name = pdb_path.stem
    output_name = f"{pdb_name}{output_suffix}"
    output_dir = pdb_path.parent / f"{output_name}_out"
    
    print(f"\n[FPocket] Running fpocket on {pdb_path.name}")
    print(f"[FPocket] Output will be in: {output_dir}")
    
    # Check if fpocket output already exists
    if output_dir.exists():
        print(f"[FPocket] Output directory exists, skipping fpocket run")
        return output_dir
    
    # Run fpocket
    cmd = ["fpocket", "-f", str(pdb_path)]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=pdb_path.parent
        )
        
        if result.returncode != 0:
            print(f"[ERROR] fpocket failed: {result.stderr}")
            return None
            
        print(f"[FPocket] ✓ Completed successfully")
        return output_dir
        
    except FileNotFoundError:
        print(f"[ERROR] fpocket not found. Is it installed? (conda install -c bioconda fpocket)")
        return None
    except Exception as e:
        print(f"[ERROR] fpocket execution failed: {e}")
        return None


def parse_pocket_residues(pocket_pdb: Path) -> List[Dict]:
    """
    Parse residues from a pocket PDB file.
    """
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


def parse_fpocket_pockets(output_dir: Path) -> List[Dict]:
    """
    Parse fpocket output to extract pocket information.
    """
    pockets = []
    
    # Parse _info.txt file for pocket metadata
    info_file = output_dir / f"{output_dir.stem.replace('_out', '')}_info.txt"
    
    if not info_file.exists():
        print(f"[WARNING] Info file not found: {info_file}")
        return pockets
    
    with open(info_file, 'r') as f:
        content = f.read()
    
    # Parse pocket entries - fpocket format
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
    
    # Parse pocket PDB files for residue information
    for pocket in pockets:
        pocket_pdb = output_dir / f"pocket{pocket['pocket_number']}_atm.pdb"
        
        if pocket_pdb.exists():
            residues = parse_pocket_residues(pocket_pdb)
            pocket['residues'] = residues
            
            # Calculate pocket center from residues
            if residues:
                pocket['center'] = calculate_pocket_center(residues)
    
    return pockets


def check_hotspot_overlap(pocket: Dict, hotspots: Dict[str, List[int]], 
                          chain_map: Dict[str, str]) -> Dict:
    """
    Check if pocket overlaps with hotspot residues.
    """
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
            
            # Check if this is the primary target (HIS80)
            if res_num == PRIMARY_TARGET['NSP10']:
                overlap['primary_target_included'] = True
    
    # Check NSP14 hotspots
    nsp14_chain = chain_map.get('NSP14', 'B')
    for res_num in hotspots.get('NSP14', []):
        if (nsp14_chain, res_num) in pocket_residue_set:
            overlap['has_overlap'] = True
            overlap['nsp14_hits'].append(res_num)
            overlap['overlapping_residues'].append(f"NSP14:{res_num}")
            
            # Check if this is the primary target partner (ASP126)
            if res_num == PRIMARY_TARGET['NSP14']:
                overlap['primary_target_included'] = True
    
    return overlap


def calculate_docking_box(pocket: Dict, padding: float = 10.0) -> Dict:
    """
    Calculate docking box parameters from pocket residues.
    """
    if not pocket['residues']:
        return None
    
    x_coords = [r['x'] for r in pocket['residues']]
    y_coords = [r['y'] for r in pocket['residues']]
    z_coords = [r['z'] for r in pocket['residues']]
    
    # Calculate bounding box
    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)
    z_min, z_max = min(z_coords), max(z_coords)
    
    # Center
    center_x = (x_min + x_max) / 2
    center_y = (y_min + y_max) / 2
    center_z = (z_min + z_max) / 2
    
    # Dimensions with padding
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
    """
    Select the best pocket for docking based on:
    1. Overlap with primary target (HIS80-ASP126)
    2. Druggability score
    3. Overlap with other conserved hotspots
    """
    if not pockets:
        return None
    
    scored_pockets = []
    
    for pocket in pockets:
        overlap = check_hotspot_overlap(pocket, hotspots, chain_map)
        
        # Scoring criteria
        score = 0.0
        
        # Primary: includes HIS80 or ASP126 (must have for drug target)
        if overlap['primary_target_included']:
            score += 1000
        
        # Secondary: includes other conserved hotspots
        score += len(overlap['nsp10_hits']) * 10
        score += len(overlap['nsp14_hits']) * 10
        
        # Tertiary: fpocket druggability score
        score += pocket.get('druggability_score', 0) * 5
        
        # Quaternary: general fpocket score
        score += pocket.get('score', 0)
        
        scored_pockets.append({
            'pocket': pocket,
            'overlap': overlap,
            'score': score,
            'has_primary_target': overlap['primary_target_included']
        })
    
    # Sort by score (descending)
    scored_pockets.sort(key=lambda x: x['score'], reverse=True)
    
    # Return best pocket
    if scored_pockets:
        best = scored_pockets[0]
        print(f"\n[SELECT] Pocket {best['pocket']['pocket_number']} selected (score: {best['score']:.1f})")
        print(f"[SELECT] Primary target included: {best['has_primary_target']}")
        print(f"[SELECT] Overlapping hotspots: {best['overlap']['overlapping_residues']}")
        return best
    
    return None


def analyze_structure(pdb_path: Path, structure_name: str, 
                      chain_map: Dict[str, str]) -> Dict:
    """
    Complete pocket analysis pipeline for one structure.
    """
    print(f"\n{'='*60}")
    print(f"ANALYZING: {structure_name}")
    print(f"FILE: {pdb_path}")
    print(f"{'='*60}")
    
    # Run fpocket
    output_dir = run_fpocket(pdb_path, f"_{structure_name}")
    
    if not output_dir or not output_dir.exists():
        return {
            'structure': structure_name,
            'status': 'FAILED',
            'error': 'fpocket execution failed'
        }
    
    # Parse pockets
    print(f"[Parse] Reading pocket information from {output_dir}")
    pockets = parse_fpocket_pockets(output_dir)
    
    if not pockets:
        return {
            'structure': structure_name,
            'status': 'NO_POCKETS',
            'error': 'No pockets detected'
        }
    
    print(f"[Parse] Found {len(pockets)} pockets")
    
    # Print pocket summary
    print(f"\n[Pockets Summary]")
    for p in pockets:
        print(f"  Pocket {p['pocket_number']}: Score={p.get('score', 0):.3f}, "
              f"Druggability={p.get('druggability_score', 0):.3f}, "
              f"Volume={p.get('volume', 0):.1f} Å³")
    
    # Select best pocket
    best = select_best_pocket(pockets, CONSERVED_HOTSPOTS, chain_map)
    
    if not best:
        return {
            'structure': structure_name,
            'status': 'NO_SUITABLE_POCKET',
            'pockets_found': len(pockets),
            'error': 'No pocket overlaps with hotspots'
        }
    
    # Calculate docking box
    docking_box = calculate_docking_box(best['pocket'], padding=10.0)
    
    # Compile results
    result = {
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
    
    return result


def main():
    """Main execution function."""
    
    print("="*70)
    print("SCRIPT 07: Pocket Detection and Druggability Analysis")
    print("Target: NSP10-NSP14")
    print("="*70)
    
    # Ensure output directory exists
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Define structures to analyze
    structures = []
    
    # 7DIY crystal structure (SARS-CoV-2)
    pdb_7diy = REFERENCE_DIR / "7DIY.pdb"
    if pdb_7diy.exists():
        structures.append({
            'path': pdb_7diy,
            'name': '7DIY',
            'chain_map': {'NSP10': 'A', 'NSP14': 'B'}
        })
    else:
        print(f"[WARNING] 7DIY.pdb not found at {pdb_7diy}")
    
    # 5C8T crystal structure (SARS-CoV-1)
    pdb_5c8t = REFERENCE_DIR / "5C8T.pdb"
    if pdb_5c8t.exists():
        structures.append({
            'path': pdb_5c8t,
            'name': '5C8T',
            'chain_map': {'NSP10': 'A', 'NSP14': 'B'}
        })
    else:
        print(f"[WARNING] 5C8T.pdb not found at {pdb_5c8t}")
    
    # AF3 model
    af3_pdb = AF3_DIR / "NSP10_NSP14_best_model.pdb"
    if af3_pdb.exists():
        structures.append({
            'path': af3_pdb,
            'name': 'AF3',
            'chain_map': {'NSP10': 'A', 'NSP14': 'B'}
        })
    else:
        print(f"[WARNING] AF3 model not found at {af3_pdb}")
    
    if not structures:
        print("[ERROR] No structures found to analyze. Exiting.")
        return 1
    
    print(f"\n[Setup] Will analyze {len(structures)} structure(s)")
    
    # Analyze each structure
    all_results = []
    
    for struct in structures:
        result = analyze_structure(struct['path'], struct['name'], struct['chain_map'])
        all_results.append(result)
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    
    successful = [r for r in all_results if r['status'] == 'SUCCESS']
    failed = [r for r in all_results if r['status'] != 'SUCCESS']
    
    print(f"\nTotal structures analyzed: {len(all_results)}")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")
    
    # Print docking boxes for successful analyses
    if successful:
        print(f"\n[Docking Box Recommendations]")
        for r in successful:
            box = r['docking_box']
            print(f"\n{r['structure']}:")
            print(f"  Center: ({box['center_x']}, {box['center_y']}, {box['center_z']})")
            print(f"  Size:   ({box['size_x']}, {box['size_y']}, {box['size_z']})")
            print(f"  Primary target included: {r['hotspot_overlap']['primary_target_included']}")
            print(f"  Druggability score: {r['selected_pocket']['druggability_score']:.3f}")
    
    # Select consensus docking box (prioritize 7DIY, then AF3)
    consensus = None
    for r in successful:
        if r['structure'] == '7DIY':
            consensus = r
            break
    if not consensus and successful:
        consensus = successful[0]
    
    if consensus:
        print(f"\n[Consensus Docking Box] Based on {consensus['structure']}")
        box = consensus['docking_box']
        print(f"  Center: ({box['center_x']}, {box['center_y']}, {box['center_z']})")
        print(f"  Size:   ({box['size_x']}, {box['size_y']}, {box['size_z']})")
    
    # Save results
    output_file = OUTPUT_DIR / "pocket_analysis.json"
    
    # Get date
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
    
    print(f"\n[Output] Results saved to: {output_file}")
    
    # Also save simplified docking config for next script
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
            'anchor_residues': {
                'NSP10': 80,
                'NSP14': 126
            }
        }
        
        docking_file = OUTPUT_DIR / "docking_config.json"
        with open(docking_file, 'w') as f:
            json.dump(docking_config, f, indent=2)
        print(f"[Output] Docking config saved to: {docking_file}")
    
    print(f"\n{'='*70}")
    print("Script 07 completed successfully")
    print("Next step: Script 08 - Docking preparation")
    print(f"{'='*70}\n")
    
    return 0


if __name__ == "__main__":
    exit(main())
