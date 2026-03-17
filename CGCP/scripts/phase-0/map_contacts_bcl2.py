#!/usr/bin/env python3
"""
Contact mapping for BCL-2/BAX validation
Detects interface residues and validates known hotspots
"""

from Bio import PDB
import numpy as np
import csv
import sys


def map_contacts(pdb_file, chain_a, chain_b, cutoff=5.0):
    """
    Map all residue-residue contacts between two chains
    
    Returns list of contacts with minimum distance
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('complex', pdb_file)
    
    model = structure[0]
    chain_a_residues = list(model[chain_a].get_residues())
    chain_b_residues = list(model[chain_b].get_residues())
    
    contacts = []
    
    for res_a in chain_a_residues:
        # Skip heteroatoms and water
        if res_a.get_id()[0] != ' ':
            continue
            
        res_a_id = res_a.get_id()[1]
        res_a_name = res_a.get_resname()
        
        for res_b in chain_b_residues:
            if res_b.get_id()[0] != ' ':
                continue
                
            res_b_id = res_b.get_id()[1]
            res_b_name = res_b.get_resname()
            
            # Calculate minimum distance between any atoms
            min_dist = float('inf')
            for atom_a in res_a.get_atoms():
                for atom_b in res_b.get_atoms():
                    dist = np.linalg.norm(atom_a.coord - atom_b.coord)
                    if dist < min_dist:
                        min_dist = dist
            
            if min_dist <= cutoff:
                contacts.append({
                    'chain_a': chain_a,
                    'res_a': f"{res_a_name}{res_a_id}",
                    'chain_b': chain_b,
                    'res_b': f"{res_b_name}{res_b_id}",
                    'distance': round(min_dist, 2)
                })
    
    return contacts


def check_known_hotspots(contacts, chain_id='A'):
    """
    Validate that known BCL-2 hotspots are detected
    
    Known hotspots: PHE105, VAL126, PHE146 (BCL-XL)
    """
    known_hotspots = ['PHE105', 'VAL126', 'PHE146']
    
    print(f"\nKnown Hotspot Validation (Chain {chain_id}):")
    print("-" * 50)
    
    for hotspot in known_hotspots:
        hotspot_contacts = [c for c in contacts if hotspot in c['res_a']]
        print(f"  {hotspot}: {len(hotspot_contacts)} contacts")
        
        # Show top 3 closest contacts
        sorted_contacts = sorted(hotspot_contacts, key=lambda x: x['distance'])
        for c in sorted_contacts[:3]:
            print(f"    -> {c['res_b']}: {c['distance']}A")
        
        if len(hotspot_contacts) == 0:
            print(f"    ⚠ WARNING: {hotspot} has no contacts!")


def main():
    pdb_file = '2YXJ.pdb'
    chain_a, chain_b = 'A', 'B'
    
    print("=" * 60)
    print("BCL-2/BAX Contact Mapping")
    print("=" * 60)
    print(f"PDB: {pdb_file}")
    print(f"Chains: {chain_a} (BCL-XL) + {chain_b} (BIM peptide)")
    print(f"Distance cutoff: 5.0A")
    print("")
    
    # Map contacts
    print("Mapping contacts...")
    contacts = map_contacts(pdb_file, chain_a, chain_b)
    
    # Save to file
    output_file = 'bcl2_contacts.tsv'
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, 
                               fieldnames=['chain_a', 'res_a', 'chain_b', 'res_b', 'distance'],
                               delimiter='\t')
        writer.writeheader()
        writer.writerows(contacts)
    
    print(f"Total contacts found: {len(contacts)}")
    print(f"Saved to: {output_file}")
    
    # Validate known hotspots
    check_known_hotspots(contacts, chain_a)
    
    # Summary statistics
    print(f"\nSummary:")
    print(f"  Unique residues in chain {chain_a} at interface: "
          f"{len(set(c['res_a'] for c in contacts))}")
    print(f"  Unique residues in chain {chain_b} at interface: "
          f"{len(set(c['res_b'] for c in contacts))}")
    
    # Check if PHE105 has contacts (critical validation)
    phe105_contacts = [c for c in contacts if 'PHE105' in c['res_a']]
    if len(phe105_contacts) > 0:
        print(f"\n✓ VALIDATION PASSED: PHE105 has {len(phe105_contacts)} contacts")
    else:
        print(f"\n✗ VALIDATION FAILED: PHE105 has no contacts!")
        sys.exit(1)


if __name__ == "__main__":
    main()
