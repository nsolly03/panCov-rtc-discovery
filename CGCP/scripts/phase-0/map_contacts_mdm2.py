#!/usr/bin/env python3
"""
Contact mapping for MDM2/p53 validation
Known hotspots: PHE19, TRP23, LEU26 on chain B (p53 peptide)
"""

from Bio import PDB
import numpy as np
import csv
import sys


def map_contacts(pdb_file, chain_a, chain_b, cutoff=5.0):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('complex', pdb_file)
    model = structure[0]

    chain_a_residues = list(model[chain_a].get_residues())
    chain_b_residues = list(model[chain_b].get_residues())

    contacts = []
    for res_a in chain_a_residues:
        if res_a.get_id()[0] != ' ':
            continue
        res_a_id = res_a.get_id()[1]
        res_a_name = res_a.get_resname()

        for res_b in chain_b_residues:
            if res_b.get_id()[0] != ' ':
                continue
            res_b_id = res_b.get_id()[1]
            res_b_name = res_b.get_resname()

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


def check_known_hotspots(contacts):
    """
    Known p53 hotspots on chain B: PHE19, TRP23, LEU26
    These insert into MDM2 hydrophobic cleft
    """
    known_hotspots = ['PHE19', 'TRP23', 'LEU26']

    print(f"\nKnown Hotspot Validation (Chain B — p53 peptide):")
    print("-" * 50)

    all_passed = True
    for hotspot in known_hotspots:
        hotspot_contacts = [c for c in contacts if hotspot in c['res_b']]
        print(f"  {hotspot}: {len(hotspot_contacts)} contacts")

        sorted_contacts = sorted(hotspot_contacts, key=lambda x: x['distance'])
        for c in sorted_contacts[:3]:
            print(f"    -> {c['res_a']}: {c['distance']}A")

        if len(hotspot_contacts) == 0:
            print(f"    WARNING: {hotspot} has no contacts!")
            all_passed = False

    return all_passed


def main():
    pdb_file = '1YCR.pdb'
    chain_a, chain_b = 'A', 'B'

    print("=" * 60)
    print("MDM2/p53 Contact Mapping")
    print("=" * 60)
    print(f"PDB: {pdb_file}")
    print(f"Chains: {chain_a} (MDM2) + {chain_b} (p53 peptide)")
    print(f"Distance cutoff: 5.0A")
    print(f"Known hotspots: PHE19, TRP23, LEU26 (chain B)")
    print("")

    print("Mapping contacts...")
    contacts = map_contacts(pdb_file, chain_a, chain_b)

    output_file = 'mdm2_contacts.tsv'
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f,
                               fieldnames=['chain_a', 'res_a', 'chain_b', 'res_b', 'distance'],
                               delimiter='\t')
        writer.writeheader()
        writer.writerows(contacts)

    print(f"Total contacts found: {len(contacts)}")
    print(f"Saved to: {output_file}")

    passed = check_known_hotspots(contacts)

    print(f"\nSummary:")
    print(f"  Unique MDM2 residues at interface: {len(set(c['res_a'] for c in contacts))}")
    print(f"  Unique p53 residues at interface:  {len(set(c['res_b'] for c in contacts))}")

    print("")
    if passed:
        print("VALIDATION PASSED: All three p53 hotspots detected")
    else:
        print("VALIDATION FAILED: Missing hotspot contacts")
        sys.exit(1)


if __name__ == "__main__":
    main()
