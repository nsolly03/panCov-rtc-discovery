#!/usr/bin/env python3
"""
Negative control: IL-2 / IL-2Ralpha contact mapping
Expected result: large flat interface, no dominant anchor,
no compact cluster — confirms pipeline correctly identifies
undruggable PPIs
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


def analyze_interface(contacts, chain_a, chain_b):
    """
    Negative control analysis — we are NOT looking for hotspots.
    We are characterizing the interface as flat and featureless.
    """

    # Count contacts per residue
    res_a_counts = {}
    for c in contacts:
        res = c['res_a']
        res_a_counts[res] = res_a_counts.get(res, 0) + 1

    # Sort by contact count
    sorted_residues = sorted(res_a_counts.items(),
                             key=lambda x: x[1], reverse=True)

    print(f"\nTop 10 chain A residues by contact count:")
    print("-" * 50)
    for res, count in sorted_residues[:10]:
        print(f"  {res}: {count} contacts")

    # Key negative control metrics
    total_res_a = len(set(c['res_a'] for c in contacts))
    total_res_b = len(set(c['res_b'] for c in contacts))
    top_residue_contacts = sorted_residues[0][1] if sorted_residues else 0
    total_contacts = len(contacts)

    # Dominance ratio: top residue contacts / total contacts
    # Low ratio = flat interface (no dominant anchor)
    # High ratio = hotspot-driven interface (druggable)
    dominance_ratio = top_residue_contacts / total_contacts if total_contacts > 0 else 0

    print(f"\nInterface Profile:")
    print("-" * 50)
    print(f"  Total contacts:              {total_contacts}")
    print(f"  Unique chain A residues:     {total_res_a}")
    print(f"  Unique chain B residues:     {total_res_b}")
    print(f"  Top residue contacts:        {top_residue_contacts}")
    print(f"  Dominance ratio:             {dominance_ratio:.3f}")
    print(f"  Avg contacts per residue:    {total_contacts/total_res_a:.1f}")

    # Negative control verdict
    print(f"\nNegative Control Assessment:")
    print("-" * 50)

    flags = []

    if total_contacts > 50:
        print(f"  ✓ Large interface ({total_contacts} contacts) — consistent with undruggable PPI")
    else:
        flags.append(f"WARNING: Fewer contacts than expected ({total_contacts})")

    if total_res_a > 20:
        print(f"  ✓ Many interface residues ({total_res_a}) — extended flat surface")
    else:
        flags.append(f"WARNING: Few interface residues ({total_res_a})")

    if dominance_ratio < 0.15:
        print(f"  ✓ Low dominance ratio ({dominance_ratio:.3f}) — no single anchor residue")
    else:
        flags.append(f"WARNING: High dominance ratio ({dominance_ratio:.3f}) — unexpected hotspot")

    if not flags:
        print(f"\n  NEGATIVE CONTROL PASSED: Interface is large, flat, featureless")
        print(f"  Pipeline correctly characterizes undruggable PPI")
    else:
        print(f"\n  NEGATIVE CONTROL WARNINGS:")
        for f in flags:
            print(f"  - {f}")

    return dominance_ratio, total_contacts, total_res_a


def main():
    pdb_file = '1Z92.pdb'
    chain_a, chain_b = 'A', 'B'

    print("=" * 60)
    print("Negative Control: IL-2 / IL-2Ralpha Contact Mapping")
    print("=" * 60)
    print(f"PDB: {pdb_file}")
    print(f"Chains: {chain_a} (IL-2) + {chain_b} (IL-2Ralpha)")
    print(f"Distance cutoff: 5.0A")
    print(f"Expected: large flat interface, no dominant anchor")
    print("")

    print("Mapping contacts...")
    contacts = map_contacts(pdb_file, chain_a, chain_b)

    output_file = 'il2_contacts.tsv'
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f,
                               fieldnames=['chain_a', 'res_a', 'chain_b',
                                           'res_b', 'distance'],
                               delimiter='\t')
        writer.writeheader()
        writer.writerows(contacts)

    print(f"Total contacts found: {len(contacts)}")
    print(f"Saved to: {output_file}")

    analyze_interface(contacts, chain_a, chain_b)


if __name__ == "__main__":
    main()
