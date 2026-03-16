import requests

pdb_id = '9HH5'

# Get full entry
entry = requests.get(
    f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}',
    timeout=30).json()

print('='*60)
print(f'Detailed inspection: {pdb_id}')
print('='*60)

# Resolution and method
res    = entry.get('rcsb_entry_info',{}).get('resolution_combined',[None])
method = entry.get('rcsb_entry_info',{}).get('experimental_method','?')
print(f'Resolution : {res[0] if res else "N/A"} A')
print(f'Method     : {method}')

# Deposition date
dep = entry.get('rcsb_accession_info',{}).get('deposit_date','?')
print(f'Deposited  : {dep}')

# Title
title = entry.get('struct',{}).get('title','?')
print(f'Title      : {title}')

# All polymer entities in detail
entity_ids = entry.get('rcsb_entry_container_identifiers',{}) \
                  .get('polymer_entity_ids', [])

print(f'\nPolymer entities ({len(entity_ids)} total):')
for eid in entity_ids:
    url = f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{eid}'
    e = requests.get(url, timeout=30).json()

    desc   = e.get('rcsb_polymer_entity',{}).get('pdbx_description','?')
    count  = e.get('rcsb_entity_instance_count','?')
    length = e.get('entity_poly',{}).get('rcsb_sample_sequence_length','?')
    chains = e.get('rcsb_entity_instance_container_identifiers',
                   [{}])
    chain_ids = [c.get('auth_asym_id','?') for c in
                 e.get('rcsb_polymer_entity_container_identifiers',{})
                  .get('auth_asym_ids',[])] if \
                 e.get('rcsb_polymer_entity_container_identifiers') else ['?']

    ptype  = e.get('entity_poly',{}).get('rcsb_entity_polymer_type','?')
    src    = e.get('rcsb_entity_source_organism',[{}])
    org    = src[0].get('ncbi_scientific_name','?') if src else '?'

    print(f'  Entity {eid}:')
    print(f'    Description  : {desc}')
    print(f'    Type         : {ptype}')
    print(f'    Length       : {length} aa')
    print(f'    Chains       : {chain_ids}')
    print(f'    Copies       : {count}')
    print(f'    Organism     : {org}')

# Check for ligands
nonpoly = entry.get('rcsb_entry_container_identifiers',{}) \
               .get('non_polymer_entity_ids', [])
if nonpoly:
    print(f'\nNon-polymer entities (ligands): {nonpoly}')
    for lid in nonpoly[:5]:
        url = f'https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{lid}'
        l = requests.get(url, timeout=30).json()
        lname = l.get('pdbx_entity_nonpoly',{}).get('name','?')
        lcomp = l.get('pdbx_entity_nonpoly',{}).get('comp_id','?')
        print(f'  {lid}: {lname} ({lcomp})')

print()
print('INTERPRETATION:')
print('  If entity 1 = NSP15 AND entity 2 = different NSP protein')
print('  -> NEW heterodimer interface candidate for pipeline')
print('  If entity 2 = fragment/domain of polyprotein including NSP15')
print('  -> same protein, not a new interface')
