import requests

pdb_ids = ['8D34','8U2X','9MRU','9MRW','9MRY','7DW6','9HH5','7TJ2','9BIH','9M48','9M49']

print('{:<8} {:<10} {:<22} {}'.format('PDB','Res(A)','Method','Entities'))
print('-' * 100)

for pdb_id in pdb_ids:
    try:
        # Entry info
        entry = requests.get(
            f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}',
            timeout=30).json()
        res    = entry.get('rcsb_entry_info',{}).get('resolution_combined',[None])
        res    = round(res[0],2) if res and res[0] else 'N/A'
        method = entry.get('rcsb_entry_info',{}).get('experimental_method','?')

        # Get polymer entity IDs from entry
        entity_ids = entry.get('rcsb_entry_container_identifiers',{}) \
                          .get('polymer_entity_ids', [])

        descs = []
        for eid in entity_ids:
            url = f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{eid}'
            e = requests.get(url, timeout=30).json()
            desc  = e.get('rcsb_polymer_entity',{}).get('pdbx_description','?')
            count = e.get('rcsb_entity_instance_count', '?')
            # Get source organism
            src = e.get('rcsb_entity_source_organism',[{}])
            org = src[0].get('ncbi_scientific_name','?') if src else '?'
            descs.append('{} (x{}, {})'.format(desc, count, org))

        print('{:<8} {:<10} {:<22} {}'.format(
            pdb_id, str(res), str(method), ' | '.join(descs)))

    except Exception as ex:
        print('{:<8} ERROR: {}'.format(pdb_id, ex))

print()
print('DECISION RULE:')
print('  1 entity (NSP15 only)     -> homodimer/hexamer, no new interface')
print('  2+ entities, both NSPs    -> NEW heterodimer interface candidate')
print('  2+ entities, NSP15+small  -> inhibitor-bound, useful for Script 16_9 validation')
