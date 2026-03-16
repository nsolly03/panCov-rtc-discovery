import requests

pdb_ids = ['8D34','8U2X','9MRU','9MRW','9MRY','7DW6','9HH5','7TJ2','9BIH','9M48','9M49']

print('{:<8} {:<12} {:<22} {}'.format('PDB','Resolution','Method','Chains / Entities'))
print('-' * 90)

for pdb_id in pdb_ids:
    try:
        # Basic entry info
        r = requests.get(
            f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}',
            timeout=30).json()

        res    = r.get('rcsb_entry_info',{}).get('resolution_combined',[None])
        res    = res[0] if res else 'N/A'
        method = r.get('rcsb_entry_info',{}).get('experimental_method','?')

        # Polymer entity descriptions via GraphQL
        gql = '''
        {
          entry(entry_id: "%s") {
            polymer_entities {
              rcsb_polymer_entity { pdbx_description }
              rcsb_entity_instance_count
            }
          }
        }
        ''' % pdb_id

        g = requests.get(
            'https://data.rcsb.org/graphql',
            params={'query': gql},
            timeout=30).json()

        entities = (g.get('data') or {}) \
                    .get('entry', {}) \
                    .get('polymer_entities', [])

        descs = []
        for e in entities:
            desc  = (e.get('rcsb_polymer_entity') or {}).get('pdbx_description','?')
            count = e.get('rcsb_entity_instance_count', 1)
            descs.append('{} (x{})'.format(desc, count))

        entity_str = ' | '.join(descs)
        print('{:<8} {:<12} {:<22} {}'.format(
            pdb_id, str(res), str(method), entity_str))

    except Exception as ex:
        print('{:<8} ERROR: {}'.format(pdb_id, ex))

print()
print('Key: look for entries where entity list contains TWO different proteins')
print('     e.g. NSP15 + NSP7 or NSP15 + NSP8 = potential new pipeline target')
