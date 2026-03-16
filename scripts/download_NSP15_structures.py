import os
import requests

structures = ["6VWW", "6W01", "6WLC", "6WXC", "2H85"]
outdir = "00-reference/pdb_structures/NSP15"
os.makedirs(outdir, exist_ok=True)

for pdb_id in structures:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    r = requests.get(url, timeout=60)
    if r.status_code == 200:
        path = f"{outdir}/{pdb_id}.pdb"
        with open(path, "w") as f:
            f.write(r.text)
        atom_lines = [l for l in r.text.split("\n")
                      if l.startswith("ATOM") or l.startswith("HETATM")]
        print(f"  {pdb_id}: OK — {len(atom_lines)} ATOM/HETATM lines -> {path}")
    else:
        print(f"  {pdb_id}: FAILED status {r.status_code}")

print("\nDownload complete. Files in 00-reference/pdb_structures/NSP15/")
