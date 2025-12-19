from pymol import cmd, stored
from psico.exporting import save_pdb
from collections import defaultdict
import json
from os import path
from tempfile import NamedTemporaryFile

pdb = snakemake.wildcards['pdb']
atom_map_file = snakemake.input['atom_map']
output_file = str(snakemake.output)

with open(atom_map_file) as file:
    atom_map = json.load(file)

cmd.fetch(pdb)

# Expect one chain
stored.chains = set()
cmd.iterate("(all)", "stored.chains.add(chain)")
assert len(stored.chains) == 1, f"Expected exactly one chain in {pdb}, got {len(stored.chains)}"

# Cleanup
cmd.remove("(hydro)")
cmd.remove("not alt ''+A")
cmd.remove("not resi 1-") # remove residues before the start of the protein
cmd.alter("(all)", "segi, alt, chain = '', '', 'A'")

# Rename the atoms of the ligand (and of the residue connected to it)
for resn_from, components in atom_map.items():
    for resn_to, res_data in components.items():
        res_type = res_data['type']
        for atom_from, atom_to in res_data['atoms'].items():
            stored.resn_name_type = resn_to, atom_to, res_type
            cmd.alter(f"resn {resn_from} and name {atom_from}", "resn, name, type = stored.resn_name_type")
    stored.extra_atoms = []
    cmd.iterate(f"resn {resn_from}", "stored.extra_atoms.append(name)")
    assert not stored.extra_atoms, f"Found additional atoms for ligand {resn_from}: {stored.extra_atoms}"

# Fix resi's in structures with gaps (see 8XX8)
stored.polymer_resis = set()
cmd.iterate('polymer', 'stored.polymer_resis.add(int(resi))')
min_resi = min(stored.polymer_resis)
stored.update_resis = {}
for i, resi in enumerate(sorted(stored.polymer_resis)):
    if resi != i + min_resi:
        stored.update_resis[str(resi)] = i + min_resi
if stored.update_resis:
    cmd.alter('polymer', 'resi = stored.update_resis[resi] if resi in stored.update_resis else resi')

# Remove all non-covalent ligands
cmd.remove("hetatm and not byres bound_to polymer")

# Re-number non-polymer resi's
stored.polymer_residues = set()
stored.not_polymer_residues = set()
cmd.iterate("polymer", "stored.polymer_residues.add(int(resi))")
cmd.iterate("not polymer", "stored.not_polymer_residues.add(int(resi))")

stored.target_resi = max(stored.polymer_residues)

# temporally re-number them to put outside of the target range
stored.offset = stored.target_resi + len(stored.not_polymer_residues) + 1000
cmd.alter("not polymer", "resi = stored.offset + int(resi)")
for resi in stored.not_polymer_residues:
    from_resi = resi + stored.offset
    stored.target_resi += 1
    cmd.alter(f"(not polymer) and resi {from_resi}", "resi = stored.target_resi")

cmd.sort()
cmd.set('pdb_conect_nodup', 0)
save_pdb(output_file, seqres = True)
