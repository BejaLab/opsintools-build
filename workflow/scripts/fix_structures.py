from pymol import cmd, stored
from collections import defaultdict
import json
from os import path
from tempfile import NamedTemporaryFile

model_file = snakemake.input['model']
atom_map_file = snakemake.input['atom_map']
output_file = str(snakemake.output)

with open(atom_map_file) as file:
    atom_map = json.load(file)

cmd.load(model_file)
cmd.alter("(all)", "segi = ''")
cmd.remove("not chain A")
cmd.remove("(hydro)")

for resn_from, components in atom_map.items():
    for resn_to, res_data in components.items():
        res_type = res_data['type']
        for atom_from, atom_to in res_data['atoms'].items():
            cmd.alter(f"resn {resn_from} and name {atom_from}", f"resn, name, type = '{resn_to}', '{atom_to}', '{res_type}'")
    stored.extra_atoms = []
    cmd.iterate(f"resn {resn_from}", "stored.extra_atoms.append(name)")
    assert len(stored.extra_atoms) == 0, f"Found additional atoms for ligand {resn_from}: {stored.extra_atoms}"

stored.het_residues = defaultdict(set)
cmd.iterate("hetatm", "stored.het_residues[chain].add(resi)")

# re-assign hetatms to polymer chains they are bound to
# and re-number resi
for het_chain, resis in stored.het_residues.items():
    for resi in resis:
        stored.main_chain = None
        cmd.iterate(f"neighbor (hetatm and chain {het_chain} and resi {resi})", "stored.main_chain = chain")
        if stored.main_chain:
            stored.chain_resis = []
            cmd.iterate(f"chain {stored.main_chain}", "stored.chain_resis.append(int(resi))")
            het_resi = 1 + max(stored.chain_resis)
            cmd.alter(f"hetatm and chain {het_chain} and resi {resi}", f"resi, chain = {het_resi}, stored.main_chain")

# Re-open the file to fix the order of the atoms
# since .sort() is not enough

with NamedTemporaryFile(suffix = '.pdb', delete_on_close = False) as temp:
    cmd.save(temp.name)
    temp.close()
    cmd.reinitialize()
    cmd.load(temp.name)

with open(output_file, 'w') as file:
    for line in cmd.get_pdbstr().splitlines():
        line_type = line.split()[0]
        if line_type in [ 'ATOM', 'TER', 'END' ]:
            file.write(line + '\n')
