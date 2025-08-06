from Bio import PDB
from classes.TrimmedChain import TrimmedChain

input_file = str(snakemake.input)
output_file = str(snakemake.output)
padding_n = snakemake.params['padding_n']
padding_c = snakemake.params['padding_c']

parser = PDB.PDBParser(QUIET = True)
structure = parser.get_structure("protein", input_file)
dssp = PDB.DSSP(structure[0], input_file)

def find_helix(structure, pdb_file):
    """
    Locate first and last helix positions inside the structure
    """
    helix_codes = [ 'H', 'I', 'G' ]
    first_pos = last_pos = None
    for key in dssp.keys():
        pos, res, sec_struct, rASA, phi, psi, NH__O_1_i, NH__O_1_e, O__NH_1_i, O__NH_1_e, NH__O_2_i, NH__O_2_e, O__NH_2_i, O__NH_2_e = dssp[key]
        if sec_struct in helix_codes:
            if first_pos is None:
                first_pos = pos
            last_pos = pos
    return first_pos, last_pos

def get_seqres(pdb_file):
    """
    Extract the SEQRES records
    """
    seqres = ''
    with open(input_file) as file:
        for line in file:
            if line.startswith('SEQRES '):
                seqres += line
    return seqres

def trim_pdb(structure, file, first_pos, last_pos, padding_n, padding_c):
    """
    Trim a pdb structure to the region [start;stop] 1-based
    and save in the the structure in output pdb file
    """
    io = PDB.PDBIO()
    io.set_structure(structure)
    trimmed_chain = TrimmedChain(model = 0, chain = 'A', start = first_pos - padding_n, end = last_pos + padding_c)
    io.save(file, trimmed_chain)

first_pos, last_pos = find_helix(structure, input_file)
assert first_pos is not None, "Did not find helices in the structure"

seqres = get_seqres(input_file)

with open(output_file, 'w') as file:
    file.write(seqres)
    trim_pdb(structure, file, first_pos, last_pos, padding_n, padding_c)
