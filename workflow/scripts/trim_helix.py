from Bio import PDB

input_file = str(snakemake.input)
output_file = str(snakemake.output)

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

def trim_pdb(structure, output_pdb, first_pos, last_pos):
    """
    Trim a pdb structure to the region [start;stop] 1-based
    and save in the the structure in output pdb file
    """
    io = PDB.PDBIO()
    io.set_structure(structure)
    class TrimmedStructure(PDB.Select):
        def __init__(self, model, chain):
            self.pos = 0
            self.model = model
            self.chain = chain
        def accept_residue(self, residue):
            protein, model, chain, (hetflag, res_id, icode) = residue.full_id
            is_hetatom = hetflag.strip()
            if not is_hetatom and model == self.model and chain == self.chain:
                self.pos += 1
                return first_pos <= self.pos <= last_pos
    io.save(output_pdb, TrimmedStructure(0, 'A'))

first_pos, last_pos = find_helix(structure, input_file)
assert first_pos is not None, "Did not find helices in the structure"

trim_pdb(structure, output_file, first_pos, last_pos)
