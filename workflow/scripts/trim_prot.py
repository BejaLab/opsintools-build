from Bio import SeqIO, PDB

# Assign files and ref_id from the snakefile
alignment_file = snakemake.input.aln
pdb_file = snakemake.input.pdb
trimmed_fasta_file = snakemake.output.fasta
trimmed_pdb_file = snakemake.output.pdb

ref_id = snakemake.params.ref
ref_lysine_pos = snakemake.params.ref_lysine_pos
max_missing_n = snakemake.params.max_missing_n
max_missing_c = snakemake.params.max_missing_c
br_prefix_n = snakemake.params.br_prefix_n
br_prefix_c = snakemake.params.br_prefix_c

amino_acids = { 'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',  'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M' }

def get_struct_seq(structure, target_model, target_chain):
    """
    Get the sequence of the PDB structure for a particular chain/model
    and the id (1-based) of the first residue
    """
    seq = ''
    # Build the sequence residue-by-residue
    for res in structure.get_residues():
        protein, model, chain, (hetflag, res_id, icode) = res.full_id
        if chain == target_chain and model == target_model and not hetflag.strip():
            assert res.resname in amino_acids, f"Non-standard amino acid: {res.resname}"
            # Add 'X' for gaps in the structure
            seq += 'X' * (res_id - len(seq) - 1)
            seq += amino_acids[res.resname]
    assert seq, f"Chain {target_chain} does not exist or its sequence is empty"
    return seq

def trim_pdb(structure, output_pdb, start, stop):
    """
    Trim a pdb structure to the region [start;stop) 0-based
    and save in the the structure in output pdb file
    """
    io = PDB.PDBIO()
    io.set_structure(structure)
    class TrimmedStructure(PDB.Select):
        def accept_residue(self, residue):
            res_index = residue.id[1] - 1
            return start <= res_index < stop
    io.save(output_pdb, TrimmedStructure())

def parse_alignment(aln_ref_seq, aln_query_seq, full_query_seq, ref_lysine_pos):
    """
    Parse the alignment keeping track of the positions of the aligned region
    and of the reference non-gap positions
    """
    aln_start = aln_stop = aln_ref_start = -1
    ref_pos = 0
    query_has_lysine = False
    # Locate the aligned region
    for aln_pos, (ref_res, query_res) in enumerate(zip(aln_ref_seq, aln_query_seq)):
        ref_gap = ref_res == '-'
        query_gap = query_res == '-'
        ref_pos += not ref_gap
        if not ref_gap and not query_gap:
            if aln_start < 0:
                aln_start = aln_pos
            aln_stop = aln_pos + 1
            if ref_pos == ref_lysine_pos:
                assert ref_res == 'K', f"Expected lysine at ref position {ref_pos}, got {ref_res}"
                query_has_lysine = query_res == 'K'
        if not ref_gap and aln_ref_start < 0:
            aln_ref_start = aln_pos
    query_start = query_stop = -1
    aln_query_pos = 0
    # Match the alignment of the query with the full-length sequence
    # taking into account eventual gaps
    for full_pos, full_res in enumerate(full_query_seq):
        if full_res != 'X':
            while aln_query_seq[aln_query_pos] == '-':
                aln_query_pos += 1
            assert full_res == aln_query_seq[aln_query_pos], f"The aligned sequence does not match the sequence in PDB at {full_pos}: {full_res} != {aln_query_seq[aln_query_pos]}"
            if aln_query_pos == aln_start:
                query_start = full_pos
            if aln_query_pos == aln_stop - 1:
                query_stop = full_pos + 1
            aln_query_pos += 1
    assert aln_start >= 0, "No aligned region"
    return aln_start, aln_stop, aln_ref_start, query_start, query_stop, query_has_lysine

# Find the reference and the query sequences in the alignment
aln_ref_seq = aln_query_seq = None
query_name = None
for record in SeqIO.parse(alignment_file, "fasta"):
    if record.id == ref_id and aln_ref_seq is None:
        aln_ref_seq = record.seq
    else:
        query_name = record.id
        aln_query_seq = record.seq

assert aln_ref_seq, "Reference sequence not found in the alignment"
assert aln_query_seq and query_name, "Query sequence not found in the alignment"

parser = PDB.PDBParser(QUIET = True)
structure = parser.get_structure("protein", pdb_file)
full_query_seq = get_struct_seq(structure, 0, 'A')

aln_start, aln_stop, aln_ref_start, query_start, query_stop, query_has_lysine = parse_alignment(aln_ref_seq, aln_query_seq, full_query_seq, ref_lysine_pos)

# Obtain the lengths of the reference sequence outside of the
# aligned region
missing_ref_n = len(aln_ref_seq[aln_ref_start:aln_start].ungap())
missing_ref_c = len(aln_ref_seq[aln_stop:].ungap())

# If missing too much, discard the query
if not query_has_lysine or missing_ref_n >= max_missing_n or missing_ref_c >= max_missing_c:
    trimmed_seq = ""
    trim_pdb(structure, trimmed_pdb_file, -1, -1)
else:
    # Add paddinging based on the missing parts of the reference
    # and the additional default padding 
    padding_n = br_prefix_n + missing_ref_n
    padding_c = br_prefix_c + missing_ref_c
    # The part of the query aligned to the reference
    full_query_start = max(query_start - padding_n, 0)
    full_query_stop = min(query_stop + padding_c, len(full_query_seq))
    trimmed_seq = full_query_seq[full_query_start:full_query_stop].strip('X')
    trim_pdb(structure, trimmed_pdb_file, full_query_start, full_query_stop)

# Save the file as a fasta file according to snakemake logic
with open(trimmed_fasta_file, "w") as out_f:
    out_f.write(f">{query_name}\n{trimmed_seq}\n")
