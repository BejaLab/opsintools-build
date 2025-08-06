from Bio import SeqIO, PDB
from Bio.Seq import Seq
from classes.USalign import USalign
from classes.TrimmedChain import TrimmedChain
import json
from math import inf

'''
This script is doing the following: 
 - Filter out proteins with length less than set in config file
 - Filter out proteins with RMSD score higher than set in config file
 - Filter out proteins with alignment length higher than set in config file
 - Filter out proteins with an amino acid that is not Lysine in the 216th position
 - Trim the pdb file of a query protein before and after the TM1 and TM7 according to config file
 - Create a fasta file out of the aln file
'''

# Assign files and params from the snakefile
aln_file = snakemake.input['aln']
query_pdb_file = snakemake.input['query_pdb']
ref_pdb_file = snakemake.input['ref_pdb']
trimmed_fasta_file = snakemake.output['fasta']
trimmed_pdb_file = snakemake.output['pdb']

query_id = snakemake.wildcards['pdb']
pad_n = snakemake.params.get('pad_n', 0)
pad_c = snakemake.params.get('pad_c', 0)

max_missing_n = snakemake.params.get('max_missing_n', inf)
max_missing_c = snakemake.params.get('max_missing_c', inf)

max_rmsd = snakemake.params.get('max_rmsd', inf)
min_aln_len = snakemake.params.get('min_aln_len', 0)
max_aln_len = snakemake.params.get('max_aln_len', inf)

min_len = snakemake.params.get('min_len', 0)
ref_lys_pos = snakemake.params.get('ref_lys_pos', None)

# Function to check if query needs to be filtered out. Save aligned sequence for both referance and query
def filter_aln_seq(aln_file, max_rmsd, min_len, min_aln_len, max_aln_len, ref_lys_pos):
    aln = USalign(aln_file)
    
    assert aln.alignment, f"Alignment is empty in {aln_file}"

    # If criteria not met return empty
    if aln.rmsd > max_rmsd or len(aln.alignment) < min_aln_len or aln.seq_lens[1] < min_len:
        return None

    # Check if for the referance lysine position there is a lysine in the query 
    lys_reached = False
    for pos_r, pos_q, res_r, res_q, distance in aln.alignment:
        if pos_r == ref_lys_pos:
            assert res_r == 'K', f"Lysine not found at reference position {ref_lys_pos}"
            if res_q != "K":
                return None
            lys_reached = True

    if ref_lys_pos is not None and not lys_reached:
        return None

    return aln

# Function to trim the pdb file according to start and stop positions
def trim_struct(structure, seq_id, query_pdb, trimmed_pdb, trimmed_fasta, start, end):
    io = PDB.PDBIO()
    io.set_structure(structure)

    seqres = ''
    with open(query_pdb_file) as file:
        for line in file:
            if line.startswith('SEQRES '):
                seqres += line

    with open(trimmed_pdb, 'w') as file:
        file.write(seqres)
        trimmed_chain = TrimmedChain(model = 0, chain = 'A', start = start, end = end)
        io.save(file, trimmed_chain)

    trimmed_record = SeqIO.SeqRecord(Seq(trimmed_chain.get_seq()), id = seq_id, description = "")
    SeqIO.write(trimmed_record, trimmed_fasta, "fasta")

def get_first_and_last(structure):
    first, *_, last = structure[0]['A'].get_residues()
    return first.id[1], last.id[1]

# Run the filter_aln_seq function on the aln file with all filter params
aln = filter_aln_seq(aln_file, max_rmsd, min_len, min_aln_len, max_aln_len, ref_lys_pos)
success = False
if aln is not None:
    parser = PDB.PDBParser(QUIET = True)
    ref_structure = parser.get_structure('ref', ref_pdb_file)
    query_structure = parser.get_structure('query', query_pdb_file)

    ref_first, ref_last = get_first_and_last(ref_structure)

    aln_first, *_, aln_last = aln.alignment
    ref_first_aln_pos, query_first_aln_pos, *_ = aln_first
    ref_last_aln_pos,  query_last_aln_pos, *_  = aln_last

    missing_n = ref_first - ref_first_aln_pos
    missing_c = ref_last  - ref_last_aln_pos

    if missing_n <= max_missing_n and missing_c <= max_missing_c:
        # Trim the structure and save it to output pdb and its trimmed sequence to fasta
        trim_struct(query_structure, query_id, query_pdb_file, trimmed_pdb_file, trimmed_fasta_file, start = query_first_aln_pos - missing_n - pad_n, end = query_last_aln_pos + missing_c + pad_c)
        success = True

if not success:
    with open(trimmed_pdb_file, 'w'), open(trimmed_fasta_file, 'w'):
        pass
