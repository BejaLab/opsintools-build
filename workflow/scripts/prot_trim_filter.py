from Bio import SeqIO, PDB
import re

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
aln_file = snakemake.input.aln
pdb_file = snakemake.input.pdb
trimmed_fasta_file = snakemake.output.fasta
trimmed_pdb_file = snakemake.output.pdb

ref_id = snakemake.params.ref
query_id = snakemake.wildcards.pdb
br_prefix_n = snakemake.params.br_prefix_n
br_suffix_c = snakemake.params.br_suffix_c
max_missing_n = snakemake.params.max_missing_n
max_missing_c = snakemake.params.max_missing_c
ref_lys_pos = snakemake.params.ref_lys_pos
max_rmsd = snakemake.params.max_rmsd
min_aln_len = snakemake.params.min_aln_len
min_len = snakemake.params.min_len
max_len = snakemake.params.max_len

# Function to check if query needs to be filtered out. Save aligned sequence for both referance and query
def filter_aln_seq(aln_file, max_rmsd, min_aln_len, min_len, ref_lys_pos):
    len_query = rmsd = aln_len = None
    
    # Read the aln file
    with open(aln_file, "r") as out_file:
        lines = out_file.readlines()
        
        # Filter criteria for each line in the aln file
        for line in lines:
            if "Length of Chain_2:" in line: # Save the query length
                len_query = int(re.search(r'Length\sof\sChain_2:\s*([0-9]+)', line).group(1))
            if "RMSD=" in line: # Save the RMSD score
                rmsd = float(re.search(r'.*RMSD=\s*([0-9]+(?:\.[0-9]*)?)', line).group(1))
            if "Aligned length=" in line: # Save the aligned length
                aln_len = int(re.search(r'.*Aligned\s*length=\s*([0-9]+)', line).group(1))

    # If criteria not met return empty    
    if (len_query is None or len_query < min_len or len_query > max_len or
        rmsd is None or rmsd > max_rmsd or 
        aln_len is None or aln_len < min_aln_len):
        return None, "", ""
    
    # Save the sequences unless wrong file type (not enough lines)
    aln_ref_seq = lines[-4].strip() if len(lines) >= 4 else ""
    aln_query_seq = lines[-2].strip() if len(lines) >= 4 else ""
    
    # Create an indexed dictionary for alinged sequences
    aln_dict = {i: (ref, query) for i, (ref, query) in enumerate(zip(aln_ref_seq, aln_query_seq))}

    # Check if for the referance lysine position there is a lysine in the query 
    nongap_count = 0
    for pos, (ref, query) in aln_dict.items():
        if ref != "-":
            nongap_count += 1
            if nongap_count == ref_lys_pos:
                if query != "K" or ref != "K":
                    return None, "", ""
                break
    
    # If lysine position not reached
    if nongap_count < ref_lys_pos:
        return None, "", ""
    
    # Return the aligned dictionary and sequences whether they are empty or not    
    return aln_dict, aln_ref_seq, aln_query_seq

# Parse alignment function, return first and last nongap position for both fer and query
def parse_aln(aln_dict):
    if aln_dict is None: # If no dictionary return empty
        return None, None, None, None
    
    first_ref_pos = next((pos for pos, (ref, _) in aln_dict.items() if ref != "-"), None)
    first_query_pos = next((pos for pos, (_, query) in aln_dict.items() if query != "-"), None)
    last_ref_pos = next((pos for pos, (ref, _) in reversed(aln_dict.items()) if ref != "-"), None)
    last_query_pos = next((pos for pos, (_, query) in reversed(aln_dict.items()) if query != "-"), None)
    return first_ref_pos, first_query_pos, last_ref_pos, last_query_pos

# Get start and stop positions for trimming
def start_stop(aln_dict, first_ref_pos, first_query_pos, last_ref_pos, last_query_pos):
    # If query is no dictionary or missing too much of TM1 or TM7 return empty 
    if aln_dict is None or first_query_pos - first_ref_pos > max_missing_n or last_ref_pos - last_query_pos > max_missing_c: 
        return -1, -1, ""
    
    # Create a full query sequence and a mapping for ungapped positions
    full_query = ""
    full_query_map = {}
    for pos, (_, query) in aln_dict.items():
        if query != "-":
            full_query_map[pos] = len(full_query)
            full_query += query

    # Start position is either br_prefix_n th position BEFORE the first ref aligned position or first aligned query pos
    prefix_seq_n = [pos for pos, (_, query) in aln_dict.items() if pos < first_ref_pos and query != "-"]
    gap_start = prefix_seq_n[-br_prefix_n] if prefix_seq_n and len(prefix_seq_n) > br_prefix_n else first_query_pos
    start = full_query_map[gap_start]

    # Stop position is either br_suffix_c th position AFTER the last ref aligned position or last aligned query pos
    suffix_seq_c = [pos for pos, (_, query) in aln_dict.items() if pos > last_ref_pos and query != "-"]
    gap_stop = suffix_seq_c[br_suffix_c - 1] if suffix_seq_c and len(suffix_seq_c) > br_suffix_c else last_query_pos
    stop = full_query_map[gap_stop]

    # Return the start and end position of the trim
    return start, stop, full_query

# Function to trim the pdb file according to start and stop positions
def trim_pdb(pdb_file, trimmed_pdb, start, stop):
    parser = PDB.PDBParser(QUIET = True)
    structure = parser.get_structure(query_id, pdb_file)
    first_res_id = next(iter(structure[0]["A"])).id[1]
    io = PDB.PDBIO()
    io.set_structure(structure)
    class TrimmedStructure(PDB.Select):
        def accept_model(self, model):
            return model.id == 0
        def accept_chain(self, chain):
            return chain.id == "A"
        def accept_residue(self, residue):
            res_index = residue.id[1] - 1
            return start + first_res_id - 1 <= res_index <= stop + first_res_id - 1
    io.save(trimmed_pdb, TrimmedStructure())

# Run the filter_aln_seq function on the aln file with all filter params
aln_dict, aln_ref_seq, aln_query_seq = filter_aln_seq(aln_file, max_rmsd, min_aln_len, min_len, ref_lys_pos)

# Run parse alignment function on the alignment dictionary
first_ref_pos, first_query_pos, last_ref_pos, last_query_pos = parse_aln(aln_dict)

# Get start and stop positions for the trim
start, stop, full_query = start_stop(aln_dict, first_ref_pos, first_query_pos, last_ref_pos, last_query_pos)

# Trim the PDB file and save it to output pdb
trim_pdb(pdb_file, trimmed_pdb_file, start, stop)

# Trim the full query sequence and save as fasta
with open(trimmed_fasta_file, "w") as out_file:
    if start != -1:
        trimmed_query_seq = full_query[start: stop + 1]
        out_file.write(f">{query_id}\n{trimmed_query_seq}\n")
    else:
        trimmed_query_seq = ""
        out_file.write(f">{query_id}\n{trimmed_query_seq}\n")
