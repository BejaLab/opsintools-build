ref_pdb = 'workflow/pdb.txt'
fasta = snakemake.input[0]
clstr = clstr_file = f"{snakemake.input[0]}.clstr"

# Function that takes a clstr file and a name of a protein in string, finds it in the clstr file and returns the
# representative sequence in that cluster
def representative_protein(file_path, query_protein):
    with open(file_path, 'r') as clstr_file:
        lines = clstr_file.readlines()
        i = 0
        for line_number, line in enumerate(lines, start=1):
            if line.startswith(">"):
                i = line_number
            elif query_protein in line:
                for l in lines[i:]:
                    if "*" in l:
                        return l.split('>', 1)[1].split('...', 1)[0].strip()
                    elif l.startswith(">"):
                        break
    return None

# Create a list called must_seq with all the reference sequences we want to add to the clustered files
with open(ref_pdb, 'r') as ref_seq_raw:
    must_seq = []
    for l in ref_seq_raw:
        must_seq.append(l.split()[0])

# Check what files from the list of all reference sequences are actually missing from the current clustered file
# Add the missing files to a list called missing_seq
with open(fasta, 'r') as clustered_fasta:
    clustered_fasta_text = clustered_fasta.read()
    missing_seq = []
    for must in must_seq:
        if must not in clustered_fasta_text:
            missing_seq.append(must)
    # For each missing seq run the function representative_protein and save it in a txt file.
    representative_delete =[]
    for missing in missing_seq:
        rep_to_del = representative_protein(clstr, missing)
        if rep_to_del is not None:
            representative_delete.append(rep_to_del)

# Create a text file that has all the id's of the wrong representatives that need to be removed
with open(snakemake.output[0], 'w') as f:
    for i in representative_delete:
        f.write(f"{i}\n")

# Create a text file with the missing reference sequences
with open(snakemake.output[1], 'w') as f:
    for i in missing_seq:
        f.write(f"{i}\n")
