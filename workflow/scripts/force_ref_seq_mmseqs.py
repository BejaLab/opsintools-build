ref_pdb = 'workflow/pdb.txt'
fasta = snakemake.input[0]
tsv = snakemake.input[1]

# Read must have referance sequence file and save in must_seq
with open(ref_pdb, 'r') as ref_seq_raw:
    must_seq = {line.strip().split()[0] for line in ref_seq_raw}

# Read the fasta file of representatives created by mmseqs save the ids
with open(fasta, 'r') as fasta_file:
    clustered_reps_mmseqs = {line.strip().split()[0][1:] for line in fasta_file if line.startswith(">")}
 
# Save the missing IDs from the must have referance
missing_seq = must_seq - clustered_reps_mmseqs

# Create a dictionary mapping the current clustering done by mmseqs
cluster_reps = {}
with open(tsv, 'r') as tsv_file:
    for line in tsv_file:
        rep, member = line.strip().split("\t")
        cluster_reps[member] = rep

# Save the representatives that need to be deleted to be replaced by referance IDs
rep_to_del = {cluster_reps[missing] for missing in missing_seq
              if missing in cluster_reps and cluster_reps[missing] not in must_seq}

# Create a text file that has all the IDs of the wrong representatives that need to be removed
with open(snakemake.output[0], 'w') as out_file:
    for i in rep_to_del:
        out_file.write(f"{i}\n")

# Create a text file with the missing reference sequences
with open(snakemake.output[1], 'w') as out_file:
    for i in missing_seq:
        out_file.write(f"{i}\n")
