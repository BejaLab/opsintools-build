from Bio import SeqIO

pdb_list_file = snakemake.input['pdb_list']
fasta_file = snakemake.input['fasta']
tsv_file = snakemake.input['tsv']

del_file = snakemake.output['rep_del']
add_file = snakemake.output['rep_add']

# Read must have referance sequence file and save in must_seq
must_seq = set()
with open(pdb_list_file, 'r') as file:
    for line in file:
        must_seq.add(line.split()[0])

# Read the fasta file of representatives created by mmseqs save the ids
clustered_reps_mmseqs = set()
for record in SeqIO.parse(fasta_file, 'fasta'):
    clustered_reps_mmseqs.add(record.id)
 
# Save the missing IDs from the must have referance
missing_seq = must_seq - clustered_reps_mmseqs

# Create a dictionary mapping the current clustering done by mmseqs
cluster_reps = {}
with open(tsv_file, 'r') as file:
    for line in file:
        rep, member = line.strip().split("\t")
        cluster_reps[member] = rep

# Save the representatives that need to be deleted to be replaced by referance IDs
rep_to_del = { cluster_reps[missing] for missing in missing_seq
              if missing in cluster_reps and cluster_reps[missing] not in must_seq }

# Create a text file that has all the IDs of the wrong representatives that need to be removed
with open(del_file, 'w') as file:
    for i in rep_to_del:
        file.write(f"{i}\n")

# Create a text file with the missing reference sequences
with open(add_file, 'w') as file:
    for i in missing_seq:
        file.write(f"{i}\n")
