from Bio import AlignIO

# Assign files and ref_id from the snakefile
alignment = snakemake.input[0]
trimmed = snakemake.output[0]
ref_id = snakemake.params.ref
br_prefix_n = snakemake.params.br_prefix_n
br_prefix_c = snakemake.params.br_prefix_c


# Create the alignment file in BioPython
aln = AlignIO.read(alignment, "fasta")

# Find and save the reference sequence from the aligned fasta file
ref_seq = next(record for record in aln if record.id == ref_id)

# Save the index for the first and last amino acids in the reference sequence
first = next(i for i, res in enumerate(ref_seq.seq) if res != "-")
last = len(ref_seq.seq) - next(i for i, res in enumerate(reversed(ref_seq.seq)) if res != "-") - 1

# Create the correct trimmed sequence to be added later to the fasta file
trimmed_seq = []

for record in aln:
    if record.id != ref_id:
        if first >= br_prefix_n: # If there are more than br_prefix_n amino acids at the start of the sequence trim at that mark
            trimmed_seq.append(record.seq[first - br_prefix_n: last + br_prefix_c + 1])
        else: # If there are less than br_prefix_n extra amino acids than the start of the sequence is not trimmed
            trimmed_seq.append(record.seq[0: last + br_prefix_c + 1])
        trimmed_seq[0] = trimmed_seq[0].strip('-') # Remove all the - at the end and beginning of the trimmed sequence

# Create a list with the accession number of the trimmed protein and the trimmed sequence
trimmed_fasta = []
for record in aln:
    if record.id != ref_id:
         trimmed_fasta.append((record.id, trimmed_seq[0]))

# Save the file as a fasta file according to snakemake logic
with open(trimmed, "w") as out_f:
    for record_id, seq in trimmed_fasta:
        out_f.write(f">{record_id}\n{seq}\n")

