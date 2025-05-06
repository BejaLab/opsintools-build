from Bio import AlignIO
import json

aln = snakemake.input[0]
out_file = snakemake.output[0]
query = snakemake.wildcards.query
ref = snakemake.config["ref_pdb"]
offset = 9 # The offset for reference 7Z09 first residue in alignment
membrane_range = snakemake.config["membrane_range"]

# Load the positions of the protein which are in the membrane
membrane_pos = []
for start, end, tm_label in membrane_range:
    for i in range(start, end + 1):
        membrane_pos.append((i, tm_label))

# Load the alignment file and filter out all irrelevant sequences
alignment = AlignIO.read(aln, "clustal")
filtered_seq_records = [record for record in alignment if record.id in [query, ref]]
residues = list(zip(filtered_seq_records[1].seq, filtered_seq_records[0].seq))

# Create a dictonary with a map of the aligned sequences with the correct residue numeration
counter = 1
mapping = {}
for i, (r, q) in enumerate(residues):
    if r != '-':
        mapping[counter + offset] = (r, q)
        counter += 1

# Validate the lysine in the 7th TM is in the right position
if mapping[216][0] == 'K' and mapping[216][1] == 'K':
    # Create a new dictionary only with the positions inside the membrane
    tm_map = {}
    for pos, tm_label in membrane_pos:
        if pos in mapping:
            tm_map[pos] = {"residue_ref": mapping[pos][0], "residue_query": mapping[pos][1], "helix": tm_label}
else:
    tm_map = "Warning: No lysine at position 216 in the 7th TM."

# Add a description to the data that will be the json file
out_data = {
    "description": f"Mapping for rhodopsin residues inside the transmembrane. Key is position number: (residue in {ref}, residue in {query} and the number of the TM helix",
    "alignment_map": tm_map
}

# Save as json file
with open(out_file, 'w') as out_f:
    json.dump(out_data, out_f, indent=2)
