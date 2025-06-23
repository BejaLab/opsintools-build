import re
from classes.USalign import USalign

# Asign file names from snakemake
alignment_files = snakemake.input
output_file = str(snakemake.output)

must_reps = snakemake.params['must_reps']
all_reps = snakemake.params['all_reps']

n_reps = int(snakemake.wildcards['n_reps'])
min_seq_id = snakemake.params['min_seq_id']

# Create an empty dictionary for RMSD valuse per representative
rmsd_dict = {}

def get_value(line, prefix):
    if match := re.search(prefix + '=\\s*([0-9.]+)', line):
        return float(match.group(1))

# Load all the input files
for rep, aln_file in zip(all_reps, alignment_files):
    # Read the aln file
    aln = USalign(aln_file)
    if aln.seq_id >= min_seq_id:
        rmsd_dict[rep] = aln.rmsd

# Sort the rmsd dictrionary by rmsd score
sorted_rmsd = dict(sorted(rmsd_dict.items(), key = lambda x: x[1]))

with open(output_file, 'w') as out_file:
    for rep in must_reps:
        out_file.write(f"{rep}\n")
    n = len(must_reps)
    for rep in sorted_rmsd.keys():
        if n >= n_reps: break
        if rep not in must_reps:
            out_file.write(f"{rep}\n")
            n += 1
