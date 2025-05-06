import re

# Asign file names from snakemake
alignments = snakemake.input.ref_aln
list_reps_to_aln = snakemake.output.rep_lists
must_reps = snakemake.params.must_reps
n_reps = snakemake.params.n_reps

# Create an empty dictionary for RMSD valuse per representative
rmsd_dict = {}

# Load all the input files
for aln in alignments:
    # Read the aln file
    with open(aln, 'r') as in_file:
        # For each line in the alignment file
        for line in in_file:
            # Find and save the representative
            if "Name of Chain_1" in line:
                rep = re.search(r'/([^/]+)/[^/]+$', line).group(1)
            # Find and save the RMSD score
            if "RMSD=" in line:
                score = float(re.search(r'RMSD=\s*([\d.]+)', line).group(1))
                # Save the representative and the score it got to the dictionary
                rmsd_dict[rep] = score
                break

# Sort the rmsd dictrionary by rmsd score
sorted_rmsd = dict(sorted(rmsd_dict.items(), key=lambda x: x[1]))

# For each n create an list with the top n files in sorted dict
for n, out_f in zip(n_reps, list_reps_to_aln):
    output_list = list(sorted_rmsd.keys())[:n]
    
    # Remove any Must Have reps in order to add them manually later. Bugged out if stay in file.
    output_list = [rep for rep in output_list if rep not in must_reps]

    # For each Must Have Rep add it at the top of the list
    for rep in must_reps:
        output_list.insert(0, rep)

    # Take only the top n items in the list
    output_list = output_list[:n]

    # Create and save the list of reps to use in the alignment
    with open(out_f, 'w') as out_file:
        for rep in output_list:
            out_file.write(f"{rep}\n")
