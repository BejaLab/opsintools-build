# Read the TM-align output file
tm = str(snakemake.input)
fasta = str(snakemake.output)
ref = snakemake.config["ref_pdb"]
query = snakemake.wildcards.pdb

with open(tm, 'r') as f:
    lines = f.readlines()

with open(fasta, 'w') as out:
    out.write(f">{ref}\n")
    out.write(f"{lines[-4].strip()}\n")
    out.write(f">{query}\n")
    out.write(f"{lines[-2].strip()}\n")
