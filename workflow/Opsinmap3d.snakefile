from Bio import SeqIO
import json, re, os

configfile: "workflow/build_config.yaml"

# Set up the file with all the accession numbers of proteins to download from AlphaFold
with open("workflow/alphafold.txt") as ids:
    alphafold_list = [ line.strip() for line in ids ]

# Set up the file with all the accession numbers of reference Rhodopsins to download from RCSB
with open("workflow/pdb.txt") as ids:
    pdb_list = [ line.split()[0] for line in ids ]

ref_pdb = config["ref_pdb"]
ident = config["ident"]

# Define the final targets of the workflow
rule all:
    input:
        "data/microbial-rhodopsins/ref.json",
        "data/microbial-rhodopsins/reps",
        "data/microbial-rhodopsins/must.txt"

# Download AlphaFold pdb files according to the full list provided
rule dload_alphafold:
    output:
       "analysis/data/alphafold/{pdb}/struct.pdb"
    resources:
        af_dload = 1
    threads: 0
    shell:
        "wget -qO {output} https://alphafold.ebi.ac.uk/files/AF-{wildcards.pdb}-F1-model_v4.pdb"

# Download RCSB pdb files according to the Rhodopsin reference list provided
rule dload_pdb:
    output:
        "analysis/data/pdb/{pdb}/struct_raw.pdb"
    resources:
        pdb_dload = 1
    threads: 0
    shell:
        "wget -qO {output} https://files.rcsb.org/download/{wildcards.pdb}.pdb"

# Fix the RCSB pdb files
rule fix_pdb:
    input:
        model = "analysis/data/pdb/{pdb}/struct_raw.pdb",
        atom_map = "workflow/atom_map.json"
    output:
        "analysis/data/pdb/{pdb}/struct.pdb"
    conda:
        "envs/pymol.yaml"
    script:
        "scripts/fix_structures.py"

# Trim helix
rule trim_helix:
    input:
        "analysis/data/{database}/{pdb}/struct.pdb"
    output:
        "analysis/data/{database}/{pdb}/struct_helices.pdb"
    params:
        padding_n = config["helix_padding_n"],
        padding_c = config["helix_padding_c"]
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/trim_helix.py"

# Align each protein to the reference Rhodopsin (BR) using us-align
rule us_align:
    input:
        ref_p = f"analysis/data/pdb/{ref_pdb}/struct_helices.pdb",
        query_p = "analysis/data/{database}/{pdb}/struct.pdb"
    output:
        "analysis/data/{database}/{pdb}/aln_to_ref_tm.txt"
    conda:
        "envs/us_align.yaml"
    shell:
        "USalign {input.ref_p} {input.query_p} -do > {output}"

# Trim each fasta file according to the BR alignment and filter according to RMSD & alignment length.
rule prot_trim_filter:
    input:
        aln = "analysis/data/{database}/{pdb}/aln_to_ref_tm.txt",
        query_pdb = "analysis/data/{database}/{pdb}/struct.pdb",
        ref_pdb = f"analysis/data/pdb/{ref_pdb}/struct.pdb"
    output:
        fasta = "analysis/data/{database}/{pdb}/{pdb}.fasta",
        pdb = "analysis/data/{database}/{pdb}/{pdb}.pdb"
    params:
        pad_n = config["pad_n"],
        pad_c = config["pad_c"],
        max_missing_n = config["max_missing_n"],
        max_missing_c = config["max_missing_c"],
        ref_lys_pos = config["ref_lysine_pos"],
        max_rmsd = config["max_rmsd"],
        min_aln_len = config["min_align_len"],
        max_aln_len = config["max_align_len"],
        min_len = config["min_len"]
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/prot_trim_filter.py"

# Add all the protein fasta files into a single file for clustering, reference proteins as is while the rest as trimmed versions.
rule all_fasta:
    input:
        expand("analysis/data/alphafold/{pdb}/{pdb}.fasta", pdb = alphafold_list),
        expand("analysis/data/pdb/{pdb}/{pdb}.fasta", pdb = pdb_list)
    output:
        "analysis/clustering/all_sequences.fasta"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit seq {input} -m 1 > {output}"

# As the clusters created by cd-hit were shown to be too large for fast MSA, this rule will cluster with mmseqs to get larger clusters
rule cluster_mmseqs:
    input: 
        "analysis/clustering/all_sequences.fasta"
    output:
        tsv = "analysis/clustering/mmseqs_{ident}_cluster.tsv",
        rep = "analysis/clustering/mmseqs_{ident}_rep_seq.fasta",
        all = "analysis/clustering/mmseqs_{ident}_all_seqs.fasta"
    params:
        prefix = lambda w, output: output[0].replace("_cluster.tsv", ""),
        seq_id = lambda w: float(w.ident) / 100,
        cov_mode = 0
    conda:
        "envs/mmseqs2.yaml"
    threads:
        10
    shell:
        "mmseqs easy-cluster {input} {params.prefix} tmp --threads {threads} --min-seq-id {params.seq_id} --cov-mode {params.cov_mode}"

# Create lists of sequences to remove and to add so that all the reference Rhodopsins from the list above will be representatives of their cluster - mmseqs
rule ready_rep_w_ref_mmseqs:
    input:
        "analysis/clustering/mmseqs_{ident}_rep_seq.fasta",
        "analysis/clustering/mmseqs_{ident}_cluster.tsv"
    output:
        "analysis/clustering/temp/mmseqs_{ident}_rep_del.txt",
        "analysis/clustering/temp/mmseqs_{ident}_ref_add.txt"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/force_ref_seq_mmseqs.py"

# Replace the representatives with the refrerence Rhodopsins according to previous rule - mmseqs
checkpoint final_fasta_mmseqs:
    input:
        original_fasta = f"analysis/clustering/mmseqs_{ident}_rep_seq.fasta",
        all_fasta = f"analysis/clustering/all_sequences.fasta",
        to_remove = f"analysis/clustering/temp/mmseqs_{ident}_rep_del.txt",
        to_add = f"analysis/clustering/temp/mmseqs_{ident}_ref_add.txt"
    output:
        "analysis/clustering/temp/mmseqs_{ident}_reps.fasta"
    conda:
        "envs/tools.yaml"
    shell:
        "cat <(seqkit grep -vf {input.to_remove} {input.original_fasta}) <(seqkit grep -f {input.to_add} {input.all_fasta}) > {output}"

def rep_pdb_files(wildcards):
    fasta_file = checkpoints.final_fasta_mmseqs.get(ident = config['ident']).output[0]
    pdb_files = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        database = "pdb" if record.id in pdb_list else "alphafold"
        pdb_files.append(f"analysis/data/{database}/{record.id}/{record.id}.pdb")
    return pdb_files

def parse_opm(data):
    matches = re.findall('(\\d+)\\(\\s*(\\d+)-\\s*(\\d+)\\)', data['subunits'][0]['segment'])
    return { f"TM{tm}": [ int(start), int(end) ] for tm, start, end in matches }

rule database_pdb:
    input:
        rep_pdb_files
    output:
        directory("data/microbial-rhodopsins/reps")
    shell:
        "mkdir -p {output} && cp {input} {output}/"

rule fetch_opm:
    output:
        "analysis/data/pdb/{pdb}/opm.json"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/fetch_opm.py"

rule write_ref_data:
    input:
        pdb = f"analysis/data/pdb/{ref_pdb}/struct.pdb",
        opm = f"analysis/data/pdb/{ref_pdb}/opm.json"
    output:
        "data/microbial-rhodopsins/ref.json"
    run:
        res_pos = {}
        with open(input.pdb) as file:
            for line in file:
                if line.startswith('ATOM'):
                    res_pos[int(line[22:26])] = 1
        with open(input.opm) as file:
            opm = json.load(file)
        data = {
            "id": ref_pdb,
            "res_pos": list(res_pos.keys()),
            "lysine": config['ref_lysine_pos'],
            "tms": parse_opm(opm)
        }
        with open(output[0], 'w') as file:
            json.dump(data, file)

rule copy_pdb_list:
    input:
        "workflow/pdb.txt"
    output:
        "data/microbial-rhodopsins/must.txt"
    shell:
        "cut -f1 {input} > {output}"
