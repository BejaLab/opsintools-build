configfile: "config.yaml"

# Set up the file with all the accession numbers of proteins to download from AlphaFold
with open("workflow/alphafold.txt") as ids:
    alphafold_list = [line.strip() for line in ids]

# Set up the file with all the accession numbers of reference Rhodopsins to download from RCSB
with open("workflow/pdb.txt") as ids:
    pdb_list = [line.split()[0] for line in ids]

# Set up lists according to the config file
c_list = config["c"]
methods = config["prws_aln_meth"]
ref_pdb = config["ref_pdb"]

# Define the final targets of the workflow
rule all:
    input:
        expand("analysis/data/alphafold/{pdb}/struct.pdb", pdb = alphafold_list),
        expand("analysis/data/alphafold/{pdb}/sequence.fasta", pdb = alphafold_list),
        expand("analysis/clustering/{c}_clustered_sequences.fasta", c = c_list),
        expand("analysis/clustering/{c}_clustered_sequences.fasta.clstr", c = c_list),
        "analysis/clustering/clustered_sequences.fasta"

# Download AlphaFold pdb files according to the full list provided
rule dload_alphafold:
    output:
       "analysis/data/alphafold/{pdb}/struct.pdb"
    resources:
        af_dload = 1
    shell:
        "wget -qO {output} https://alphafold.ebi.ac.uk/files/AF-{wildcards.pdb}-F1-model_v4.pdb"

# Download RCSB pdb files according to the Rhodopsin reference list provided
rule dload_pdb:
    output:
        "analysis/data/pdb/{pdb}/struct_raw.pdb"
    resources:
        pdb_dload = 1
    shell:
        "wget -qO {output} https://files.rcsb.org/download/{wildcards.pdb}.pdb"

# Fix the RCSB pdb files - remove the ligand atoms
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

# Rename the pdb files for later rules to the accession number
rule rename_pdb:
    input:
        "analysis/data/{database}/{pdb}/struct.pdb"
    output:
        "analysis/data/{database}/{pdb}/{pdb}.pdb"
    shell:
        "cp {input} {output}"

# Align each protein to the reference Rhodopsin (BR) using FoldMason
rule foldmason_aln:
    input:
        ref = f"analysis/data/pdb/{ref_pdb}/{ref_pdb}.pdb",
        query = "analysis/data/{database}/{pdb}/{pdb}.pdb"
    output:
        "analysis/data/{database}/{pdb}/aln_to_br_aa.fa"
    params:
        prefix = lambda w, output: str(output).replace('_aa.fa', '')
    conda:
        "envs/foldmason.yaml"
    shell:
        "foldmason easy-msa {input.query} {input.ref} {params.prefix} {resources.tmpdir} --report-mode 2"

# Convert pdb files to fasta
rule pdb_to_fasta:
    input:
        "analysis/data/{database}/{pdb}/{pdb}.pdb"
    output:
        "analysis/data/{database}/{pdb}/raw_sequence.fasta"
    shell:
        "pdb_tofasta {input} > {output} && sed -i '1s/.*/>{wildcards.pdb}/' {output}"

# Fix fasta files that still got ligands in the sequence as X from previous rule
rule fix_fasta:
    input:
        "analysis/data/{database}/{pdb}/raw_sequence.fasta"
    output:
        "analysis/data/{database}/{pdb}/{pdb}.fasta"
    shell:
        "sed '2,$ s/X//g' {input} | sed '/^$/d' > {output}"

# Align each protein to the reference Rhodopsin (BR) using t-coffee
rule t_coffee_aln:
    input:
        ref = f"analysis/data/pdb/{ref_pdb}/{ref_pdb}.fasta",
        query = "analysis/data/{database}/{pdb}/{pdb}.fasta"
    output:
        "analysis/data/{database}/{pdb}/aln_to_br_t_coffee.fa"
    conda:
        "envs/t_coffee.yaml"
    shell:
        "t_coffee -in {input.ref} {input.query} -output fasta_aln -outfile {output}"

rule t_coffee_pairwise_template:
    input:
        ref_p = f"analysis/data/pdb/{ref_pdb}/{ref_pdb}.pdb",
        query_p = "analysis/data/{database}/{pdb}/{pdb}.pdb"
    output:
        "analysis/data/{database}/{pdb}/t_coffee_input.txt"
    run:
        with open(str(output), 'w') as out:
            out.write(f">{ref_pdb}  _P_ {input.ref_p}\n")
            out.write(f">{wildcards.pdb}  _P_ {input.query_p}\n")

# Align each protein to the reference Rhodopsin (BR) using t-coffee 3d alignment
rule t_coffee_3d_aln:
    input:
        ref_p = f"analysis/data/pdb/{ref_pdb}/{ref_pdb}.pdb",
        query_p = "analysis/data/{database}/{pdb}/{pdb}.pdb",
        ref_s = f"analysis/data/pdb/{ref_pdb}/{ref_pdb}.fasta",
        query_s = "analysis/data/{database}/{pdb}/{pdb}.fasta",
        template = "analysis/data/{database}/{pdb}/t_coffee_input.txt"
    output:
        "analysis/data/{database}/{pdb}/aln_to_br_coffee_3d.fa"
    params:
        methods = "sap_pair,mustang_pair,t_coffee_msa,probcons_msa"
    conda:
        "envs/t_coffee.yaml"
    shell:
        "t_coffee {input.ref_s} {input.query_s} -output fasta_aln -outfile {output} -method {params.methods} -template_file {input.template} -pdb_min_sim 0 -pdb_min_cov 0 -n_core 1"

# Align each protein to the reference Rhodopsin (BR) using tm-align
rule tm_aln:
    input:
        ref_p = f"analysis/data/pdb/{ref_pdb}/{ref_pdb}.pdb",
        query_p = "analysis/data/{database}/{pdb}/{pdb}.pdb"
    output:
        "analysis/data/{database}/{pdb}/aln_to_br_tm.aln"
    conda:
        "envs/tm_aln.yaml"
    shell:
        "TMalign {input.ref_p} {input.query_p} > {output}"

# Create a fasta file from the aln file that was created by the tm_align rule
rule tm_aln_to_fasta:
    input:
        "analysis/data/{database}/{pdb}/aln_to_br_tm.aln"
    output:
        "analysis/data/{database}/{pdb}/aln_to_br_tm.fa"
    script:
        "scripts/tm_aln_to_fasta.py"

# For the list of reference rhodopsins run tm-algin with benchmark directive to check resource management
rule resource_check_tm_align:
    input:
        ref_p = f"analysis/data/pdb/{ref_pdb}/{ref_pdb}.pdb",
        query_p = "analysis/data/pdb/{pdb}/{pdb}.pdb"
    output:
        "analysis/data/pdb/{pdb}/bnch_aln_to_br_tm.aln"
    conda:
        "envs/tm_aln.yaml"
    benchmark:
        repeat("analysis/data/pdb/{pdb}/tm_benchmakr.tsv", 100)
    shell:
        "TMalign {input.ref_p} {input.query_p} > {output}"

# For the list of reference rhodopsins run foldmason with benchmark directive to check resource management.
rule resource_check_foldmason:
    input:
        ref = f"analysis/data/pdb/{ref_pdb}/{ref_pdb}.pdb",
        query = "analysis/data/pdb/{pdb}/{pdb}.pdb"
    output:
        "analysis/data/pdb/{pdb}/bnch_aln_to_br_aa.fa"
    params:
        prefix = lambda w, output: str(output).replace('_aa.fa', '')
    conda:
        "envs/foldmason.yaml"
    benchmark:
        repeat("analysis/data/pdb/{pdb}/fm_benchmark.tsv", 100)
    shell:
        "foldmason easy-msa {input.query} {input.ref} {params.prefix} {resources.tmpdir} --report-mode 2"

# Trim each fasta file according to the BR alignment. br_prefix_n amino acids before TM1 (see in config.yaml) amino acids and br_prefix_c amino acids after TM7
rule trim_prot:
    input:
        "analysis/data/{database}/{pdb}/aln_to_br_tm.fa"
    output:
        "analysis/data/{database}/{pdb}/trimmed.fasta"
    params:
        ref=ref_pdb,
        br_prefix_n=config["br_prefix_n"],
        br_prefix_c=config["br_prefix_c"]
    script:
        "scripts/trim_prot.py"

# Add all the protein fasta files into a single file for clustering, reference proteins as is while the rest as trimmed versions.
rule all_fasta:
    input:
        expand("analysis/data/alphafold/{pdb}/trimmed.fasta", pdb = alphafold_list),
        expand("analysis/data/pdb/{pdb}/trimmed.fasta", pdb = pdb_list)
    output:
        "analysis/clustering/all_raw_sequences.fasta"
    shell:
        "cat {input} >> {output}"

# Remove from the all_raw_sequences.fasta file all the sequences that have less than 180 amino acids
rule remove_short:
    input:
         "analysis/clustering/all_raw_sequences.fasta"
    output:
         "analysis/clustering/all_sequences.fasta"
    params:
         n=config["min_len"]
    shell:
         "seqkit seq {input} -g -m {params.n} -o {output}"

# Cluster all the sequences using cd-hit with diffrerent -c values for analysis
rule cluster_fasta:
    input:
        "analysis/clustering/all_sequences.fasta"
    output:
        "analysis/clustering/{c}_clustered_sequences.fasta",
        "analysis/clustering/{c}_clustered_sequences.fasta.clstr"
    shell:
        "cd-hit -i {input} -o {output[0]} -c 0.{wildcards.c} -n 2"

# Create a list of accession numbers recieved from previous rule
rule cluster_analysis:
    input:
        "analysis/clustering/{c}_clustered_sequences.fasta.clstr"
    output:
        "analysis/clustering/cluster_analysis.txt"
    params:
        c=c_list
    shell:
        "grep -c '^>' {input} > {output}"

# Create lists of sequences to remove and to add so that all the reference Rhodopsins from the list above will be representatives of their cluster
rule ready_rep_w_ref:
    input:
        "analysis/clustering/5_clustered_sequences.fasta"
    output:
        "analysis/clustering/temp/representative_to_del.txt",
        "analysis/clustering/temp/references_to_add.txt"
    script:
        "scripts/force_ref_seq.py"

# Replace the representatives with the refrerence Rhodopsins according to previous rule
rule replace_rep_with_ref:
    input:
        original_fasta="analysis/clustering/5_clustered_sequences.fasta",
        all_fasta="analysis/clustering/all_sequences.fasta",
        to_remove="analysis/clustering/temp/representative_to_del.txt",
        to_add="analysis/clustering/temp/references_to_add.txt"
    output:
        "analysis/clustering/clustered_sequences.fasta"
    shell:
        """
        seqkit grep -v -f {input.to_remove} {input.original_fasta} -o temp_removed.fasta
        seqkit grep -f {input.to_add} {input.all_fasta} -o temp_added.fasta
        cat temp_removed.fasta temp_added.fasta > {output}
        rm -f temp_removed.fasta temp_added.fasta
        """
# Analysis of the length of the clusters
rule cluster_length:
    input:
        "analysis/clustering/clustered_sequences.fasta"
    output:
        "analysis/clustering/clustered_length_analysis.txt"
    shell:
        "seqkit fx2tab -l {input} | tr -s ' ' '\t' | cut -f1,3 | sort -k2,2n -t$'\t' > {output}"

# Analysis of the length of all the fasta files
rule all_fasta_length:
    input:
        "analysis/clustering/all_sequences.fasta"
    output:
        "analysis/clustering/all_length_analysis.txt"
    shell:
        "seqkit fx2tab -l {input} | tr -s ' ' '\t' | cut -f1,3 | sort -k2,2n -t$'\t' > {output}"

# Create a file with all the accession numbers of the clustered proteins
rule all_aln_proteins:
    input:
        "analysis/clustering/clustered_sequences.fasta"
    output:
        "analysis/clustering/clustered_proteins.txt"
    shell:
        "seqkit fx2tab {input} | tr -s ' ' '\t' | cut -f1 > {output}"

