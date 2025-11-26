# Add all the protein fasta files into a single file for clustering
rule all_fasta:
    input:
        expand("analysis/data/alphafold/{pdb}/{{dataset_name}}/{pdb}.fasta", pdb = alphafold_list),
        expand("analysis/data/pdb/{pdb}/{{dataset_name}}/{pdb}.fasta", pdb = pdb_list)
    output:
        "analysis/clustering/{dataset_name}/all_sequences.fasta"
    conda:
        "envs/seqkit.yaml"
    shell:
        "seqkit seq {input} -gm 1 > {output}"

# As the clusters created by cd-hit were shown to be too large for fast MSA, this rule will cluster with mmseqs to get larger clusters
rule cluster_mmseqs:
    input: 
        "analysis/clustering/{dataset_name}/all_sequences.fasta"
    output:
        tsv = "analysis/clustering/{dataset_name}/mmseqs_{ident}_cluster.tsv",
        rep = "analysis/clustering/{dataset_name}/mmseqs_{ident}_rep_seq.fasta",
        all = "analysis/clustering/{dataset_name}/mmseqs_{ident}_all_seqs.fasta"
    log:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}.log"
    params:
        prefix = lambda w, output: output[0].replace("_cluster.tsv", ""),
        seq_id = lambda w: float(w.ident) / 100,
        cov_mode = 0
    conda:
        "envs/mmseqs2.yaml"
    shadow:
        "minimal"
    threads:
        10
    shell:
        "mmseqs easy-cluster {input} {params.prefix} tmp --threads {threads} --min-seq-id {params.seq_id} --cov-mode {params.cov_mode} &> {log}"

# Create lists of sequences to remove and to add so that all the reference Rhodopsins from the list above will be representatives of their cluster - mmseqs
rule ready_rep_w_ref_mmseqs:
    input:
        pdb_list = pdb_list_file,
        fasta = "analysis/clustering/{dataset_name}/mmseqs_{ident}_rep_seq.fasta",
        tsv = "analysis/clustering/{dataset_name}/mmseqs_{ident}_cluster.tsv"
    output:
        rep_del = "analysis/clustering/{dataset_name}/mmseqs_{ident}_rep_del.txt",
        rep_add = "analysis/clustering/{dataset_name}/mmseqs_{ident}_ref_add.txt"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/force_ref_seq_mmseqs.py"

# Replace the representatives with the reference Rhodopsins according to previous rule - mmseqs
checkpoint final_fasta_mmseqs:
    input:
        original_fasta = "analysis/clustering/{dataset_name}/mmseqs_{ident}_rep_seq.fasta",
        all_fasta = "analysis/clustering/{dataset_name}/all_sequences.fasta",
        to_remove = "analysis/clustering/{dataset_name}/mmseqs_{ident}_rep_del.txt",
        to_add = "analysis/clustering/{dataset_name}/mmseqs_{ident}_ref_add.txt"
    output:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.fasta"
    conda:
        "envs/seqkit.yaml"
    shell:
        "cat <(seqkit grep -vf {input.to_remove} {input.original_fasta}) <(seqkit grep -f {input.to_add} {input.all_fasta}) > {output}"
