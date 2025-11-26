# Download AlphaFold pdb files according to the full list provided
rule dload_alphafold:
    output:
       "analysis/data/alphafold/{pdb}/struct.pdb"
    params:
        version = 'v6'
    resources:
        af_dload = 1,
    threads: 0
    shell:
        "wget -qO {output} https://alphafold.ebi.ac.uk/files/AF-{wildcards.pdb}-F1-model_{params.version}.pdb"

# Fix the RCSB pdb files
rule dload_pdb:
    input:
        atom_map = "workflow/atom_map.json"
    output:
        "analysis/data/pdb/{pdb}/struct.pdb"
    conda:
        "envs/pymol.yaml"
    shadow:
        "minimal"
    script:
        "scripts/dload_and_clean_pdb.py"

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

# Align each protein to the reference rhodopsin (e.g. BR) using us-align
rule us_align:
    input:
        ref_p = "analysis/data/pdb/{ref_pdb}/struct_helices.pdb",
        query_p = "analysis/data/{database}/{pdb}/struct.pdb"
    output:
        "analysis/data/{database}/{pdb}/aln_to_ref_tm_{ref_pdb}.txt"
    conda:
        "envs/us_align.yaml"
    shell:
        "USalign {input.ref_p} {input.query_p} -do > {output}"

# Trim each fasta file according to the ref alignment and filter according to RMSD & alignment length.
rule prot_trim_filter:
    input:
        query_pdb = "analysis/data/{database}/{pdb}/struct.pdb",
        aln = f"analysis/data/{{database}}/{{pdb}}/aln_to_ref_tm_{ref_pdb}.txt",
        ref_pdb = f"analysis/data/pdb/{ref_pdb}/struct_helices.pdb"
    output:
        fasta = "analysis/data/{database}/{pdb}/{dataset_name}/{pdb}.fasta",
        pdb = "analysis/data/{database}/{pdb}/{dataset_name}/{pdb}.pdb"
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

# Fetch data from OPM
rule fetch_opm:
    output:
        "analysis/data/pdb/{pdb}/opm.json"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/fetch_opm.py"

rule dload_uniclust:
    output:
        directory("analysis/uniclust")
    params:
        url = "https://wwwuser.gwdguser.de/~compbiol/uniclust/2023_02/UniRef30_2023_02_hhsuite.tar.gz"
    shell:
        "mkdir -p {output} && wget -O- {params.url} | tar xvz -C {output}"

rule dload_bfd:
    output:
        directory("analysis/BFD")
    params:
        url = "https://bfd.mmseqs.com/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz"
    shell:
        "mkdir -p {output} && wget -O- {params.url} | tar xvz -C {output}"
