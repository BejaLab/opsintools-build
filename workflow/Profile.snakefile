num_chunks = 100
chunks = [ f"{i:02d}" for i in range(num_chunks) ]

rule templates:
    input:
        lambda w: rep_pdb_files(w.ident_templ)
    output:
        "analysis/clustering/{dataset_name}/templates_{ident_templ}.txt"
    run:
        with open(str(output), 'w') as file:
            for pdb_file in input:
                file.write(f">{Path(pdb_file).stem} _P_ {pdb_file}\n")

# Not used
rule align_reps_t_coffee:
    input:
        fasta = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.fasta",
        templates = "analysis/clustering/{dataset_name}/templates_{ident_templ}.txt",
        pdb = lambda w: rep_pdb_files(w.ident_templ)
    output:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_{ident_templ}_reps.aln"
    params:
        methods = 'mustang_pair,sap_pair,mafftlinsi_msa'
    log:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_{ident_templ}_reps.aln.log"
    conda:
        "envs/t_coffee.yaml"
    shadow:
        "minimal"
    threads:
        workflow.cores
    shell:
        "t_coffee {input.fasta} -outfile {output} -output aln,score_ascii -method {params.methods} -template_file {input.templates} -pdb_min_sim 90 -pdb_min_cov 0 -thread {threads} 2>{log}"

rule align_reps:
    input:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.fasta"
    output:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.fas"
    log:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.fas.log"
    conda:
        "envs/mafft.yaml"
    threads:
        workflow.cores
    shell:
        "mafft --localpair --maxiterate 1000 --thread {threads} {input} > {output} 2> {log}"

rule trim_aln:
    input:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.fas"
    output:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_trimmed.fas"
    conda:
        "envs/trimal.yaml"
    shell:
        """
        readarray -t COLS < <(trimal -in {input} -automated1 -out /dev/null -colnumbering | grep -o '[[:digit:]]*')
        trimal -selectcols {{ ${{COLS[0]}}-${{COLS[-1]}} }} -complementary -in {input} -out {output}
        """

rule aln_to_a3m:
    input:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.fas"
    output:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.a3m"
    params:
        M = 60
    conda:
        "envs/hhsuite.yaml"
    shell:
        "reformat.pl fas a3m {input} - -v 0 -M {params.M} | hhconsensus -i /dev/stdin -o /dev/stdout | grep -v '^#' | sed '2,$s/^>/>@/' > {output}"

rule hhsearch_chunk:
    input:
        a3m_ffdata  = "analysis/databases_split/{db}_{chunk}_{of}_a3m.ffdata",
        a3m_ffindex = "analysis/databases_split/{db}_{chunk}_{of}_a3m.ffindex",
        cs219_ffdata  = "analysis/databases_split/{db}_{chunk}_{of}_cs219.ffdata",
        cs219_ffindex = "analysis/databases_split/{db}_{chunk}_{of}_cs219.ffindex",
        hhm_ffdata  = "analysis/databases_split/{db}_{chunk}_{of}_hhm.ffdata",
        hhm_ffindex = "analysis/databases_split/{db}_{chunk}_{of}_hhm.ffindex",
        a3m = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.a3m"
    output:
        "analysis/hhsearch/{dataset_name}/chunks/mmseqs_{ident}_reps_{db}_{chunk}_{of}.hhr"
    log:
        "analysis/hhsearch/{dataset_name}/chunks/mmseqs_{ident}_reps_{db}_{chunk}_{of}.log"
    params:
        db = "analysis/databases_split/{db}_{chunk}_{of}",
        e = 1,
        BZ = 100000
    conda:
        "envs/hhsuite.yaml"
    threads:
        10
    shell:
        "hhsearch -all -glob -i {input.a3m} -d {params.db} -o {output} -e {params.e} -cpu {threads} -B {params.BZ} -Z {params.BZ} -maxmem inf &> {log}"

rule filter_hhsearch:
    input:
        hhr = expand("analysis/hhsearch/{{dataset_name}}/chunks/mmseqs_{{ident}}_reps_{{db}}_{chunk}_{of}.hhr", chunk = chunks, of = num_chunks),
        names = expand("analysis/databases_split/{{db}}_{chunk}_{of}_names.txt", chunk = chunks, of = num_chunks),
        ffindex = expand("analysis/databases_split/{{db}}_{chunk}_{of}_{{ext}}.ffindex", chunk = chunks, of = num_chunks)
    output:
        "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}_{ext}.ffindex"
    params:
        probab = config["probab_hhsearch"],
        prefix = "{db}_"
    conda:
        "envs/bioformats.yaml"
    script:
        "scripts/filter_hhsearch_records.py"

rule merge_ffindex:
    input:
        datas = expand("analysis/databases/{db}/{prefix}_{{ext}}.ffdata", zip, db = config['databases'], prefix = [ prefixes[d] for d in config['databases'] ]),
        indexes = expand("analysis/hhsearch/{{dataset_name}}/filtered/mmseqs_{{ident}}_reps_{db}_{{ext}}.ffindex", db = config['databases'])
    output:
        data = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps_{ext}.ffdata",
        index = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps_{ext}.ffindex"
    conda:
        "envs/ffindex_py.yaml"
    shell:
        "ffindex_merge_py -k -d {output.data} -i {output.index} -- {input.datas} {input.indexes}"

rule hhsearch_filtered:
    input:
        a3m_ffindex = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps_a3m.ffindex",
        cs219_ffindex = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps_cs219.ffindex",
        hhm_ffindex = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps_hhm.ffindex",
        a3m_ffdata = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps_a3m.ffdata",
        cs219_ffdata = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps_cs219.ffdata",
        hhm_ffdata = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps_hhm.ffdata",
        a3m = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.a3m"
    output:
        hhr = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps.hhr",
        a3m = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps.a3m"
    log:
        "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps.log"
    params:
        db = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps",
        e = 1,
        M = 1000000
    conda:
        "envs/hhsuite.yaml"
    threads:
        10
    shell:
        "hhsearch -all -glob -i {input.a3m} -d {params.db} -o {output.hhr} -oa3m {output.a3m} -e 1 -cpu {threads} -B {params.M} -Z {params.M} -maxseq {params.M} &> {log}"

rule hhfilter:
    input:
        "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps.a3m"
    output:
        "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps_hhfilter.a2m"
    params:
        ident = config["ident_hmmbuild"],
        M = 1000000
    conda:
        "envs/hhsuite.yaml"
    shell:
        "hhfilter -id {params.ident} -maxseq {params.M} -i {input} -o /dev/stdout | addss.pl /dev/stdin /dev/stdout | reformat.pl a3m a2m /dev/stdin {output}"

rule a2m_trim:
    input:
        "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps_hhfilter.a2m"
    output:
        a2m = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps_hhfilter_trimmed.a2m",
        a3m = "analysis/hhsearch/{dataset_name}/merged/mmseqs_{ident}_reps_hhfilter_trimmed.a3m"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/trim_helix_a2m.py"

rule copy_a3m:
    input:
        lambda w: "analysis/hhsearch/{{dataset_name}}/merged/mmseqs_{ident}_reps_hhfilter_trimmed.a3m".format(ident = config['ident_hhsearch'])
    output:
        "output/{dataset_name}/profile.a3m"
    shell:
        "cp {input} {output}"

rule hmmbuild:
    input:
        lambda w: "analysis/hhsearch/{{dataset_name}}/merged/mmseqs_{ident}_reps_hhfilter_trimmed.a2m".format(ident = config['ident_hhsearch'])
    output:
        "output/{dataset_name}/profile.hmm"
    params:
        args = '--wnone'
    conda:
        "envs/hmmer.yaml"
    shell:
        "hmmbuild {params.args} -n {wildcards.dataset_name} {output} {input}"

rule hhsearch_ref:
    input:
        profile = "output/{dataset_name}/profile.hmm",
        fasta = f"analysis/data/pdb/{ref_pdb}/struct.pdb.fasta"
    output:
        "output/{dataset_name}/profile_ref.txt"
    conda:
        "envs/hmmer.yaml"
    shell:
        "hmmsearch --nobias -o {output} {input.profile} {input.fasta}"

rule hmm_consensus:
    input:
        "output/{dataset_name}/profile.hmm"
    output:
        "output/{dataset_name}/profile.fasta"
    conda:
        "envs/hmmer.yaml"
    shell:
        "hmmemit -c -o {output} {input}"
