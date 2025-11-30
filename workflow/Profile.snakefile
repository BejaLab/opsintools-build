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
        # opsinalign3d -i {input} -o {output} -t {threads} --methods {params.methods}

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
        a3m = expand("analysis/databases_split/{{db}}_{chunk}_{of}_a3m.ffindex", chunk = chunks, of = num_chunks),
        cs219 = expand("analysis/databases_split/{{db}}_{chunk}_{of}_cs219.ffindex", chunk = chunks, of = num_chunks),
        hhm = expand("analysis/databases_split/{{db}}_{chunk}_{of}_hhm.ffindex", chunk = chunks, of = num_chunks)
    output:
        a3m = "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}_a3m.ffindex",
        cs219 = "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}_cs219.ffindex",
        hhm = "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}_hhm.ffindex"
    params:
        probab = 92
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/filter_hhsearch_records.py"

rule link_filtered_ffdata:
    input:
        lambda w: "analysis/databases/{db}/{prefix}_{ext}.ffdata".format(db = w.db, prefix = dbs[w.db], ext = w.ext)
    output:
        "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}_{ext}.ffdata"
    shell:
        "ln -sr {input} {output}"

rule hhsearch_filtered:
    input:
        a3m_ffindex = "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}_a3m.ffindex",
        cs219_ffindex = "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}_cs219.ffindex",
        hhm_ffindex = "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}_hhm.ffindex",
        a3m_ffdata = "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}_a3m.ffdata",
        cs219_ffdata = "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}_cs219.ffdata",
        hhm_ffdata = "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}_hhm.ffdata",
        a3m = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.a3m"
    output:
        hhr = "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}.hhr",
        a3m = "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}.a3m"
    log:
        "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}.log"
    params:
        db = "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}",
        e = 1,
        BZ = 100000
    conda:
        "envs/hhsuite.yaml"
    threads:
        10
    shell:
        "hhsearch -all -glob -i {input.a3m} -d {params.db} -o {output.hhr} -oa3m {output.a3m} -e 1 -cpu {threads} -B {params.BZ} -Z {params.BZ} &> {log}"

rule hhfilter:
    input:
        "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}.a3m"
    output:
        "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}.a2m"
    params:
        ident = config["ident_hmmbuild"]
    conda:
        "envs/hhsuite.yaml"
    shell:
        "hhfilter -id {params.ident} -i {input} -o /dev/stdout | reformat.pl a3m a2m - - -v 0 | seqkit grep -vrp _consensus$ -o {output}"

rule hmmbuild:
    input:
        lambda w: "analysis/hhsearch/{dataset_name}/filtered/mmseqs_{ident}_reps_{db}.a2m".format(dataset_name = w.dataset_name, db = w.db, ident = config['ident_hhsearch'])
    output:
        "output/{dataset_name}/profile_{db}_{type}.hmm"
    params:
        args = lambda w: { 'opt': '--wnone', 'def': '' }[w.type]
    conda:
        "envs/hmmer.yaml"
    shell:
        "hmmbuild {params.args} -n {wildcards.dataset_name} {output} {input}"
