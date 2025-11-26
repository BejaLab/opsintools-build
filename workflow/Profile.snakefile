dbs = { 'BFD': 'bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt', 'uniclust': 'UniRef30_2023_02' }
chunks = 100

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

rule hhblits_bfd:
    input:
        db_dir = "analysis/databases/BFD",
        a3m = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.a3m"
    output:
        hhr = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_BFD_hhblits.hhr",
        a3m = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_BFD_hhblits.a3m"
    log:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_BFD_hhblits.log"
    params:
        db = lambda w, input: f"{input.db_dir}/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
        e = 1e-5,
        n = 3
    conda:
        "envs/hhsuite.yaml"
    threads:
        10
    shell:
        "hhblits -i {input.a3m} -d {params.db} -o {output.hhr} -oa3m {output.a3m} -n {params.n} -e {params.e} -id {wildcards.ident} -cpu {threads} &> {log}"

rule hhsearch_uniclust:
    input:
        db_dir = "analysis/databases/uniclust",
        a3m = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.a3m"
    output:
        hhr = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_uniclust.hhr",
        a3m = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_uniclust.a3m"
    log:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_uniclust.log"
    params:
        db = lambda w, input: f"{input.db_dir}/UniRef30_2023_02",
        e = 1e-5,
        BZ = 1000,
        maxres = 10000,
        maxmem = 150
    conda:
        "envs/hhsuite.yaml"
    threads:
        1
    shell:
        "hhsearch -mark -glob -i {input.a3m} -d {params.db} -o {output.hhr} -oa3m {output.a3m} -e {params.e} -cpu {threads} -B {params.BZ} -Z {params.BZ} -maxres {params.maxres} -maxmem {params.maxmem} &> {log}"

rule hhsearch_bfd:
    input:
        db_dir = "analysis/databases/BFD",
        a3m = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.a3m"
    output:
        hhr = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_BFD.hhr",
        a3m = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_BFD.a3m"
    log:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_BFD.log"
    params:
        db = lambda w, input: f"{input.db_dir}/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
        e = 1e-5,
        BZ = 1000000
    conda:
        "envs/hhsuite.yaml"
    threads:
        5
    shell:
        "hhsearch -all -glob -i {input.a3m} -d {params.db} -o {output.hhr} -oa3m {output.a3m} -e {params.e} -id {wildcards.ident} -cpu {threads} -B {params.BZ} -Z {params.BZ} &> {log}"

rule a3m_to_a2m:
    input:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_BFD.a3m"
    output:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_BFD.a2m"
    conda:
        "envs/hhsuite.yaml"
    shell:
        "reformat.pl a3m a2m {input} {output}"

rule trim_a2m:
    input:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_BFD.a2m"
    output:
        "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_BFD_trimmed.a2m"
    run:
        with open(str(output), 'w') as file:
            for record in SeqIO.parse(str(input), "fasta"):
                if record.id != 'ss_pred' and record.id != 'ss_conf' and not record.id.endswith('_consensus'):
                    SeqIO.write(record, file, 'fasta')

rule database_num_records:
    input:
        "analysis/databases/{db}"
    output:
        "analysis/databases/{db}_num_records.txt"
    params:
        ffindex = lambda w, input: Path(str(input)) / (dbs[w.db] + '_a3m.ffindex')
    shell:
        "wc -l < {params.ffindex} > {output}"

rule split_ffindex:
    input:
        db_dir = "analysis/databases/{db}",
        num_records = "analysis/databases/{db}_num_records.txt"
    output:
        directory("analysis/databases_subsets/{db}_{chunk}_{of}")
    params:
        input_prefix = lambda w, input: Path(input.db_dir) / dbs[w.db],
        output_prefix = lambda w, output: Path(str(output)) / dbs[w.db]
    script:
        "scripts/split_ffindex.py"

rule hhsearch_chunk:
    input:
        db_dir = "analysis/databases_subsets/{db}_{chunk}_{of}",
        a3m = "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps.a3m"
    output:
        hhr = "analysis/hhsearch/{dataset_name}/mmseqs_{ident}_reps_{db}_{chunk}_{of}.hhr",
        a3m = "analysis/hhsearch/{dataset_name}/mmseqs_{ident}_reps_{db}_{chunk}_{of}.a3m"
    log:
        "analysis/hhsearch/{dataset_name}/mmseqs_{ident}_reps_{db}_{chunk}_{of}.log"
    params:
        db = lambda w, input: Path(input.db_dir) / dbs[w.db],
        e = 1e-5,
        BZ = 1000000
    conda:
        "envs/hhsuite.yaml"
    threads:
        10
    resources:
        hhsearch = 1
    shell:
        "hhsearch -all -glob -i {input.a3m} -d {params.db} -o {output.hhr} -oa3m {output.a3m} -e {params.e} -id {wildcards.ident} -cpu {threads} -B {params.BZ} -Z {params.BZ} -maxmem inf &> {log}"

rule hhsearch_collect:
    input:
        expand("analysis/hhsearch/{{dataset_name}}/mmseqs_{{ident}}_reps_{{db}}_{chunk}_{of}.a3m", chunk = [ f"{i:02d}" for i in range(chunks) ], of = chunks)
    output:
        "analysis/hhsearch/{dataset_name}/mmseqs_{ident}_reps_{db}.a3m"
    shell:
        "cat {input} > {output}"

rule hmmbuild:
    input:
        lambda w: "analysis/clustering/{dataset_name}/mmseqs_{ident}_reps_BFD_trimmed.a2m".format(dataset_name = w.dataset_name, ident = config['ident_align'])
    output:
        "output/{dataset_name}/profile_{type}.hmm"
    params:
        args = lambda w: { 'opt': '--fragthresh 0 --wblosum', 'def': '' }[w.type]
    conda:
        "envs/hmmer.yaml"
    shell:
        "hmmbuild {params.args} -n {wildcards.dataset_name} {output} {input}"
