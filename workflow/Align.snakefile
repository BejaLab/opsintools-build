canon_microb_rhodopsin_reps ,= glob_wildcards("output/canon_microb_rhodopsin/reps/{rep}.pdb")
heliorhodopsin_reps ,= glob_wildcards("output/heliorhodopsin/reps/{rep}.pdb")

rule opsinalign3d:
    input:
        expand("output/canon_microb_rhodopsin/reps/{rep}.pdb", rep = canon_microb_rhodopsin_reps),
        expand("output/heliorhodopsin/reps/{rep}.pdb", rep = heliorhodopsin_reps)
    output:
        directory("analysis/microb_rhodopsins")
    params:
        methods = "mtm_align_msa"
    conda:
        "envs/opsintools.yaml"
    threads:
        10
    shell:
        "opsinalign3d -i {input} -o {output} --methods {params.methods} -t {threads}"

rule parse_align:
    input:
        aln = "analysis/microb_rhodopsins",
        refs = expand("output/{profile}/ref.json", profile = [ "canon_microb_rhodopsin", "heliorhodopsin" ])
    output:
        "output/microb_rhodopsins_align.tsv"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/parse_align.py"
