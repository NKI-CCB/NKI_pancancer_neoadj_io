rule run_pvac_strict:
    input:
        config["output_folder"]+"/annotated/{tumor}-vs-{normal}_variants.annotated.vcf" 
    output:
        config["output_folder"]+"/pvactools/{tumor}-vs-{normal}_strict/MHC_Class_I/{tumor}.filtered.tsv"
    params:
        path_table = get_hla,
        mhc_folder = config["output_folder"]+"/pvactools/{tumor}-vs-{normal}_strict/",
#        iedb = "/DATA/j.traets/Tools/pvactools/"
    log:
        config["output_folder"]+"/logs/pvac/{tumor}-vs-{normal}_strict.log"
    conda:
        "../envs/pvac.yaml"
    threads:
        4
    shell:
        "pvacseq run {input} {wildcards.tumor} {params.path_table} MHCflurry {params.mhc_folder} -t 4"

rule run_pvac_binding_f:
    input:
        config["output_folder"]+"/pvactools/{tumor}-vs-{normal}/MHC_Class_I/{tumor}.filtered.tsv"
    output:
        config["output_folder"]+"/pvactools/{tumor}-vs-{normal}/MHC_Class_I/{tumor}.filtered_bf.tsv"
    conda:
        "../envs/pvac.yaml"
    threads:
        4
    shell:
        "pvacseq binding_filter {input} {output}"

rule run_pvac_relaxt:
    input:
        vcf = config["output_folder"] + "/annotated/{tumor}-vs-{normal}_variants.annotated.vcf",
        hla = config["output_folder"]+"/HLA-LA/hlaI_{normal}.tsv"
    output:
        config["output_folder"] + "/pvactools/{tumor}-vs-{normal}_relaxt/MHC_Class_I/{tumor}.filtered.tsv"
    params:
        alleles = get_hla,   # MUST return comma-separated HLA alleles
        mhc_folder = config["output_folder"] + "/pvactools/{tumor}-vs-{normal}_relaxt/",
        vaf_threshold = float(config["params"]["pvactools"]["vaf_threshold"])
    log:
        config["output_folder"] + "/logs/pvac/{tumor}-vs-{normal}.log"
    conda:
        "../envs/pvac.yaml"
    threads:
        4
    shell:
        "pvacseq run "
        "{input.vcf} "
        "{wildcards.tumor} "
        "{params.alleles} "
        "MHCflurry "
        "{params.mhc_folder} "
        "-t 4 --tdna-vaf {params.vaf_threshold}"

