rule CNVkit_access:
    input:
        bam_files = expand(config["output_folder"]+"/mapped/{sample}_sorted_hg38_ARRG_dedup_recal.bam",sample=samples["sample_ID"]),
        genome = config["params"]["CNVkit"]["ref"]
    output:
        access = config["output_folder"]+"/CNVkit/access.hg38.bed"
    singularity:
        config["SIF"]["cnvkit"]
    log:
        config["output_folder"] + "/logs/CNVkit/access_CNVkit_run.log"
    shell:
        """
        cnvkit.py access {input.genome} -o {output.access} 2> {log}
        """

rule CNVkit_all:
    input:
        bam_tumor = expand(config["output_folder"]+"/mapped/{tumor}_sorted_hg38_ARRG_dedup_recal.bam",tumor=samples_matched["Tumor"]),
        bam_normal = expand(config["output_folder"]+"/mapped/{normal}_sorted_hg38_ARRG_dedup_recal.bam",normal=samples_matched["Normal"]),
        genome = config["params"]["CNVkit"]["ref"],
        ref = config["output_folder"]+"/CNVkit/results_normal/normals_reference.cnn",
        access = config["output_folder"]+"/CNVkit/access.hg38.bed",
    params:
        bait = config["params"]["CNVkit"]["bait_bed"],
        refflat = config["params"]["CNVkit"]["refflat"],
        output_f = config["output_folder"]+"/CNVkit_test/results/"
    output:
        out1 = config["output_folder"]+"/CNVkit_test/results"
    singularity:
        config["SIF"]["cnvkit"]
    log:
        config["output_folder"]+"/logs/CNVkit_test/CNVkit_run.log"
    shell:
        """
        cnvkit.py batch {input.bam_tumor} --normal {input.bam_normal} --targets {params.bait} --annotate {params.refflat} --fasta {input.genome} --access {input.access} --output-reference {input.ref} --output-dir {params.output_f} --diagram --scatter  2> {log}
        """

rule CNVkit_purity:
    input:
        cns = config["output_folder"]+"/CNVkit/results/{tumor}_sorted_hg38_ARRG_dedup_recal.cns",
        facets = config["output_folder"] + "/FACETS/fitted/{tumor}-vs-{normal}.snppile.csv.gz_fitted.csv"
    output:
        config["output_folder"]+"/CNVkit/results/{tumor}_{normal}/{tumor}_sorted_hg38_ARRG_dedup_recal_purity_thr.call.cns"
    params:
        purity = get_tumor_purity,
        ploidy = get_tumor_ploidy,
        thresholds = get_threshold
    singularity:
        config["SIF"]["cnvkit"]
    shell:
        """
        cnvkit.py call {input.cns} -y -m threshold -t={params.thresholds} --purity {params.purity} -o {output}
        """

rule export_CNVkit:
    input:
        config["output_folder"]+"/CNVkit/results/{tumor}_sorted_hg38_ARRG_dedup_recal.cnr"
    output:
        config["output_folder"]+"/CNVkit/results/{tumor}_sorted_hg38_ARRG_dedup_recal.seg"
    singularity:
        config["SIF"]["cnvkit"]
    shell:
        """
        cnvkit.py export seg {input} --enumerate-chroms -o {output}
        """
