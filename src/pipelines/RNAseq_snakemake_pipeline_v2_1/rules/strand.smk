rule index_bam:
    input:
        bam = config["output_folder"] + "STAR_mapped/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        config["output_folder"] + "STAR_mapped/{sample}_Aligned.sortedByCoord.out.bam.bai"
    singularity: "sif/samtools_1.9_staphb_2305.sif"
    shell:
        """
        samtools index {input.bam}
        """

rule check_strand:
    input:
        bam = expand(config["output_folder"] + "STAR_mapped/{sample}_Aligned.sortedByCoord.out.bam",sample=samples["sample_ID"]),
        bai = expand(config["output_folder"] + "STAR_mapped/{sample}_Aligned.sortedByCoord.out.bam.bai",sample=samples["sample_ID"])
    output:
        config["output_folder"] + "ngsderive/strandedness.txt"
    singularity: "sif/ngsderive_v1.2.0_stjudecloud_2305.sif"
    params:
        genome_dir = config["ngsderive"]["reference"]
    shell:
        """
        ngsderive strandedness {input.bam} -g {params.genome_dir} > {output}
        """
