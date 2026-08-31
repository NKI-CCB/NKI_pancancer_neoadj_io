rule copy_fastq:
    input:
        R1 = config["fastq_folder"] + "/{sample}" + config["read_name"]["R1"] + ".fastq.gz",
        R2 = config["fastq_folder"] + "/{sample}" + config["read_name"]["R2"] + ".fastq.gz"
    output:
        R1_cp = temp(config["output_folder"] + "/fastq_files/{sample}" + config["read_name"]["R1"] + ".fastq.gz"),
        R2_cp = temp(config["output_folder"] + "/fastq_files/{sample}" + config["read_name"]["R2"] + ".fastq.gz")
    shell:
        """
        cp {input.R1} {output.R1_cp}
        cp {input.R2} {output.R2_cp}
        """

rule bwa_mem:
    input:
        reads = temp([config["output_folder"] + "/fastq_files/{sample}" + config["read_name"]["R1"] + ".fastq.gz",config["output_folder"] + "/fastq_files/{sample}" + config["read_name"]["R2"] + ".fastq.gz"]),
        genome = config["params"]["hlala"]["ref"]
    output:
        temp(config["output_folder"]+"/mapped_hla/{sample}_hg38.bam") # temporary files
    conda:
        "../envs/mapping.yaml"
    threads:
        10
    shell:
        "bwa mem -M -t {threads} -T 0 {input.genome} {input.reads} | samtools view -Shb - > {output}"

rule samtools_sort: # coordinate
    input:
        temp(config["output_folder"]+"/mapped_hla/{sample}_hg38.bam")
    output:
        temp(config["output_folder"]+"/mapped_hla/{sample}_sorted_hg38.bam")
    conda:
        "../envs/mapping.yaml"
    threads:
        2
    shell:
        "samtools sort {input} > {output}"
