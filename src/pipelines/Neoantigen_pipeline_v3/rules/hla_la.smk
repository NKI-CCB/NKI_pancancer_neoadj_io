
rule index_bam:
    input:
        bam=config["output_folder"]+"/mapped_hla/{normal}_sorted_hg38.bam",
    output:
        bai=temp(config["output_folder"]+"/mapped_hla/{normal}_sorted_hg38.bam.bai")
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools index {input}"


rule download_HLALA_graph:
   output:
       directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/PRG"),
       directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome"),
       directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/knownReferences"),
       directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/mapping"),
       directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/mapping_PRGonly"),
       directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/referenceGenomeSimulations"),
       directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/sampledReferenceGenomes"),
       directory("resources/graphs/PRG_MHC_GRCh38_withIMGT/translation"),
       "resources/graphs/PRG_MHC_GRCh38_withIMGT/sequences.txt",
   shell:
       "cd resources/graphs && wget  http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz "
       "&& tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz"


rule index_HLALA:
   input:
       "resources/graphs/PRG_MHC_GRCh38_withIMGT/sequences.txt",
   output:
       "resources/graphs/PRG_MHC_GRCh38_withIMGT/serializedGRAPH",
       "resources/graphs/PRG_MHC_GRCh38_withIMGT/serializedGRAPH_preGapPathindex",
   conda:
       "../envs/hla_la.yaml"
   params:
       path=lambda wc, input: os.path.dirname(os.path.dirname(input[0])),
       graph=lambda wc, input: os.path.basename(os.path.dirname(input[0]))
   shell:
       "HLA-LA.pl --prepareGraph 1 --customGraphDir {params.path} --graph {params.graph}"

# Copy graph - get HLA-LA running without conda
#rule copy_graph:
#    input:
#        "resources/graphs/PRG_MHC_GRCh38_withIMGT/sequences.txt"
#    output:
#        "copied_graph.txt"
#    conda:
#        "../envs/hla_la.yaml"
#    shell:
#        #"cp -r ./resources/graphs $CONDA_PREFIX/opt/hla-la && touch copied_graph.txt"
#        r"""
#        mkdir -p $CONDA_PREFIX/opt/hla-la/graphs        
#        cp -r resources/graphs/PRG_MHC_GRCh38_withIMGT $CONDA_PREFIX/opt/hla-la/graphs/ && touch copied_graph.txt
#        """

rule HLA_LA:
    input:
        bam=config["output_DNA_folder"]+"/mapped/{normal}_sorted_hg38_ARRG_dedup_recal.bam",
        bai=config["output_DNA_folder"]+"/mapped/{normal}_sorted_hg38_ARRG_dedup_recal.bam.bai",
        graph="copied_graph.txt",
 #       index="/DATA/j.traets/Pipelines/Neoantigen_snakemake_pipeline_draft/concat_draft/resources/graphs/PRG_MHC_GRCh38_withIMGT/serializedGRAPH_preGapPathIndex"
    output:
        config["output_folder"]+"/HLA-LA/output/{normal}/hla/R1_bestguess.txt",
    threads: 10
    params:
        out_dir = config["output_folder"]+"/HLA-LA/output",
        sample_id = replace_sample_id,
        #replace_folder =  replace_folder_id,
        #new_folder =  new_folder_id,
    conda:
        "../envs/hla_la.yaml"
    priority: 50
    shell:
        #"HLA-LA.pl --BAM {input.bam} --sampleID {params.sample_id} --graph PRG_MHC_GRCh38_withIMGT --workingDir {params.out_dir} --maxThreads {threads}"
        r"""
        # Run HLA-LA into the normalized directory
        HLA-LA.pl \
            --BAM {input.bam} \
            --sampleID {params.sample_id} \
            --graph PRG_MHC_GRCh38_withIMGT \
            --workingDir {params.out_dir} \
            --maxThreads {threads}

        # After HLA-LA finishes, copy the result into the expected Snakemake output path
        sample_dir="{params.out_dir}/{params.sample_id}/hla"	
        final="{output}"
	      echo "Sample dir: where HLA output is located, where to copy bestguess from"
        echo "$$sample_dir"

        # Try to find the full-resolution bestguess file
	if [ -f "$$sample_dir/R1_bestguess.txt" ]; then
	    bestguess="$$sample_dir/R1_bestguess.txt"
	    echo "INFO: Using full resolution"
            cp "$$bestguess" "{output}"
        # If missing, fall back to G-group typing
        elif [ -f "$$sample_dir/R1_bestguess_G.txt" ]; then
            bestguess="$$sample_dir/R1_bestguess_G.txt"
            echo "WARNING: Using G-group typing for {params.sample_id}"
            cp "$$bestguess" "{output}"
         # If neither exists, skip gracefully
        else
            echo "WARNING: No bestguess file found for {params.sample_id}, skipping"
            touch "$$final"
        fi
        """

rule parse_HLA_LA:
    input:
        config["output_folder"]+"/HLA-LA/output/{normal}/hla/R1_bestguess.txt"
    output:
        config["output_folder"]+"/HLA-LA/hlaI_{normal}.tsv",
        config["output_folder"]+"/HLA-LA/hlaII_{normal}.tsv"
    priority: 50
    script:
        "parse_HLA_types.py"
