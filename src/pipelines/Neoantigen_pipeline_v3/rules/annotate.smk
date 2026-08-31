
rule annotate_variants:
    input:
        calls = config["output_DNA_folder"] + "/GATK_out/{tumor}-vs-{normal}-GATKFiltered-pass.vcf",
        plugins = "vep/plugins",
        # optionally add reference genome fasta
        # fasta = config["params"]["bwa"]["ref"]
    output:
        calls = config["output_folder"] + "/annotated/{tumor}-vs-{normal}_variants.annotated.vcf"
    params:
        cache = config["params"]["vep"]["cache_dir"],
        fasta = config["params"]["vep"]["fasta"],
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins = ["Frameshift", "Wildtype"],
        plugin_dir = config["params"]["vep"]["plugin_dir"]
    conda:
        "../envs/annotate.yaml"
    threads: 4
    shell:
      """
	export PERL5LIB={params.plugin_dir}:${{PERL5LIB:-}}

        vep \
            --format vcf \
            --input_file {input.calls} \
            --output_file {output.calls} \
            --dir_cache {params.cache} \
            --fasta {params.fasta} \
            --offline \
            --vcf \
            --symbol \
            --terms SO \
            --tsl \
            --hgvs \
            --plugin Frameshift,{params.plugin_dir}/Frameshift.pm \
            --plugin Wildtype,{params.plugin_dir}/Wildtype.pm \
            --force_overwrite
            """
