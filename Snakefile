import yaml
config = yaml.safe_load(open("config.yaml"))

# Check snakemake dry run of what would be run
# snakemake -n

# Create a directed acyclic graph (DAG) showing the scripts, order of execution and their relationships
# snakemake --forceall --dag | dot -Tpdf > snakemake_execution_dag.pdf

# Run all
rule targets:
  input:
    "reports/parse_trial_metadata.html",
    "reports/gene_signatures.html",
    "reports/gene_signatures_bellini.html",
    "reports/gene_signatures_imcision.html",
    "reports/gene_signatures_nabucco.html",
    "reports/gene_signatures_niche_msi.html",
    "reports/gene_signatures_niche_mss.html",
    "reports/gene_signatures_opacin_neo.html",
    "reports/gene_signatures_prado.html",
    "reports/gene_signatures_meta.html",
    "reports/pancan_WES_analysis.html",
    "reports/pancan_TMB_response_rate_regression.html",
    "reports/pancan_WES_meta.html",
    "reports/pancan_WES_RNA_integration.html",
    "reports/pancan_glmer_sensitivity.html",
    "reports/pancan_figures.html"

# -------------------------------------------
# RNAseq rules

rule gene_signatures:
  input:
    raw_data = config["raw_data_dir"]
  output:
    "reports/gene_signatures.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/gene_signatures.Rmd $PWD/reports/gene_signatures.html \
    --output_dir "gene_signatures"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: ""
    """
    
rule gene_signatures_bellini:
  input:
    raw_data = config["raw_data_dir"]
  output:
    "reports/gene_signatures_bellini.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/gene_signatures.Rmd $PWD/reports/gene_signatures_bellini.html \
    --output_dir "gene_signatures_bellini"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: "BELLINI"
    """
    
rule gene_signatures_imcision:
  input:
    raw_data = config["raw_data_dir"]
  output:
    "reports/gene_signatures_imcision.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/gene_signatures.Rmd $PWD/reports/gene_signatures_imcision.html \
    --output_dir "gene_signatures_imcision"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: "IMCISION"
    """
 
rule gene_signatures_matisse:
  input:
    raw_data = config["raw_data_dir"]
  output:
    "reports/gene_signatures_matisse.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/gene_signatures.Rmd $PWD/reports/gene_signatures_matisse.html \
    --output_dir "gene_signatures_matisse"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: "MATISSE"
    """
    
rule gene_signatures_nabucco:
  input:
    raw_data = config["raw_data_dir"]
  output:
    "reports/gene_signatures_nabucco.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/gene_signatures.Rmd $PWD/reports/gene_signatures_nabucco.html \
    --output_dir "gene_signatures_nabucco"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: "NABUCCO"
    """
    
    
rule gene_signatures_niche_msi:
  input:
    raw_data = config["raw_data_dir"]
  output:
    "reports/gene_signatures_niche_msi.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/gene_signatures.Rmd $PWD/reports/gene_signatures_niche_msi.html \
    --output_dir "gene_signatures_niche_msi"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: "NICHE-MSI"
    """
    
rule gene_signatures_niche_mss:
  input:
    raw_data = config["raw_data_dir"]
  output:
    "reports/gene_signatures_niche_mss.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/gene_signatures.Rmd $PWD/reports/gene_signatures_niche_mss.html \
    --output_dir "gene_signatures_niche_mss"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: "NICHE-MSS"
    """
    
rule gene_signatures_opacin_neo:
  input:
    raw_data = config["raw_data_dir"]
  output:
    "reports/gene_signatures_opacin_neo.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/gene_signatures.Rmd $PWD/reports/gene_signatures_opacin_neo.html \
    --output_dir "gene_signatures_opacin_neo"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: "OpACIN_neo"
    """
    
rule gene_signatures_prado:
  input:
    raw_data = config["raw_data_dir"]
  output:
    "reports/gene_signatures_prado.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/gene_signatures.Rmd $PWD/reports/gene_signatures_prado.html \
    --output_dir "gene_signatures_prado"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: "PRADO"
    """
    
rule gene_signatures_meta:
  input:
    "reports/gene_signatures.html",
    "reports/gene_signatures_bellini.html",
    "reports/gene_signatures_imcision.html",
    "reports/gene_signatures_nabucco.html",
    "reports/gene_signatures_niche_msi.html",
    "reports/gene_signatures_niche_mss.html",
    "reports/gene_signatures_opacin_neo.html",
    "reports/gene_signatures_prado.html"
  output:
    "reports/gene_signatures_meta.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/gene_signatures_meta.Rmd $PWD/reports/gene_signatures_meta.html \
    --output_dir "gene_signatures_meta"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: ""
    """

# -------------------------------------------

# WES rules
  
rule WES_analysis:
  output:
    "reports/pancan_WES_analysis.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/pancan_WES_analysis.Rmd $PWD/reports/pancan_WES_analysis.html \
    --output_dir "WES_analysis"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: ""
    """

rule TMB_response_rate_regression:
  input:
    "reports/pancan_WES_analysis.html",
  output:
    "reports/pancan_TMB_response_rate_regression.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/pancan_TMB_response_rate_regression.Rmd $PWD/reports/pancan_TMB_response_rate_regression.html \
    --output_dir "WES_analysis_TMB_response_regression"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: ""
    """
    
rule WES_analysis_meta:
  input:
    "reports/pancan_WES_analysis.html",
  output:
    "reports/pancan_WES_meta.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/pancan_WES_meta.Rmd $PWD/reports/pancan_WES_meta.html \
    --output_dir "WES_analysis_meta"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: ""
    """
    
# -------------------------------------------

# WES RNA integration rules
    
rule WES_RNA_integration:
  input:
    "reports/pancan_WES_analysis.html",
    "reports/gene_signatures.html"
  output:
    "reports/pancan_WES_RNA_integration.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/pancan_WES_RNA_integration.Rmd $PWD/reports/pancan_WES_RNA_integration.html \
    --output_dir "WES_RNA_integration"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: ""
    """

rule pancan_glmer_sensitivity:
  input:
    "reports/pancan_WES_analysis.html",
    "reports/gene_signatures.html"
  output:
    "reports/pancan_glmer_sensitivity.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/pancan_glmer_sensitivity.Rmd $PWD/reports/pancan_glmer_sensitivity.html \
    --output_dir "GLMM_analysis"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: ""
    """

# -------------------------------------------

# Figures rules
rule pancan_figures:
  input:
    "reports/gene_signatures_meta.html",
    "reports/pancan_WES_meta.html",
    "reports/pancan_TMB_response_rate_regression.html",
    "reports/pancan_WES_RNA_integration.html",
    "reports/pancan_glm_sensitivity.html",
    "reports/pancan_glmer_sensitivity.html"
  output:
    "reports/pancan_figures.html"
  shell:
    """
    Rscript src_knit_master.Rmd src/pancan_figures.Rmd $PWD/reports/pancan_figures.html \
    --output_dir "_pancan_figures"   \
    --metadata "combined_rna_metadata.csv"   \
    --export TRUE \
    --analyze_specific_trial: ""
    """

