# Genomic and transcriptomic correlates of neoadjuvant immunotherapy response across cancer types
# Author: Pedro Batista Tan
This repository contains analysis code for the pancancer neoadjuvant immune checkpoint blockade (ICB) project from the Netherlands Cancer Institute (NKI).
A renv.lock file is provided to track packages and environment used in this R project.
R Analysis code is provided in the source "src" folder.
A snakemake workflow is provided to track which corresponding scripts should be run and the order of execution for the main analysis.
The snakemake workflow knits each script in order, using the main 

An overview of the execution of scripts is provided in the "snakemake_execution_dag.pdf"
DNA and RNA preprocessing pipelines are written as snakemake workflows provided in /src/pipelines/.

System Requirements:
Analyses were performed on a Linux server environment (Ubuntu 24.04.3 LTS).
There are no requirements for specialized, non-standard hardware.
Analyses were executed on a server with 56 CPU cores and 128 GB RAM.

Computations were carried out using R version 4.2.3 within RStudio.
The workflow requires the following R packages (specific versions also listed in renv lock): 
- meta (v8.1), lme4 (v1.1-37). ggplot2 (v3.5.2), tidyverse (v2.0.0), ggpubr (v0.6.0), ggpattern (v1.1.4), pheatmap (v1.0.13), corrplot (v0.95), export (v0.3.0)
- DNA sequencing: bcl2fastq (v2.20), Skewer (v0.2.2), BWA (v0.7.17), Picard (v2.27.5), FastQC (v0.11.9), Qualimap (v2.2.1), Mosdepth (v0.3.3), NGSCheckMate (v1.0.1), MultiQC (v1.13), GATK4 (v4.2.2.0), GATK4 Mutect2 (v4.4.0), Variant Effect Predictor (VEP) (v105), vcf2maf (v1.6.17), FACETS (v0.6.2), MutationalPatterns (v3.8.1), Maftools (v2.14.0), VariantAnnotation (v1.44.1)
- RNA sequencing: bcl2fastq (v2.20), Skewer (v0.2.2), STAR (v2.7.9), HTSeq (v2.0.2), Kallisto (v0.48.0), FastQC (v0.11.9), ngsderive (v1.2.0), and MultiQC (v1.8), DESeq2 (v.1.38.3), GSVA package (v1.46), matrixStats (v1.1.0)


Instalation guide:
Instalation of R and R packages can be tipically done within a day.
To install R and Rstudio, see https://www.r-project.org/ and https://posit.co/download/rstudio-desktop/
For packages, run install.packages("[package_name]") or BiocManager::install("[package_name]") for Bioconductor packages.
The code expects an .Rproj which is and treated as the origin path with the package "here".
This project should be activated with "renv" and can be rebuilt with information from the renv.lock file provided.
The main project structure should contain a "data/processed/" and "reports/" folder to receive the output from source (src/) scripts.

Demo:
The snakemake workflow contains the information to run the main analysis.
We have provided a demo for signature analysis and meta-analysis of this study.
*This demo does not contain any real data, only mock data*
This is available for the gene_signatures.Rmd and gene_signatures_meta.Rmd, that would source the "data/processed_demo/" folder, provided with the code. 
When loading objects, the paths have to be adapted to use this data/processed_demo/ as input.
This can be toggled with the "use_demo" parameter in YAML headers, and will change the directory from which the data would be sourced.

The output folder can be changed within scripts, or declared as parameters within the knit snakemake command line calls.
A data/processed folder should be created manually to store the expected output.
Expected output is folders within the data/processed/ folder, which would contain the output for corresponding scripts in the snakefile, and a reports/ folder.
Expected run time for the demo data should be complete within 1 day.

## License
This project is licensed under the MIT License – see the LICENSE file for details.