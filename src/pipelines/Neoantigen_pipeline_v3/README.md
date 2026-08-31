# Neoantigen pipeline written in Snakemake 
# updated 2026-05

## Draft version v3
### Includes rules for HLA-typing (HLA-LA), annotating vcf files (VEP), neoantigen calling on SNVs (GATK) and indels (GATK/Strelka) using pVACtools (MHCflurry)

### Requires:
- GATK output
- STAR mapped bam files

### Pvactools has two modes:
- strict: >0.25 VAF (disabled RNA transcripts)
- relaxt: >0.05 VAF (disabled RNA transcripts)

### Before running:
mamba env create -f environment.yml
mamba create
conda activate neoantigen

And adjust the config file accordingly

### Run pipeline:
conda run -n neoantigen --no-capture-output python Run_Neoantigenpipeline.py


