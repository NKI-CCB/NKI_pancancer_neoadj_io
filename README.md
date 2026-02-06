# Genomic and transcriptomic correlates of neoadjuvant immunotherapy response across cancer types
# Author: Pedro Batista Tan

This repository contains analysis code for the pancancer neoadjuvant immune checkpoint blockade (ICB) project from the Netherlands Cancer Institute (NKI).
A renv.lock file is provided to track packages and environment used in this R project.
R Analysis code is provided in the source "src" folder.
A snakemake workflow is provided to track which corresponding scripts should be run and the order of execution for the main analysis.

DNA and RNA preprocessing pipelines are written as snakemake workflows provided in /src/pipelines/