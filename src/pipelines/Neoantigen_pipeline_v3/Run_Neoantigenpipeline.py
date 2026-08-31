#                              Neoantigen pipeline                             #
#                    ### script to run Neoantigen pipeline ###                 #
# 
# Run in the corresponding neoantigen conda env 
# mamba env create -f ./environment.yml
# conda activate neoantigen

import sys
import subprocess
import argparse
import os
import yaml
import pandas as pd
import logging
from datetime import datetime

# config file location fixed
def load_yaml():
    with open("./config/config.yaml", "r") as stream:
        try:
            safe_config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return safe_config

def check_path_format():
    if yaml_config["fastq_folder"][-1] != "/":
        yaml_config["fastq_folder"] = yaml_config["fastq_folder"]+"/"
    if yaml_config["output_folder"][-1] != "/":
        yaml_config["output_folder"] = yaml_config["output_folder"]+"/"
  
# run snakemakepipeline with OS command based on user input
def run_snakemake(yaml_config,config):
    # prepare bind volumes
    bind_d1 = yaml_config["fastq_folder"]
    bind_d1 = bind_d1.split("/")
    bind_d1.pop()
    bind_d1 = "/".join(bind_d1)
    
    bind_d2 = yaml_config["output_DNA_folder"]
    bind_d2 = bind_d2.split("/")
    bind_d2.pop()
    bind_d2 = "/".join(bind_d2)
    
    bind_d3 = yaml_config["output_folder"]
    bind_d3 = bind_d3.split("/")
    bind_d3.pop()
    bind_d3 = "/".join(bind_d3)

    CHECK_FOLDER = os.path.isdir(bind_d3)
    if not CHECK_FOLDER:
        os.makedirs(bind_d3)
    
    bind_d4 = yaml_config["params"]["hlala"]["ref"]
    bind_d4 = bind_d4.split("/")
    bind_d4.pop()
    bind_d4 = "/".join(bind_d4)

    log.info("bind singularity:")
    log.info(bind_d1)
    log.info(bind_d2)
    log.info(bind_d3)
    log.info(bind_d4)
    
    logging.shutdown()
    
    # save log file and config file to output folder
    log_config = bind_d3+"/logs/"+file_path[:-4]
    CHECK_FOLDER = os.path.isdir(log_config)
    if not CHECK_FOLDER:
        os.makedirs(log_config)
    os.system("mv " + file_path + " " + log_config + "/" + file_path)
    os.system("cp config/config.yaml " + log_config + "/")

    # run snakemake
    os.system('snakemake --use-conda --conda-frontend conda --use-singularity --singularity-args "-B ' + bind_d1 + ' -B ' + bind_d2  +  ' -B ' +  bind_d3 + ' -B ' +  bind_d4 + ' " --cores ' + str(config["cores"]) + ' ' + str(yaml_config["snakemake"]))
#snakemake --use-conda --conda-frontend conda --use-singularity --singularity-args " -B /DATA/Datasets/PanCancer/Preprocessed/pipelines/Neoantigen_pipeline_draft_v2/test_data -B/DATA/Datasets/PanCancer/Preprocessed/pipelines/Neoantigen_pipeline_draft_v2/test_output -B /DATA/Datasets/PanCancer/Preprocessed/PRADO/DNAseq_output_rerun" --cores 1

if __name__ == "__main__":
  # prepare logging, stdout and log file
  now = datetime.now()
  file_path = now.strftime("%Y-%m-%d_%H:%M:%S.log")
  print(file_path)

  log = logging.getLogger()
  logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
  log.setLevel(logging.DEBUG)

  fh = logging.FileHandler(file_path)
  fh.setFormatter(logFormatter)
  log.addHandler(fh)

  ch = logging.StreamHandler()
  ch.setFormatter(logFormatter)
  log.addHandler(ch)
  
  parser = argparse.ArgumentParser(description="Snakemake Neoantigen - Pipeline",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-c", "--cores",default=1, help="Number of cores to be used (int)")
  args = parser.parse_args()
  config = vars(args)
  
  # print input user
  log.info("\n \n ### User input ### ")
  log.info("Nr of cores: " + str(config["cores"]))
  
  yaml_config = load_yaml()
  check_path_format()
  
  # print most important parameters in config file
  log.info("\n \n ### Parameters in config.yaml file ### ")
  log.info("Sample file: " + yaml_config["samples"])
  log.info("Read name 1: *" + yaml_config["read_name"]["R1"] + ".fastq.gz")
  log.info("Read name 2: *" + yaml_config["read_name"]["R2"] + ".fastq.gz")
  log.info("Reference: " + yaml_config["params"]["hlala"]["ref"])

  run_snakemake(yaml_config,config)
