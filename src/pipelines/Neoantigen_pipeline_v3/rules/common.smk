# common script
import glob
import itertools
import pandas as pd

configfile: "config/config.yaml"

samples = pd.read_table(config["samples"],sep=";")
#print(set(samples["sample_ID"]))

# function to create new data frame, tumor vs normal
def creat_matched(samples):
    norm_list = []
    tumor_list = []
    rna_list = []
    for x_p in set(samples["Patient"]):
        sel_pat = samples[samples["Patient"]==x_p]
        norm_list.append(list(list(sel_pat[sel_pat["Tumor_yes"]=="NO"]["sample_ID"])*len(sel_pat[sel_pat["Tumor_yes"]=="YES"]["sample_ID"])))
        tumor_list.append(list(sel_pat[sel_pat["Tumor_yes"]=="YES"]["sample_ID"]))
#        rna_list.append(list(sel_pat[sel_pat["Tumor_yes"]=="YES"]["RNA"]))
    norm_list = list(itertools.chain(*norm_list))
    tumor_list = list(itertools.chain(*tumor_list))
#    rna_list = list(itertools.chain(*rna_list))
#    matched_col = pd.DataFrame(list(zip(norm_list,tumor_list,rna_list)),columns =['Normal', 'Tumor','RNA'])
    matched_col = pd.DataFrame(list(zip(norm_list,tumor_list)),columns =['Normal', 'Tumor'])
    return(matched_col)

samples_matched = creat_matched(samples)
print("Check data frame!")
print(samples_matched)

# collect output
def output_rules_neoantigen():
    samples_matched = creat_matched(samples)

    prepare_graph_hlala = "copied_graph.txt"

    bam_output = expand(config["output_folder"]+"/mapped_hla/{sample}_sorted_hg38.bam",sample=set(samples["sample_ID"])),
    hla_guess = expand(config["output_folder"]+"/HLA-LA/output/{normal}/hla/R1_bestguess.txt",normal=samples_matched["Normal"]),
    hla_i = expand(config["output_folder"]+"/HLA-LA/hlaI_{normal}.tsv",normal=samples_matched["Normal"]),
    hla_ii = expand(config["output_folder"]+"/HLA-LA/hlaII_{normal}.tsv",normal=samples_matched["Normal"]),

    annot = expand(config["output_folder"]+"/annotated/{tumor}-vs-{normal}_variants.annotated.vcf",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"]),
    annot_indel = expand(config["output_folder"]+"/annotated/{tumor}-vs-{normal}-indels-variants.annotated.vcf",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"]),

    mhc_relaxt = expand(config["output_folder"]+"/pvactools/{tumor}-vs-{normal}_relaxt/MHC_Class_I/{tumor}.filtered_bf.tsv",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"]),
    mhc_strict = expand(config["output_folder"]+"/pvactools/{tumor}-vs-{normal}_strict/MHC_Class_I/{tumor}.filtered_bf.tsv",zip,tumor=samples_matched["Tumor"],normal=samples_matched["Normal"])
    return mhc_relaxt

def get_hla(wildcards):
    hla_file = f"{config['output_folder']}/HLA-LA/hlaI_{wildcards.normal}.tsv"
    hla_table = pd.read_table(hla_file)

    alleles = []
    for x in hla_table["Allele"]:
        # Convert HLA-LA format (e.g. HLA-A0205) → pvacseq format (HLA-A*02:05)
        gene = x[:5]          # HLA-A
        field1 = x[5:7]       # 02
        field2 = x[7:9]       # 05
        alleles.append(f"{gene}*{field1}:{field2}")

    return ",".join(alleles)

def replace_sample_id(wildcards):
    print(wildcards)
    temp_name = str(wildcards)
    return(temp_name.replace("-","_"))

def replace_folder_id(wildcards):
    temp_name = str(wildcards)
    temp_name = temp_name.replace("-","_")
    return(config["output_folder"]+"/HLA-LA/output/"+temp_name+"/*")

def new_folder_id(wildcards):
    temp_name = str(wildcards)
    return(config["output_folder"]+"/HLA-LA/output/"+temp_name)

