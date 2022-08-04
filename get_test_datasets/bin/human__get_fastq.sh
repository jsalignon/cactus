#!/bin/bash

specie="human"
prepro_dir="preprocessing/${specie}"
source $get_test_datasets_dir/get_test_datasets_functions.sh

mkdir -p $specie/data/mrna $specie/data/atac $specie/conf $specie/design

# creating the run.config file
cp run.config $specie/conf
sed -i "3s/^/\n  specie = 'human'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chromatin_state = 'ENCFF941SVR'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chip_ontology = 'cell_type.fibroblast'\n/" $specie/conf/run.config
cat $specie/conf/run.config

## details on the cell line:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98758
 # we measured the chromatin accessibility landscape using ATAC-seq following mock treatment, SSRP1 knockdown, or SUPT16H knockdown in human fibroblasts 
# cell type: human secondary fibroblasts
# genotype/variation: hiF-T cells carrying DOX-inducible, polycistronic human OCT4/KLF4/c-MYC/SOX2 (OKMS) cassette
# passages/stage: 13-18
# hiF-T cells -> derive from hBJ fibroblasts (= cell line established from skin taken from the normal foreskin of a neonatal male) https://www.sciencedirect.com/science/article/pii/S009286741500700X

## details on the chromatin state file:
# ENCFF941SVR: ChromHMM 18-state model of BSS00066: AG09309 from donor(s) ENCDO002AAA, cell line,	fibroblast,	skin of body, connective tissue

## alternative parameters if the cell are considered more like stem cell:
# sed -i "5s/^/\n  chip_ontology = 'cell_type.stem_cell'\n/" $specie/conf/run.config
# sed -i "5s/^/\n  chromatin_state = 'ENCFF676VUR'\n/" $specie/conf/run.config
# details on the chromatin state file: ENCFF676VUR, ChromHMM 18-state model of BSS00735: iPS-11a male adult (36 years) from donor(s) ENCDO632AGT,iPS-11a	cell line,	stem cell, induced pluripotent stem cell, skin of body


# downloading our fastq samples of interest and subsampling them
nextflow run nf-core/fetchngs --input "$samples_ids_dir/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 
# --force_sratools_download # => using this options results in files of the format "SRR7101009_R1.fastq.gz" instead of "SRX2794538_SRR5521297_R1.fastq.gz" which crash my parsing script

# creating the sample_info file
make_samples_info_file ${prepro_dir}

# renaming files
rename -v 's/SRX/atac_SRX/' ${prepro_dir}/fastq/SRX2794*
rename -v 's/SRX/mrna_SRX/' ${prepro_dir}/fastq/SRX4029*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' ${prepro_dir}/fastq/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' ${prepro_dir}/fastq/*
