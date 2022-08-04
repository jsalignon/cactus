#!/bin/bash

specie="mouse"
prepro_dir="preprocessing/${specie}"
source $get_test_datasets_dir/get_test_datasets_functions.sh

mkdir -p $specie/data/mrna $specie/data/atac $specie/conf $specie/design

# creating the run.config file
cp run.config $specie/conf
sed -i "3s/^/\n  specie = 'mouse'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chromatin_state = 'ENCFF809HLK'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chip_ontology = 'all'\n/" $specie/conf/run.config
cat $specie/conf/run.config

# ENCFF809HLK	mm10	ChromHMM 18 state model for kidney (postnatal 0 days), mesoderm,	excretory system,	mouse

# downloading our fastq samples of interest and subsampling them
nextflow run nf-core/fetchngs --input "$samples_ids_dir/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 -resume

# creating a simple reference file
make_samples_info_file ${prepro_dir}

# renaming files
cat ${prepro_dir}/samplesheet/samples_info.csv
rename -v 's/SRX/mrna_SRX/' ${prepro_dir}/fastq/SRX117086{63..78}*
rename -v 's/SRX/atac_SRX/' ${prepro_dir}/fastq/SRX117086{79..90}*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' ${prepro_dir}/fastq/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' ${prepro_dir}/fastq/*
ls ${prepro_dir}/fastq
