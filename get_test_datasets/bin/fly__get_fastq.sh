#!/bin/bash

specie="fly"
prepro_dir="preprocessing/${specie}"

mkdir -p $specie/data/mrna $specie/data/atac $specie/conf $specie/design

# creating the run.config file
cp run.config $specie/conf
sed -i "3s/^/\n  specie = 'fly'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chromatin_state = 'iHMM.M1K16.fly_L3'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chip_ontology = 'all'\n/" $specie/conf/run.config
cat $specie/conf/run.config

# downloading our fastq samples of interest and subsampling them
nextflow run nf-core/fetchngs --input "$samples_ids_dir/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 -resume

# creating a simple reference file
make_samples_info_file ${prepro_dir}
cat ${prepro_dir}/samplesheet/samples_info.tsv

# renaming files
cat ${prepro_dir}/samplesheet/samples_info.csv
rename -v 's/SRX/mrna_SRX/' ${prepro_dir}/fastq/SRX81740{44..53}*
rename -v 's/SRX/atac_SRX/' ${prepro_dir}/fastq/SRX81740{34..43}*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' ${prepro_dir}/fastq/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' ${prepro_dir}/fastq/*
ls ${prepro_dir}/fastq

