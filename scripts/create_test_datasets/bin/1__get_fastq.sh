#!/bin/bash

species=$1

prepro_dir="preprocessing/${species}"
fastq_dir=${prepro_dir}/fastq


# making directory structure
mkdir -p $species/data/mrna $species/data/atac $species/parameters $species/design

# cleaning up the fastq folder
if [ -d $fastq_dir ]; then rm -r $fastq_dir ; fi

# downloading fastq files
nextflow run nf-core/fetchngs --input "$samples_ids_dir/srr_accession/srr_${species}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 -resume

# checking sample details to rename them
make_samples_info_file ${prepro_dir}

# renaming files
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' $fastq_dir/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' $fastq_dir/*
ls $fastq_dir
