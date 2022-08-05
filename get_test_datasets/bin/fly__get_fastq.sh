#!/bin/bash

specie="fly"
prepro_dir="preprocessing/${specie}"
fastq_dir=${prepro_dir}/fastq

source $get_test_datasets_bin_dir/get_test_datasets_functions.sh

# making directory structure
mkdir -p $specie/data/mrna $specie/data/atac $specie/conf $specie/design

# cleaning up the fastq folder
if [ -d $fastq_dir ]; then rm -r $fastq_dir ; fi

# downloading fastq files
nextflow run nf-core/fetchngs --input "$samples_ids_dir/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 -resume # -process.cache "deep"
# using the option '-process.cache "deep"' gives this error:
# No such variable: Exception evaluating property 'single_end' for java.util.ArrayList, Reason: groovy.lang.MissingPropertyException: No such property: single_end for class: sun.nio.fs.UnixPath
#  -- Check script '/home/jersal/.nextflow/assets/nf-core/fetchngs/./workflows/sra.nf' at line: 76 or see '.nextflow.log' file for more details

# checking sample details to rename them
make_samples_info_file ${prepro_dir}

# renaming files
rename -v 's/SRX/mrna_SRX/' $fastq_dir/SRX81740{44..53}*
rename -v 's/SRX/atac_SRX/' $fastq_dir/SRX81740{34..43}*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' $fastq_dir/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' $fastq_dir/*
ls $fastq_dir

