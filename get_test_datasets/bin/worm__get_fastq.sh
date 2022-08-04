#!/bin/bash

specie="worm"
prepro_dir="preprocessing/${specie}"
source $get_test_datasets_bin_dir/get_test_datasets_functions.sh

# making directory structure
mkdir -p $specie/data/mrna $specie/data/atac $specie/conf $specie/design

# cleaning up the fastq folder
fastq_dir=${prepro_dir}/fastq
if [ -d $fastq_dir ]; then rm -r $fastq_dir ; fi

# creating the run.config file
cp run.config $specie/conf
sed -i "3s/^/\n  specie = 'worm'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chromatin_state = 'iHMM.M1K16.worm_L3'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chip_ontology = 'all'\n/" $specie/conf/run.config
cat $specie/conf/run.config

# downloading fastq files
nextflow run nf-core/fetchngs --input "$samples_ids_dir/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 -resume

# creating a basic reference file
make_samples_info_file ${prepro_dir}
cat ${prepro_dir}/samplesheet/samples_info.tsv

# renaming files
rename -v 's/SRX/mrna_SRX/' ${fastq_dir}/SRX30291{12..20}*
rename -v 's/SRX/atac_SRX/' ${fastq_dir}/SRX30291{24..35}*
rename -v 's/SRX/atac_SRX/' ${fastq_dir}/SRX2333004*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' ${fastq_dir}/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' ${fastq_dir}/*
