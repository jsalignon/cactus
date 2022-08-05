#!/bin/bash

specie="mouse"
prepro_dir="preprocessing/${specie}"
source $get_test_datasets_bin_dir/get_test_datasets_functions.sh

n_reads_atac=$1
n_reads_mrna=$2


# ATAC-Seq
nextflow $get_test_datasets_bin_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_atac --experiment atac

# mRNA-Seq
nextflow $get_test_datasets_bin_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_mrna --experiment mrna
