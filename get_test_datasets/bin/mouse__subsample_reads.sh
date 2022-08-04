#!/bin/bash

specie="mouse"
prepro_dir="preprocessing/${specie}"
source $get_test_datasets_dir/get_test_datasets_functions.sh

n_reads_atac=$1
n_reads_mrna=$2

# subsampling reads
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_atac --experiment atac
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_mrna --experiment mrna

# adding the sample_id column (to edit manually)
awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(mrna_|atac_|rep)/, "", sample_id); \
  gsub(/old_/, "Old", sample_id); \
  gsub(/young_/, "Yng", sample_id); \
  gsub(/kidney/, "Kid", sample_id); \
  gsub(/liver/, "Liv", sample_id); \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv
cat ${prepro_dir}/samplesheet/samples_info_1.tsv

# making the fastq design files
make_fastq_info_file $specie $n_reads_atac $n_reads_mrna

# making the comparison design file
cat > ${specie}/design/comparisons.tsv <<EOL
YngKid OldKid
YngLiv OldLiv
YngKid YngLiv
OldKid OldLiv
EOL

# making the groups design file
cat > ${specie}/design/groups.tsv << EOL
all YngKid_vs_OldKid YngLiv_vs_OldLiv YngKid_vs_YngLiv OldKid_vs_OldLiv
age YngKid_vs_OldKid YngLiv_vs_OldLiv 
tissue YngKid_vs_YngLiv OldKid_vs_OldLiv
EOL

# regions to remove
touch ${specie}/design/regions_to_remove.tsv

