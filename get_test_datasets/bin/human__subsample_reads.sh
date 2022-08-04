#!/bin/bash

specie="human"
prepro_dir="preprocessing/${specie}"
source $get_test_datasets_dir/get_test_datasets_functions.sh

n_reads_atac=$1
n_reads_mrna=$2

# subsampling reads
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_atac --experiment atac
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_mrna --experiment mrna

# adding the sample_id column (to edit manually)
awk 'BEGIN {OFS = "\t"} { \
  condition = tolower($6) ; \
  replicate = $6 ; \
  gsub(/.*rep/, "", replicate); \
  gsub(/.*mock.*/, "ctl", condition); \
  gsub(/.*ssrp1.*/, "ssrp1", condition); \
  gsub(/.*supt16h.*/, "supt16h", condition); \
  sample_id = condition "_" replicate ; \
  gsub(/sample_title_sample_title/, "sample_id", sample_id)
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv
cat ${prepro_dir}/samplesheet/samples_info_1.tsv

# making the fastq design files
make_fastq_info_file $specie $n_reads_atac $n_reads_mrna

# making the comparison design file
cat > ${specie}/design/comparisons.tsv <<EOL
ssrp1 ctl
supt16h ctl
ssrp1 supt16h
EOL

# making the groups design file
cat > ${specie}/design/groups.tsv << EOL
all ssrp1_vs_ctl supt16h_vs_ctl ssrp1_vs_supt16h
ctl ssrp1_vs_ctl supt16h_vs_ctl
supt16h supt16h_vs_ctl ssrp1_vs_supt16h
EOL

# regions to remove
cat > ${specie}/design/regions_to_remove.tsv << EOL
ssrp1 ssrp1->chr11:57,325,986-57,335,892
supt16h supt16h->chr14:21,351,476-21,384,019
EOL
