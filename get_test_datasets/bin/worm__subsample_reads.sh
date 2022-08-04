#!/bin/bash

specie="worm"
prepro_dir="preprocessing/${specie}"
source $get_test_datasets_bin_dir/get_test_datasets_functions.sh

n_reads_atac=$1
n_reads_mrna=$2

# subsampling reads
nextflow $get_test_datasets_bin_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_atac --experiment atac
nextflow $get_test_datasets_bin_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_mrna --experiment mrna

# adding the sample_id column (to edit manually)
awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(rep|-|_control)/, "", sample_id); \
  gsub(/rluc/, "ctl", sample_id); \
  gsub(/ctrl/, "ctl", sample_id); \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv

# making the fastq design files
make_fastq_info_file $specie $n_reads_atac $n_reads_mrna

cp ${specie}/design/atac_fastq.tsv ${specie}/design/atac_fastq__with_input.tsv
grep -v input ${specie}/design/atac_fastq__with_input.tsv > ${specie}/design/atac_fastq__without_input.tsv

# making the comparison design file
cat > ${specie}/design/comparisons.tsv <<EOL
hmg4 ctl
spt16 ctl
hmg4 spt16
EOL

# making the groups design file
cat > ${specie}/design/groups.tsv << EOL
all hmg4_vs_ctl spt16_vs_ctl hmg4_vs_spt16
ctl hmg4_vs_ctl spt16_vs_ctl
spt16 spt16_vs_ctl hmg4_vs_spt16
EOL

# regions to remove
cat > ${specie}/design/regions_to_remove.tsv << EOL
hmg4 Hmg4->chrIII:7,379,143-7,381,596
spt16 Spt16->chrI:10,789,130-10,793,152
EOL
