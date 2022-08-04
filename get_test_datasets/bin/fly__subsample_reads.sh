#!/bin/bash

specie="fly"
prepro_dir="preprocessing/${specie}"

n_reads_atac=$1
n_reads_mrna=$2

# subsampling reads
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_atac --experiment atac
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_mrna --experiment mrna

# adding the sample_id column (to edit manually)
awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(atac|rna)seq_rep/, "", sample_id); \
  gsub(/bap170/, "b170", sample_id); \
  gsub(/nurf301/, "n301", sample_id); \
  gsub(/lacz/, "ctl", sample_id); \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv

# making the fastq design files
make_fastq_info_file $specie $n_reads_atac $n_reads_mrna

# making the comparison design file
cat > ${specie}/design/comparisons.tsv <<EOL
gaf ctl
b170 ctl
n301 ctl
n301b170 ctl
b170 n301b170
n301 n301b170
EOL

# making the groups design file
cat > ${specie}/design/groups.tsv << EOL
all gaf_vs_ctl b170_vs_ctl n301_vs_ctl n301b170_vs_ctl b170_vs_n301b170 n301_vs_n301b170
ctl gaf_vs_ctl b170_vs_ctl n301_vs_ctl n301b170_vs_ctl
n301b170 n301b170_vs_ctl b170_vs_n301b170 n301_vs_n301b170
EOL

# regions to remove
cat > ${specie}/design/regions_to_remove.tsv << EOL
gaf gaf->3L:14,747,929-14,761,049
b170 bap170->2R:6,636,512-6,642,358
n301 nurf301->3L:233,926-246,912
n301b170 bap170->2R:6,636,512-6,642,358
n301b170 nurf301->3L:233,926-246,912
EOL
