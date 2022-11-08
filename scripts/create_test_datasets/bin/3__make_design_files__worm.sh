#!/bin/bash

species="worm"
n_reads_atac=$1
n_reads_mrna=$2
prepro_dir="preprocessing/${species}"
source $create_test_datasets_bin_dir/create_test_datasets_functions.sh


# run.yml
cat > ${species}/parameters/run.yml << EOL
res_dir                   : 'results/test_worm'
species                    : 'worm'
chromatin_state           : 'iHMM.M1K16.worm_L3'
split__threshold_type     : 'rank' 
split__threshold_values   : [ 200, 1000 ]
EOL

# atac_fastq.tsv and mrna_fastq.tsv
make_samples_info_file ${prepro_dir}
awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(rep|-|_control)/, "", sample_id); \
  gsub(/rluc/, "ctl", sample_id); \
  gsub(/ctrl/, "ctl", sample_id); \
  gsub(/sample_title_sample_title/, "sample_id", sample_id); \
  if (NR == 1) sample_id = "sample_id"; \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv
make_fastq_info_file $species $n_reads_atac $n_reads_mrna
cp ${species}/design/atac_fastq.tsv ${species}/design/atac_fastq__with_input.tsv
grep -v input ${species}/design/atac_fastq__with_input.tsv > ${species}/design/atac_fastq__without_input.tsv

# comparisons.tsv
cat > ${species}/design/comparisons.tsv <<EOL
hmg4 ctl
spt16 ctl
hmg4 spt16
EOL

# groups.tsv
cat > ${species}/design/groups.tsv << EOL
all hmg4_vs_ctl spt16_vs_ctl hmg4_vs_spt16
ctl hmg4_vs_ctl spt16_vs_ctl
spt16 spt16_vs_ctl hmg4_vs_spt16
EOL

# regions_to_remove.tsv
cat > ${species}/design/regions_to_remove.tsv << EOL
hmg4 Hmg4->chrIII:7,379,143-7,381,596
spt16 Spt16->chrI:10,789,130-10,793,152
EOL

# genes_to_remove.tsv
touch ${species}/design/genes_to_remove.tsv


# regions_to_remove_empty.tsv
touch ${species}/design/regions_to_remove_empty.tsv

# run__no_rtr.yml
yml_file="${species}/parameters/run__no_rtr.yml"
cp ${species}/parameters/run.yml $yml_file
sed -i 's/test_worm/test_worm__no_rtr/g' $yml_file
cat >> $yml_file << EOL
design__regions_to_remove : 'design/regions_to_remove_empty.tsv'
disable_all_enrichments   : true
EOL

# run__enrich_only_genes_self.yml
yml_file="${species}/parameters/run__enrich_only_genes_self.yml"
cp ${species}/parameters/run.yml $yml_file
sed -i 's/test_worm/test_worm__enrich_only_genes_self/g' $yml_file
sed -i 's/\[ 200, 1000 \]/\[ 200 \]/g' $yml_file
cat >> $yml_file << EOL
do_genes_self_enrichment  : true
do_peaks_self_enrichment  : false
do_gene_set_enrichment    : false
do_chrom_state_enrichment : false
do_motif_enrichment       : false
do_chip_enrichment        : false
EOL

replace_spaces_by_tabs_in_the_design_tsv_files $species
