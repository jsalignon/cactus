

## NOTE: This script should be run after having run the script to create the 
#        test datasets (as the files are actually downloaded in this script)


################
# Initialization

# setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
create_test_ds_dir=$cactus_dir/scripts/create_test_datasets
create_test_datasets_bin_dir=$create_test_ds_dir/bin
samples_ids_dir=$create_test_ds_dir/samples_ids
test_datasets_dir=$cactus_dir/test_datasets
singularity_dir=$homedir/workspace/singularity_containers/
# case_study_dir=$cactus_dir/case_study
case_study_dir=$test_datasets_dir/application_note

# loading functions
source $create_test_datasets_bin_dir/create_test_datasets_functions.sh

# changing directory
cd $case_study_dir

# Initialization
################



#######
# human

species=human
case_study_specie_dir=$case_study_dir/$species

cp -r $species/design $case_study_specie_dir
sed -i 's/data.*K_reads_atac/..\/..\/preprocessing\/human\/fastq\/atac/' $case_study_specie_dir/design/atac_fastq.tsv
sed -i 's/data.*K_reads_mrna/..\/..\/preprocessing\/human\/fastq\/mrna/' $case_study_specie_dir/design/mrna_fastq.tsv

# run.yml
yml_file="$case_study_specie_dir/parameters/run.yml"
cp ${species}/parameters/run.yml $yml_file
sed -i 's/rank/FDR/g' $yml_file
sed -i 's/200, 1000/1.3, 3/g' $yml_file

# human
#######




######
# worm

species=worm

# run.yml
mkdir -p $case_study_specie_dir/parameters 
yml_file="$case_study_specie_dir/parameters/run.yml"
cp ${species}/parameters/run.yml $yml_file
sed -i 's/rank/FDR/g' $yml_file
sed -i 's/200, 1000/1.3, 3/g' $yml_file

# tsv files
cp -r $species/design $case_study_specie_dir
sed -i 's/data.*K_reads_atac/..\/..\/preprocessing\/worm\/fastq\/atac/' $case_study_specie_dir/design/atac_fastq.tsv
sed -i 's/data.*K_reads_mrna/..\/..\/preprocessing\/worm\/fastq\/mrna/' $case_study_specie_dir/design/mrna_fastq.tsv

cat $case_study_specie_dir/design/atac_fastq.tsv


## => more analysis are need to add the hmg-3 samples that are not included
#     in the test dataset

design_dir="$case_study_specie_dir/design/"
cat $design_dir/atac_fastq.tsv
cat $design_dir/mrna_fastq.tsv

cd $create_test_datasets_bin_dir/../samples_ids
gsm_to_srr worm_case_study
cd $test_datasets_dir

data_dir=$case_study_specie_dir/data
mkdir -p $data_dir/mrna $data_dir/atac
fastq_dir=${prepro_dir}/fastq
if [ -d $fastq_dir ]; then rm -r $fastq_dir ; fi

prepro_dir="preprocessing/worm_case_study"
nextflow run nf-core/fetchngs --input "$samples_ids_dir/srr_accession/srr_worm_case_study.txt" --outdir $prepro_dir -profile singularity -r 1.6 -resume
# checking sample details to rename them
make_samples_info_file ${prepro_dir}

# renaming files
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' $fastq_dir/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' $fastq_dir/*
ls $fastq_dir

cat ${prepro_dir}/samplesheet/samples_info.tsv 

rename -v 's/SRX/mrna_SRX/' $fastq_dir/SRX30291{09..11}*
rename -v 's/SRX/atac_SRX/' $fastq_dir/SRX30291{27..29}*

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

fastq_info_file="${prepro_dir}/samplesheet/fastq_info.tsv"
awk 'BEGIN {OFS=""} { \
  library_strategy = tolower($5) ; \
  gsub(/-seq/, "", library_strategy) ; \
  gsub(/rna/, "mrna", library_strategy) ; \
  library_layout = $4
  gsub(/SINGLE/, "", library_layout) ; \
  gsub(/PAIRED/, "_R1", library_layout) ; \
  if (NR != 1) print $7, " ../../preprocessing/worm_case_study/fastq/", library_strategy, "_", $1, "_", $2, library_layout, ".fastq.gz" \
}' ${prepro_dir}/samplesheet/samples_info_1.tsv > ${fastq_info_file}
cat ${fastq_info_file}

grep atac ${fastq_info_file} >> $design_dir/atac_fastq.tsv
grep mrna ${fastq_info_file} >> $design_dir/mrna_fastq.tsv
awk -i inplace -v OFS="\t" '$1=$1' $design_dir/atac_fastq.tsv
awk -i inplace -v OFS="\t" '$1=$1' $design_dir/mrna_fastq.tsv

# comparisons.tsv
cat > ${design_dir}/comparisons.tsv <<EOL
hmg3 ctl
hmg4 ctl
spt16 ctl
hmg4 spt16
hmg3 spt16
hmg3 hmg4
EOL

# groups.tsv
cat > ${design_dir}/groups.tsv << EOL
all hmg3_vs_ctl hmg4_vs_ctl spt16_vs_ctl hmg4_vs_spt16 hmg3_vs_spt16 hmg3_vs_hmg4
ctl hmg3_vs_ctl hmg4_vs_ctl spt16_vs_ctl
fact hmg4_vs_spt16 hmg3_vs_spt16 hmg3_vs_hmg4
EOL

replace_spaces_by_tabs_in_the_design_tsv_files $case_study_specie_dir

cat $design_dir/atac_fastq.tsv
cat $design_dir/mrna_fastq.tsv
cat $design_dir/groups.tsv
cat $design_dir/comparisons.tsv

# worm
######
