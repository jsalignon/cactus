

##############################################
### Converting GSM ids to SRR ids
##############################################

# this was done since there was an issue with SRA that changed its data downloading procedures making the fetchngs pipeline not working anymore when downloading via SRA-tools. 
# a workaround that consist in using ENA instead of SRA is described here: https://github.com/nf-core/fetchngs/issues/98
# however this works only with SRA ids not GEO ids (because there is only a function to convert SRA its to ENA ids not from GEO ids in fetchngs) therefore we need to convert our GEO ids.
# this conversion was done using the entrez-direct API via a BioContainer

singularity_dir="/home/jersal/workspace/singularity_containers"

cd $singularity_dir
singularity pull https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1

gsm_to_srr (){
  specie=$1
  my_gsm_ids="$(cat samples_ids/gsm_accession/gsm_$specie.txt | awk '{print}' ORS=' OR ' RS='\r\n' )"
  singularity exec $singularity_dir/entrez-direct:16.2--he881be0_1 esearch -db sra -query "$my_gsm_ids" | singularity exec $singularity_dir/entrez-direct:16.2--he881be0_1 efetch -format runinfo | cut -d ',' -f 1 - | grep -v 'Run' - > "samples_ids/srr_accession/srr_$specie.txt"
}

cd /home/jersal/workspace/cactus_test_datasets
gsm_to_srr worm
gsm_to_srr human
gsm_to_srr mice
gsm_to_srr fly


# $1 = ${prepro_dir}   => head -3 ${prepro_dir}/samplesheet/samplesheet.csv 
make_samples_info_file (){
  cut -d"," -f5,6,15,17,20,28 ${1}/samplesheet/samplesheet.csv | sed 's/\"//g' | sed 's/ /_/g' | awk 'BEGIN { FS = ","; OFS = "\t"} {print $2, $1, $3, $4, $5, $6}'| column -t > ${1}/samplesheet/samples_info.tsv
}

# $1 = ${prepro_dir}   => cat ${prepro_dir}/samplesheet/samples_info.csv 
make_fastq_info_file (){
  awk 'BEGIN {OFS=""} { \
    library_strategy = tolower($5) ; \
    gsub(/-seq/, "", library_strategy) ; \
    gsub(/rna/, "mrna", library_strategy) ; \
    library_layout = $4
    gsub(/SINGLE/, "", library_layout) ; \
    gsub(/PAIRED/, "_R1", library_layout) ; \
    if (NR != 1) print $7, " data/", library_strategy, "/sample_100K_reads_", library_strategy, "_", $1, "_", $2, library_layout, ".fastq.gz" \
  }' ${1}/samplesheet/samples_info_1.tsv > ${1}/samplesheet/fastq_info.tsv
}



##############################################
### Human (GSE98758)
##############################################

specie="human"
prepro_dir="preprocessing/${specie}"
specie_dir="specie/${specie}"

# downloading our fastq samples of interest and subsampling them
nextflow run nf-core/fetchngs --input "samples_ids/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 --force_sratools_download 

# creating the sample_info file
make_samples_info_file ${prepro_dir}
cat ${prepro_dir}/samplesheet/samples_info.tsv

# renaming files
rename -v 's/SRX2794/atac_SRX2794/' ${prepro_dir}/fastq/*
rename -v 's/SRX4029/mrna_SRX4029/' ${prepro_dir}/fastq/*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' ${prepro_dir}/fastq/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' ${prepro_dir}/fastq/*

# subsampling reads
nextflow subsample_fastq.nf --indir ${prepro_dir}/fastq -resume

# moving data to appropriate folders
mv ${prepro_dir}/fastq_100K_reads/*_atac_*.fastq.gz ${specie_dir}/data/atac/
mv ${prepro_dir}/fastq_100K_reads/*_mrna_*.fastq.gz ${specie_dir}/data/mrna/

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

# making the fastq_info file
make_fastq_info_file ${prepro_dir}
cat ${prepro_dir}/samplesheet/fastq_info.tsv

# making the fastq design files
grep atac ${prepro_dir}/samplesheet/fastq_info.tsv > ${specie_dir}/design/atac_fastq.tsv
grep mrna ${prepro_dir}/samplesheet/fastq_info.tsv > ${specie_dir}/design/mrna_fastq.tsv
truncate -s -1 ${specie_dir}/design/atac_fastq.tsv
truncate -s -1 ${specie_dir}/design/mrna_fastq.tsv
cat ${specie_dir}/design/atac_fastq.tsv
cat ${specie_dir}/design/mrna_fastq.tsv

# making the comparison design file
cat > ${specie_dir}/design/comparisons.tsv <<EOL
ssrp1 ctl
supt16h ctl
ssrp1 supt16h
EOL
truncate -s -1 ${specie_dir}/design/comparisons.tsv

# making the groups design file
cat > ${specie_dir}/design/groups.tsv << EOL
all ssrp1_vs_ctl supt16h_vs_ctl ssrp1_vs_supt16h
ctl ssrp1_vs_ctl supt16h_vs_ctl
supt16h supt16h_vs_ctl ssrp1_vs_supt16h
EOL
truncate -s -1 ${specie_dir}/design/groups.tsv

# regions to remove
cat > ${specie_dir}/design/regions_to_remove.tsv << EOL
ssrp1 ssrp1->chr11:57,325,986-57,335,892
supt16h supt16h->chr14:21,351,476-21,384,019
EOL
truncate -s -1 ${specie_dir}/design/regions_to_remove.tsv


##############################################
### Worm (GSE98758)
##############################################

specie="worm"
prepro_dir="preprocessing/${specie}"
specie_dir="specie/${specie}"

# downloading our fastq samples of interest and subsampling them
# nextflow run nf-core/fetchngs --input samples_id/sra_accession/sra_acc_worm.txt --outdir preprocessing/worm -profile singularity -r 1.6 -resume
# nextflow run nf-core/fetchngs --input samples_ids/gsm_accession/gsm_worm.txt --outdir preprocessing/worm -profile singularity -r 1.6  --force_sratools_download 
nextflow run nf-core/fetchngs --input "samples_ids/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 -resume

# creating a simple reference file
make_samples_info_file ${prepro_dir}
cat ${prepro_dir}/samplesheet/samples_info.tsv

# renaming files
rename -v 's/SRX/mrna_SRX/' ${prepro_dir}/fastq/SRX30291{12..20}*
rename -v 's/SRX/atac_SRX/' ${prepro_dir}/fastq/SRX30291{24..35}*
rename -v 's/SRX/atac_SRX/' ${prepro_dir}/fastq/SRX2333004*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' ${prepro_dir}/fastq/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' ${prepro_dir}/fastq/*

# subsampling reads
nextflow subsample_fastq.nf --indir ${prepro_dir}/fastq -resume

# moving data to appropriate folders
mv ${prepro_dir}/fastq_100K_reads/*_atac_*.fastq.gz ${specie_dir}/data/atac/
mv ${prepro_dir}/fastq_100K_reads/*_mrna_*.fastq.gz ${specie_dir}/data/mrna/

# adding the sample_id column (to edit manually)
awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(rep|-|_control)/, "", sample_id); \
  gsub(/rluc/, "ctl", sample_id); \
  gsub(/ctrl/, "ctl", sample_id); \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv

# making the fastq_info file
make_fastq_info_file ${prepro_dir}
cat ${prepro_dir}/samplesheet/fastq_info.tsv

# making the fastq design files
grep atac ${prepro_dir}/samplesheet/fastq_info.tsv > ${specie_dir}/design/atac_fastq.tsv
grep mrna ${prepro_dir}/samplesheet/fastq_info.tsv > ${specie_dir}/design/mrna_fastq.tsv

cp ${specie_dir}/design/atac_fastq.tsv ${specie_dir}/design/atac_fastq__with_input.tsv
grep -v input ${specie_dir}/design/atac_fastq__with_input.tsv > ${specie_dir}/design/atac_fastq__without_input.tsv

# making the comparison design file
cat > ${specie_dir}/design/comparisons.tsv <<EOL
hmg4 ctl
spt16 ctl
hmg4 spt16
EOL

# making the groups design file
cat > ${specie_dir}/design/groups.tsv << EOL
all hmg4_vs_ctl spt16_vs_ctl hmg4_vs_spt16
ctl hmg4_vs_ctl spt16_vs_ctl
spt16 spt16_vs_ctl hmg4_vs_spt16
EOL

# regions to remove
cat > ${specie_dir}/design/regions_to_remove.tsv << EOL
hmg4 Hmg4->chrIII:7,379,143-7,381,596
spt16 Spt16->chrI:10,789,130-10,793,152
EOL


##############################################
### Mice (GSE181797)
##############################################

specie="mice"
prepro_dir="preprocessing/${specie}"
specie_dir="specie/${specie}"

# downloading our fastq samples of interest and subsampling them
nextflow run nf-core/fetchngs --input "samples_ids/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 -resume

# creating a simple reference file
make_samples_info_file ${prepro_dir}
cat ${prepro_dir}/samplesheet/samples_info.tsv

# renaming files
cat ${prepro_dir}/samplesheet/samples_info.csv
rename -v 's/SRX/mrna_SRX/' ${prepro_dir}/fastq/SRX117086{63..78}*
rename -v 's/SRX/atac_SRX/' ${prepro_dir}/fastq/SRX117086{79..90}*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' ${prepro_dir}/fastq/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' ${prepro_dir}/fastq/*
ls ${prepro_dir}/fastq

# subsampling reads
nextflow subsample_fastq.nf --indir ${prepro_dir}/fastq -resume

# moving data to appropriate folders
mv ${prepro_dir}/fastq_100K_reads/*_atac_*.fastq.gz ${specie_dir}/data/atac/
mv ${prepro_dir}/fastq_100K_reads/*_mrna_*.fastq.gz ${specie_dir}/data/mrna/

# adding the sample_id column (to edit manually)
awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(mrna_|atac_|rep)/, "", sample_id); \
  gsub(/young/, "yng", sample_id); \
  gsub(/kidney/, "kid", sample_id); \
  gsub(/liver/, "liv", sample_id); \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv

# making the fastq_info file
make_fastq_info_file ${prepro_dir}
cat ${prepro_dir}/samplesheet/fastq_info.tsv

# making the fastq design files
grep atac ${prepro_dir}/samplesheet/fastq_info.tsv > ${specie_dir}/design/atac_fastq.tsv
grep mrna ${prepro_dir}/samplesheet/fastq_info.tsv > ${specie_dir}/design/mrna_fastq.tsv
truncate -s -1 ${specie_dir}/design/atac_fastq.tsv
truncate -s -1 ${specie_dir}/design/mrna_fastq.tsv

# making the comparison design file
cat > ${specie_dir}/design/comparisons.tsv <<EOL
yng_kid old_kid
yng_liv old_liv
yng_kid yng_liv
old_kid old_liv
EOL
truncate -s -1 ${specie_dir}/design/comparisons.tsv

# making the groups design file
cat > ${specie_dir}/design/groups.tsv << EOL
all yng_kid_vs_old_kid yng_liv_vs_old_liv yng_kid_vs_yng_liv old_kid_vs_old_liv
age yng_kid_vs_old_kid yng_liv_vs_old_liv 
tissue yng_kid_vs_yng_liv old_kid_vs_old_liv
EOL
truncate -s -1 ${specie_dir}/design/groups.tsv

# regions to remove
touch ${specie_dir}/design/regions_to_remove.tsv


##############################################
### Fly (GSE149339)
##############################################

specie="fly"
prepro_dir="preprocessing/${specie}"
specie_dir="specie/${specie}"

# downloading our fastq samples of interest and subsampling them
nextflow run nf-core/fetchngs --input "samples_ids/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 -resume

# creating a simple reference file
make_samples_info_file ${prepro_dir}
cat ${prepro_dir}/samplesheet/samples_info.tsv

# renaming files
cat ${prepro_dir}/samplesheet/samples_info.csv
rename -v 's/SRX/mrna_SRX/' ${prepro_dir}/fastq/SRX81740{44..53}*
rename -v 's/SRX/atac_SRX/' ${prepro_dir}/fastq/SRX81740{34..43}*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' ${prepro_dir}/fastq/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' ${prepro_dir}/fastq/*
ls ${prepro_dir}/fastq

# subsampling reads
nextflow subsample_fastq.nf --indir ${prepro_dir}/fastq -resume

# moving data to appropriate folders
mv ${prepro_dir}/fastq_100K_reads/*_atac_*.fastq.gz ${specie_dir}/data/atac/
mv ${prepro_dir}/fastq_100K_reads/*_mrna_*.fastq.gz ${specie_dir}/data/mrna/

# adding the sample_id column (to edit manually)
awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(atac|rna)seq_rep/, "", sample_id); \
  gsub(/bap170/, "b170", sample_id); \
  gsub(/nurf301/, "n301", sample_id); \
  gsub(/lacz/, "ctl", sample_id); \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv

# making the fastq_info file
make_fastq_info_file ${prepro_dir}
cat ${prepro_dir}/samplesheet/fastq_info.tsv

# making the fastq design files
grep atac ${prepro_dir}/samplesheet/fastq_info.tsv > ${specie_dir}/design/atac_fastq.tsv
grep mrna ${prepro_dir}/samplesheet/fastq_info.tsv > ${specie_dir}/design/mrna_fastq.tsv
truncate -s -1 ${specie_dir}/design/atac_fastq.tsv
truncate -s -1 ${specie_dir}/design/mrna_fastq.tsv

# making the comparison design file
cat > ${specie_dir}/design/comparisons.tsv <<EOL
gaf ctl
b170 ctl
n301 ctl
n301b170 ctl
b170 n301b170
n301 n301b170
EOL
truncate -s -1 ${specie_dir}/design/comparisons.tsv

# making the groups design file
cat > ${specie_dir}/design/groups.tsv << EOL
all gaf_vs_ctl b170_vs_ctl n301_vs_ctl n301b170_vs_ctl b170_vs_n301b170 n301_vs_n301b170
ctl gaf_vs_ctl b170_vs_ctl n301_vs_ctl n301b170_vs_ctl
n301b170 n301b170_vs_ctl b170_vs_n301b170 n301_vs_n301b170
EOL
truncate -s -1 ${specie_dir}/design/groups.tsv

# regions to remove
cat > ${specie_dir}/design/regions_to_remove.tsv << EOL
gaf gaf->3L:14,747,929-14,761,049
b170 bap170->2R:6,636,512-6,642,358
n301 nurf301->3L:233,926-246,912
n301b170 bap170->2R:6,636,512-6,642,358
n301b170 nurf301->3L:233,926-246,912
EOL
truncate -s -1 ${specie_dir}/design/regions_to_remove.tsv


##############################################
### Copying the datasets to another folder to test them 
##############################################

cp -r specie/* ../test_test_datasets/



##############################################
### The end
##############################################

# saving compressed objects and their sizes
du -h -d1 specie/ > test_datasets_sizes.txt
tar --use-compress-program="pigz -p 15 -k -r" -cf specie.tar.gz specie


