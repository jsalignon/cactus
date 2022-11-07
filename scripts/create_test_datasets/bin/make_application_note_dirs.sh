

############################################################
species=worm

cp -r $species/design application_note/$species
sed -i 's/data.*K_reads_atac/..\/..\/preprocessing\/worm\/fastq\/atac/' application_note/$species/design/atac_fastq.tsv
sed -i 's/data.*K_reads_mrna/..\/..\/preprocessing\/worm\/fastq\/mrna/' application_note/$species/design/mrna_fastq.tsv

# run.yml
mkdir -p application_note/$species/parameters 
yml_file="application_note/$species/parameters/run.yml"
cp ${species}/parameters/run.yml $yml_file
sed -i 's/rank/FDR/g' $yml_file
sed -i 's/200, 1000/1.3, 3/g' $yml_file



############################################################
species=human

cp -r $species/design application_note/$species
sed -i 's/data.*K_reads_atac/..\/..\/preprocessing\/human\/fastq\/atac/' application_note/$species/design/atac_fastq.tsv
sed -i 's/data.*K_reads_mrna/..\/..\/preprocessing\/human\/fastq\/mrna/' application_note/$species/design/mrna_fastq.tsv

# run.yml
yml_file="application_note/$species/parameters/run.yml"
cp ${species}/parameters/run.yml $yml_file
sed -i 's/rank/FDR/g' $yml_file
sed -i 's/200, 1000/1.3, 3/g' $yml_file


