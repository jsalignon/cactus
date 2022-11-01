

############################################################
specie=worm

cp -r $specie/design application_note/$specie
sed -i 's/data.*K_reads_atac/..\/..\/preprocessing\/worm\/fastq\/atac/' application_note/$specie/design/atac_fastq.tsv
sed -i 's/data.*K_reads_mrna/..\/..\/preprocessing\/worm\/fastq\/mrna/' application_note/$specie/design/mrna_fastq.tsv

# run.yml
mkdir -p application_note/$specie/parameters 
yml_file="application_note/$specie/parameters/run.yml"
cp ${specie}/parameters/run.yml $yml_file
sed -i 's/rank/FDR/g' $yml_file
sed -i 's/200, 1000/1.3, 3/g' $yml_file



############################################################
specie=human

cp -r $specie/design application_note/$specie
sed -i 's/data.*K_reads_atac/..\/..\/preprocessing\/human\/fastq\/atac/' application_note/$specie/design/atac_fastq.tsv
sed -i 's/data.*K_reads_mrna/..\/..\/preprocessing\/human\/fastq\/mrna/' application_note/$specie/design/mrna_fastq.tsv

# run.yml
yml_file="application_note/$specie/parameters/run.yml"
cp ${specie}/parameters/run.yml $yml_file
sed -i 's/rank/FDR/g' $yml_file
sed -i 's/200, 1000/1.3, 3/g' $yml_file


