

cactus_dir=~/workspace/cactus
singularity_dir=~/workspace/singularity_containers

get_data_dir=$cactus_dir/software/get_data
data_dir=$cactus_dir/data

export NXF_SINGULARITY_CACHEDIR=${singularity_dir}

tree -d data_archive/worm/




##############################################
### Creating all directories and getting the needed containers
##############################################

cd $data_dir ; mkdir worm ; cd worm
mkdir -p blacklisted_regions bowtie2_indexes_conta chromatin_states encode_CHIP genome/annotation genome/sequence homer_genome motifs_PWMS org_db refGene_UCSC
cd $data_dir ; cp -r worm fly ; cp -r worm mouse ; cp -r worm human
tree -d


cd $singularity_dir
singularity pull https://depot.galaxyproject.org/singularity/cvbio:3.0.0--hdfd78af_1  


##############################################
### Genome Fasta and Annotation files
##############################################

cd $data_dir 

sed -i 's/\r//g' $get_data_dir/bin/get_fasta_and_gff.sh # this command is needed to replace windows carriage return to unix
source $get_data_dir/bin/get_fasta_and_gff.sh

get_fasta_and_gff 100 dna_sm caenorhabditis_elegans WBcel235 toplevel worm
get_fasta_and_gff 100 dna_sm drosophila_melanogaster BDGP6.28 toplevel fly
get_fasta_and_gff 100 dna_sm mus_musculus GRCm38 primary_assembly mouse
get_fasta_and_gff 100 dna_sm homo_sapiens GRCh38 primary_assembly human


# getting the R objects seq_info (size of chromosome)
: '
R
library(magrittr)
v_specie_code = c(human = 'hg38', mouse = 'mm10', fly = 'dm6', worm = 'ce11')
purrr::iwalk(v_specie_code, function(code, specie){
  print(specie)
  seqinfo = rtracklayer::SeqinfoForUCSCGsenome(code)
  seqinfo@seqnames %<>% gsub('chr', '', .)
  saveRDS(seqinfo, paste0(specie, '/genome/sequence/seqinfo.rds'))
})
'




##############################################
### Blacklisted regions
##############################################



cd $data_dir 
prepro_dir="preprocessing/blacklisted_regions"
mkdir -p $prepro_dir

singularity pull cvbio:3.0.0--hdfd78af_1 

use a container for cvbio
# cvbio:3.0.0--hdfd78af_1

get_blacklisted_ensembl_file (){
  SPECIE_CODE=$1
  NCBI_CODE=$2
  BLACKLIST_FILE="${SPECIE_CODE}-blacklist.v2.bed.gz"
  MAPPING_FILE="${NCBI_CODE}_UCSC2ensembl.txt"
  MAPPING_PATH="https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master"
  BLACKLIST_PATH="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists"
  wget -O blacklist_ncbi.bed.gz $BLACKLIST_PATH/$BLACKLIST_FILE
  gunzip blacklist_ncbi.bed.gz
  wget -O NCBI_to_Ensembl.txt $MAPPING_PATH/$MAPPING_FILE
  singularity run $singularity_dir/cvbio:3.0.0--hdfd78af_1 cvbio UpdateContigNames -i blacklist_ncbi.bed -o blacklist_ensembl.bed -m NCBI_to_Ensembl.txt --comment-chars '#' --columns 0 --skip-missing true
  mv blacklist_ncbi.bed > $prepro_dir/${SPECIE_CODE}_blacklist_ncbi.bed
  mv NCBI_to_Ensembl.bed > $prepro_dir/${SPECIE_CODE}_NCBI_to_Ensembl.bed
}



# test:
SPECIE_CODE="ce11"
NCBI_CODE="WBcel235"

for MAPPING_FILE in BDGP6_UCSC2ensembl.txt WBcel235_UCSC2ensembl.txt GRCm38_UCSC2ensembl.txt GRCh38_UCSC2ensembl.txt
 do wget https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/$MAPPING_FILE
done 

cd $OUTF ; mkdir blacklisted_regions
cd $OUTF/blacklisted_regions

for BLACKLIST_FILE in ce11-blacklist.v2.bed.gz dm6-blacklist.v2.bed.gz mm10-blacklist.v2.bed.gz hg38-blacklist.v2.bed.gz
 do wget https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/$BLACKLIST_FILE
done
gunzip *
mv ce11-blacklist.v2.bed ce11_blacklist_v2.bed
mv dm6-blacklist.v2.bed dm6_blacklist_v2.bed
mv mm10-blacklist.v2.bed mm10_blacklist_v2.bed
mv hg38-blacklist.v2.bed hg38_blacklist_v2.bed

cvbio UpdateContigNames -i ce11_blacklist_v2.bed -o ce11_blacklist_v2__Ensembl_named.bed -m $CMP/WBcel235_UCSC2ensembl.txt --comment-chars '#' --columns 0 --skip-missing true
cvbio UpdateContigNames -i dm6_blacklist_v2.bed -o dm6_blacklist_v2__Ensembl_named.bed -m $CMP/BDGP6_UCSC2ensembl.txt --comment-chars '#' --columns 0 --skip-missing true
cvbio UpdateContigNames -i mm10_blacklist_v2.bed -o mm10_blacklist_v2__Ensembl_named.bed -m $CMP/GRCm38_UCSC2ensembl.txt --comment-chars '#' --columns 0 --skip-missing true
cvbio UpdateContigNames -i hg38_blacklist_v2.bed -o hg38_blacklist_v2__Ensembl_named.bed -m $CMP/GRCh38_UCSC2ensembl.txt --comment-chars '#' --columns 0 --skip-missing true

wc -l *
