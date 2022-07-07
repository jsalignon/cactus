

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
depot_galaxy = "https://depot.galaxyproject.org/singularity/"
singularity pull "$depot_galaxy/cvbio:3.0.0--hdfd78af_1"


##############################################
### Genome Fasta and Annotation files
##############################################

cd $data_dir 

release="106"
sed -i 's/\r//g' $get_data_dir/bin/get_fasta_and_gff.sh # this command is needed to replace windows carriage return to unix
source $get_data_dir/bin/get_fasta_and_gff.sh

get_fasta_and_gff $release dna_sm caenorhabditis_elegans  WBcel235 toplevel         worm
get_fasta_and_gff $release dna_sm drosophila_melanogaster BDGP6.28 toplevel         fly
get_fasta_and_gff $release dna_sm mus_musculus            GRCm38   primary_assembly mouse
get_fasta_and_gff $release dna_sm homo_sapiens            GRCh38   primary_assembly human


nextflow $get_data_dir/parse_genome_files.nf 



# TO DO: add the code below to the pipeline


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
prepro_dir="preprocessing"
prepro_bl_dir="$prepro_dir/blacklisted_regions"
mkdir -p $prepro_bl_dir

sed -i 's/\r//g' $get_data_dir/bin/get_blacklisted_ensembl_file.sh # this command is needed to replace windows carriage return to unix
source $get_data_dir/bin/get_blacklisted_ensembl_file.sh

get_blacklisted_ensembl_file worm  ce11 WBcel235 $prepro_bl_dir
get_blacklisted_ensembl_file fly   dm6  BDGP6    $prepro_bl_dir
get_blacklisted_ensembl_file mouse mm10 GRCm38   $prepro_bl_dir
get_blacklisted_ensembl_file human hg38 GRCh38   $prepro_bl_dir



depot_galaxy = "https://depot.galaxyproject.org/singularity/"
params.annotationhub = "${depot_galaxy}/bioconductor-annotationhub:3.2.0--r41hdfd78af_0"


singularity shell bioconductor-annotationhub:3.2.0--r41hdfd78af_0

R

library(AnnotationHub)

ah = AnnotationHub()

specie_initial = 'Ce'
orgdb_sqlite_file = paste0('org.', specie_initial, '.eg.db.sqlite')
mcols(query(ah, 'OrgDb', '2021-10-08'))[1,] %>% as.data.frame
mcols(query(ah, 'ncbi/standard/3.14/org.Ce.eg.db.sqlite'))[1,] %>% as.data.frame
orgdb = query(ah, 'ncbi/standard/3.14/org.Ce.eg.db.sqlite')[[1]]













