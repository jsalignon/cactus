

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





























## Chromatin states: to do later... pretty complicated...

# 
# process getting_chromatin_states {
#   tag "${specie}"
# 
#   container = params.liftover
# 
#   publishDir path: "${specie}/chromatin_state", mode: 'link'
# 
#   input:
# 		set specie, file(gff3_genes_only), file(gff3_without_pseudogenes_and_ncRNA) from channel_species_2
# 
#   output:
#     file('*')
# 
#   shell:
#   '''
#       url_liftover="https://hgdownload.cse.ucsc.edu/goldenPath/ce10/liftOver"
#       wget $url_liftover/ce10ToCe11.over.chain.gz
# 
# 		#!/usr/bin/env Rscript
# 
# 		library(AnnotationHub)
# 		specie = '!{specie}'
# 
# 		specie_initial = switch(specie, 
# 			worm  = 'Ce', 
# 			fly   = 'Dm', 
# 			mouse = 'Mm', 
# 			human = 'Hs'
# 		)
# 
# 		orgdb_sqlite_file = paste0('org.', specie_initial, '.eg.db.sqlite')
# 		orgdb <- query(ah, c("OrgDb", "maintainer@bioconductor.org"))[[orgdb_sqlite_file]]
# 
# 		orgdb = AnnotationDbi::loadDb(paste0(genomedir, '/../org_db/org_db_', specie_initial, '.sqlite'))
# 
# 		orgdb_sqlite_file = paste0('org.', specie_initial, '.eg.db.sqlite')
# 		// mcols(query(ah, 'OrgDb', '2021-10-08'))[1,] %>% as.data.frame
# 		// mcols(query(ah, 'ncbi/standard/3.14/org.Ce.eg.db.sqlite'))[1,] %>% as.data.frame
# 		// orgdb = query(ah, 'ncbi/standard/3.14/org.Ce.eg.db.sqlite')[[1]]
# 
#   '''
# 
# }
# 
# 
# 
# # 
# # cd $LOP 
# # wget https://hgdownload.cse.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz
# # wget https://hgdownload.cse.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDm6.over.chain.gz
# # wget https://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz
# # wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
# 
# 
# ## initialization
# cd $OUTF ; mkdir chromatin_states
# cd $OUTF/chromatin_states ; mkdir raw modif lifted split
# 
# 
# ## downloading raw chromatin state files 
# cd $OUTF/chromatin_states/raw
# 
# for BED in iHMM.M1K16.fly_EL.bed iHMM.M1K16.fly_L3.bed iHMM.M1K16.human_GM.bed iHMM.M1K16.human_H1.bed iHMM.M1K16.worm_EE.bed iHMM.M1K16.worm_L3.bed
#  do wget http://compbio.med.harvard.edu/modencode/webpage/hihmm/$BED
# done
# 
# ## modifying names and custom files
# cp iHMM.M1K16.worm_EE.bed ../modif/worm_EE_ce10.bed 
# cp iHMM.M1K16.worm_L3.bed ../modif/worm_L3_ce10.bed 
# cp iHMM.M1K16.fly_EL.bed ../modif/fly_EL_dm3.bed 
# cp iHMM.M1K16.fly_L3.bed ../modif/fly_L3_dm3.bed 
# cp iHMM.M1K16.human_GM.bed ../modif/human_GM12878_hg19.bed 
# cp iHMM.M1K16.human_H1.bed ../modif/human_H1_hESC_hg19.bed 
# 
# gawk -i inplace '{print "chr"$0}' ../modif/worm_EE_ce10.bed 
# gawk -i inplace '{print "chr"$0}' ../modif/worm_L3_ce10.bed 
# 
# 
# ## remaping genomic coordinates to a newer genome version
# cd $OUTF/chromatin_states/lifted
# mkdir unmapped
# 
# $LOP/../liftOver ../modif/worm_EE_ce10.bed  $LOP/ce10ToCe11.over.chain.gz worm_EE_ce11.bed unmapped/worm_EE_ce11_unmapped.bed
# $LOP/../liftOver ../modif/worm_L3_ce10.bed  $LOP/ce10ToCe11.over.chain.gz worm_L3_ce11.bed unmapped/worm_L3_ce11_unmapped.bed
# 
# $LOP/../liftOver ../modif/fly_EL_dm3.bed  $LOP/dm3ToDm6.over.chain.gz fly_EL_dm6.bed unmapped/fly_EL_dm6_unmapped.bed 
# $LOP/../liftOver ../modif/fly_L3_dm3.bed  $LOP/dm3ToDm6.over.chain.gz fly_L3_dm6.bed unmapped/fly_L3_dm6_unmapped.bed 
# 
# $LOP/../liftOver ../modif/human_GM12878_hg19.bed  $LOP/hg19ToHg38.over.chain.gz human_GM12878_hg38.bed unmapped/human_GM12878_hg38_unmapped.bed 
# $LOP/../liftOver ../modif/human_H1_hESC_hg19.bed  $LOP/hg19ToHg38.over.chain.gz human_H1_hESC_hg38.bed unmapped/human_H1_hESC_hg19_unmapped.bed 
# 
# 
# ## splitting bed files with one file for each chromatin state
# cd $OUTF/chromatin_states/split
# 
# cp ../lifted/*.bed .
# awk ' { FOLDERN=FILENAME; gsub(".bed", "", FOLDERN); "mkdir -p " FOLDERN | getline ; print >  FOLDERN"/"$4".bed" }' *.bed 
# rm *.bed
# 
# # removing the unmapped files
# find . -name "17_Unmap.bed" -delete









