


##################################################################
## Initialization
##################################################################


ATAC=/home/jersal/lluis/atac
OUTF=$ATAC/data/20.09.21

cd $ATAC/data
mkdir -p $OUTF

cd $OUTF


##################################################################
## Structure of this script
##################################################################

## Genome Fasta and Annotation files
## Blacklisted regions
## OrgDb
## Chromatin states
## PWMS
## Homer files
## CHIP-Seq files
## Trimmed down test datasets


##################################################################
## Genome Fasta and Annotation files
##################################################################

cd $OUTF ; mkdir genomes
GENF=$OUTF/genomes # GENome Folder
cd $GENF

function get_fasta_and_gff {
  
  release=$1
  dna_type=$2
  specie=$3
  genome=$4
  assembly=$5
  
  mkdir -p $genome/annotation
  mkdir -p $genome/sequence
  
  URL=ftp://ftp.ensembl.org/pub/release-$release/
  
  gff3_file=${specie^}.$genome.$release.gff3.gz
  fasta_file=${specie^}.$genome.$dna_type.$assembly.fa.gz
  
  wget -O - $URL/gff3/$specie/$gff3_file | gunzip -c > $genome/annotation/anno.gff3
  wget -O - $URL/fasta/$specie/dna/$fasta_file | gunzip -c > $genome/sequence/genome.fa
  echo "gff3 file : $gff3_file" > $genome/README
  echo "fasta file: $fasta_file" >> $genome/README
  cd $GENF
  
}

get_fasta_and_gff 100 dna_sm caenorhabditis_elegans WBcel235 toplevel
get_fasta_and_gff 100 dna_sm drosophila_melanogaster BDGP6.28 toplevel
get_fasta_and_gff 100 dna_sm mus_musculus GRCm38 primary_assembly
get_fasta_and_gff 100 dna_sm homo_sapiens GRCh38 primary_assembly


## in R
cd ~/workspace/cactus/data
R
cur_seq_info = rtracklayer::SeqinfoForUCSCGenome('ce11')
cur_seq_info@seqnames %<>% gsub('chr', '', .)
saveRDS(cur_seq_info, '~/workspace/cactus/data/worm/genome/sequence/cur_seqinfo.rds')


##################################################################
## Blacklisted regions
##################################################################

## list of blacklisted regions for the Encode species
# https://github.com/Boyle-Lab/Blacklist/tree/master/lists

## a tool to convert names between difference sources of the same genome
# https://github.com/dpryan79/ChromosomeMappings
# conda install -c bioconda cvbio

# cvbio UpdateContigNames -i ce11_blacklist_v2.bed -o ce11_blacklist_v2__Ensembl_named.bed -m WBcel235_UCSC2ensembl.txt --comment-chars '#' --columns 0 --skip-missing true
# # note on this tool:
# # --comment-chars '#': means line starting with a # will be kept as they are
# # --columns 0: means the first column contains the chromosome names to modify
# # --skip-missing true: means lines without corresponding chr names will be excluded from the new file (it may be worth to use wc -l to check if not too many lines were silently skipped when setting this argument to true)

### lastest ensembl version of each species
# https://www.ensembl.org/info/about/species.html

## correspondances between genome versions from Ensembl and UCSC
# https://hgdownload.cse.ucsc.edu/downloads.html


ATAC="/home/jersal/lluis/atac"

cd $ATAC/tools/scripts ; mkdir chromosome_mappings
CMP=$ATAC/tools/scripts/chromosome_mappings

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
# 97 ce11_blacklist_v2.bed
# 97 ce11_blacklist_v2__Ensembl_named.bed
# 182 dm6_blacklist_v2.bed
# 182 dm6_blacklist_v2__Ensembl_named.bed
# 636 hg38_blacklist_v2.bed
# 636 hg38_blacklist_v2__Ensembl_named.bed
# 3435 mm10_blacklist_v2.bed
# 3435 mm10_blacklist_v2__Ensembl_named.bed
# 8700 total

## a note about drosophila genome assemblies
# https://www.ncbi.nlm.nih.gov/assembly/organism/7227/latest/
# dm6 is probably not what we want... it is a release from 204 
# Release 6 plus ISO1 MT



## correspondances between genome versions from Ensembl and UCSC
# https://hgdownload.cse.ucsc.edu/downloads.html
# 
# Human genome
#  Dec. 2013 (GRCh38/hg38) 
#  Feb. 2009 (GRCh37/hg19)
#  Mar. 2006 (NCBI36/hg18)
# 
# Mouse genome
#  Dec. 2011 (GRCm38/mm10)
#  Jul. 2007 (NCBI37/mm9)
# 
# D. melanogaster genome
#  Aug. 2014 (BDGP Release 6 + ISO1 MT/dm6)
#  Apr. 2006 (BDGP R5/dm3)
# 
# C. elegans genome
#  Feb. 2013 (WBcel235/ce11) 
#  Oct. 2010 (WS220/ce10) 
#  May 2008  (WS190/ce6)



# NCBI=(   ce11    dm6   mm10   hg38  )




##################################################################
## OrgDb
##################################################################

# Note: I do it manually (and not in a loop or in the nextflow script), since Bioconductor ask me to update hundred of packages, and I don't know how to load a library with a variable

# orgDb packages are not really related to a specific genome version
# https://support.bioconductor.org/p/84593/#84594
# The orgDb packages don't really contain any positional annotation. They used to, but these days you will be directed to a TxDb package if you try to get positional info. And the TxDb have the build in the package name. The orgDb packages mostly contain mappings between various databases and some functional annotation, none of which is based on any build. In fact, most of that stuff is updated weekly or monthly, so the orgDb packages get outdated to a certain extent rather quickly.
# The only data in the orgDb packages that comes from hg19 is the outdated CHRLOC table, and as I already mentioned, any query to that table will result in a message saying you should use a TxDb package for that info. The remaining data are not based on any build, as things like an Entrez Gene ID are build agnostic. So to say it is based on any build is misleading.

ATAC="/home/jersal/lluis/atac"

cd $OUTF ; mkdir org_db
cd $OUTF/org_db

R

# BiocManager::install("org.Ce.eg.db")
library(org.Ce.eg.db)
# sink('org_db_ce_package_info.txt')
# 	org.Ce.eg.db
# sink()
AnnotationDbi::saveDb(org.Ce.eg.db, file = 'org_db_ce.sqlite')

# BiocManager::install("org.Dm.eg.db")
library(org.Dm.eg.db)
# sink('org_db_dm_package_info.txt')
# 	org.Dm.eg.db
# sink()
AnnotationDbi::saveDb(org.Dm.eg.db, file = 'org_db_dm.sqlite')

# BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
# sink('org_db_mm_package_info.txt')
# 	org.Mm.eg.db
# sink()
AnnotationDbi::saveDb(org.Mm.eg.db, file = 'org_db_mm.sqlite')

# BiocManager::install("org.hs.eg.db")
library(org.Hs.eg.db)
# sink('org_db_hs_package_info.txt')
# 	org.Hs.eg.db
# sink()
AnnotationDbi::saveDb(org.Hs.eg.db, file = 'org_db_hs.sqlite')



##################################################################
## Chromatin states
##################################################################

## Some are available here but not done with ChromHMM and with older genome versions, and mice is missing. But at least it was done systematically.
# https://www.encodeproject.org/comparative/chromatin/#hiHMM

# NCBI=(   ce11    dm6   mm10   hg38  )

# otherwise I could potentially create them myself, using: https://rdrr.io/bioc/coMET/man/chromatinHMMAll_UCSC.html

# or download each of them one by one on http://compbio.mit.edu/ChromHMM/

# the problem is that there are different tissues in mice and humans so different chromatin states for each tissue
# i.e. brain for mouse here: https://www.encodeproject.org/annotations/ENCSR347PMS/
# https://main.genome-browser.bx.psu.edu/cgi-bin/hgTrackUi?g=meryChromHmm7s&db=mm9

# => for humans and mice: the user should include the file for his cell line.
# => just include a file by default for fly and worm

# => in fact use this list 
https://www.nature.com/articles/nature13415
https://sites.google.com/site/kasohn/software
https://www.encodeproject.org/comparative/chromatin/#hiHMM

# there are just 6 entries but at least they are constitent and there is a nice publication to explain how and why they were created

# # here are the species available:
    human (hg19) - H1-hESC
    human (hg19) - GM12878
    fly (dm3) - LE
    fly (dm3) - L3
    worm (ce10) - EE
    worm (ce10) - L3
		
# => now to do: downloading them, and changing the plotting script accordingly!
# I wanted to add a pvalue anyway

# for mouse, well too bad, there is no file available for now. One can however genererate it if needed.

# make the liftover
# => download the homer files on the server as it crashes all the time in the containers. 

# of note: there are actually many other options of files to use, such as Enhancers,  Chromatin-based inferred topological domains and their boundaries,  Hi-C defined topological domains, Lamina Associated Domains (LADs), Heterochromatin domains... But well the stages are weird and not consistent. 
https://www.encodeproject.org/comparative/chromatin/


LOP=$ATAC/tools/scripts/liftover_chain_files

# ## downloading the liftover files ==>> Uncomment all this paragraph if needed!!
# cd $ATAC/tools/scripts
# wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
# chmod +700 liftOver
# 
# # LOP stands for Lift Over chain files Path
# 
# cd $LOP 
# wget https://hgdownload.cse.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz
# wget https://hgdownload.cse.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDm6.over.chain.gz
# wget https://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz
# wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz


## initialization
cd $OUTF ; mkdir chromatin_states
cd $OUTF/chromatin_states ; mkdir raw modif lifted split


## downloading raw chromatin state files 
cd $OUTF/chromatin_states/raw

for BED in iHMM.M1K16.fly_EL.bed iHMM.M1K16.fly_L3.bed iHMM.M1K16.human_GM.bed iHMM.M1K16.human_H1.bed iHMM.M1K16.worm_EE.bed iHMM.M1K16.worm_L3.bed
 do wget http://compbio.med.harvard.edu/modencode/webpage/hihmm/$BED
done

## modifying names and custom files
cp iHMM.M1K16.worm_EE.bed ../modif/worm_EE_ce10.bed 
cp iHMM.M1K16.worm_L3.bed ../modif/worm_L3_ce10.bed 
cp iHMM.M1K16.fly_EL.bed ../modif/fly_EL_dm3.bed 
cp iHMM.M1K16.fly_L3.bed ../modif/fly_L3_dm3.bed 
cp iHMM.M1K16.human_GM.bed ../modif/human_GM12878_hg19.bed 
cp iHMM.M1K16.human_H1.bed ../modif/human_H1_hESC_hg19.bed 

gawk -i inplace '{print "chr"$0}' ../modif/worm_EE_ce10.bed 
gawk -i inplace '{print "chr"$0}' ../modif/worm_L3_ce10.bed 


## remaping genomic coordinates to a newer genome version
cd $OUTF/chromatin_states/lifted
mkdir unmapped

$LOP/../liftOver ../modif/worm_EE_ce10.bed  $LOP/ce10ToCe11.over.chain.gz worm_EE_ce11.bed unmapped/worm_EE_ce11_unmapped.bed
$LOP/../liftOver ../modif/worm_L3_ce10.bed  $LOP/ce10ToCe11.over.chain.gz worm_L3_ce11.bed unmapped/worm_L3_ce11_unmapped.bed
 
$LOP/../liftOver ../modif/fly_EL_dm3.bed  $LOP/dm3ToDm6.over.chain.gz fly_EL_dm6.bed unmapped/fly_EL_dm6_unmapped.bed 
$LOP/../liftOver ../modif/fly_L3_dm3.bed  $LOP/dm3ToDm6.over.chain.gz fly_L3_dm6.bed unmapped/fly_L3_dm6_unmapped.bed 

$LOP/../liftOver ../modif/human_GM12878_hg19.bed  $LOP/hg19ToHg38.over.chain.gz human_GM12878_hg38.bed unmapped/human_GM12878_hg38_unmapped.bed 
$LOP/../liftOver ../modif/human_H1_hESC_hg19.bed  $LOP/hg19ToHg38.over.chain.gz human_H1_hESC_hg38.bed unmapped/human_H1_hESC_hg19_unmapped.bed 


## splitting bed files with one file for each chromatin state
cd $OUTF/chromatin_states/split

cp ../lifted/*.bed .
awk ' { FOLDERN=FILENAME; gsub(".bed", "", FOLDERN); "mkdir -p " FOLDERN | getline ; print >  FOLDERN"/"$4".bed" }' *.bed 
rm *.bed

# removing the unmapped files
find . -name "17_Unmap.bed" -delete

## checking final files
# number of files in each subdirectory
du -a | cut -d/ -f2 | sort | uniq -c | sort -nr
17 worm_L3_ce11
17 worm_EE_ce11
17 human_H1_hESC_hg38
17 human_GM12878_hg38
17 fly_L3_dm6
17 fly_EL_dm6
 1 352744  .
 
# number of lines in each subdirectory
find . -type f -print0 | wc -l --files0-from=- > nb_of_lines_in_each_bed_file.plain




##################################################################
## PWMS
##################################################################

# http://cisbp.ccbr.utoronto.ca/bulk.php


# description of homer .motif files
# http://homer.ucsd.edu/homer/motif/
# http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html

# I probably also need to convert motifs from PFMs to PWMs format 
# https://github.com/GreenleafLab/chromVARmotifs
# but no, it is pwm: http://cisbp.ccbr.utoronto.ca/faq.html

# this tool was updated recently (2019) with a Nature Genetics publication describing the updates
# http://cisbp.ccbr.utoronto.ca/update.html

## getting the scores for homer custom motifs files (email correspondance with team making CisBP)
# I would like to use cisbp motifs with Homer (for worm, fly, mouse and humans). However, for homer one needs to set the Log odds detection threshold, which specifies whether a given sequence is enough of a "match" to be considered recognized by the motif.  (see http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html)
# We usually set the cutoff to 70% of the maximum possible log odds score of each motif matrix. One thing you need to be careful is that Homer uses natural log instead of log2 in its calculation. Let me know if you have any other questions.


cd $OUTF ; mkdir motifs_PWMS
cd $OUTF/motifs_PWMS ; mkdir full_database
cd $OUTF/motifs_PWMS/full_database

## when downloading the whole data base the plus motifs are missing so we skip that for now

# wget http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/Bulk_downloads/EntireDataset/PWMs.zip
# wget http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/Bulk_downloads/EntireDataset/TF_Information_all_motifs.txt.zip
# unzip PWMs.zip
# unzip TF_Information_all_motifs.txt.zip
# 
# # removing dollars to preventing the crash of R scripts
# sed -i 's/#/|/g' TF_Information_all_motifs.txt
# 
# wc -l TF_Information_all_motifs.txt # 10,879,323 entries
# 
# head -1 TF_Information_all_motifs.txt > TF_Information_all_motifs_encode_species.txt
# grep -E "(Caenorhabditis_elegans|Drosophila_melanogaster|Mus_musculus|Homo_sapiens)" TF_Information_all_motifs.txt >> TF_Information_all_motifs_encode_species.txt
# 
# wc -l TF_Info*.txt
# # 10,879,323 TF_Information_all_motifs.txt
# #     48,152 TF_Information_all_motifs_encode_species.txt
# 
# grep Homo_sapiens TF_Information_all_motifs.txt | wc -l 
# cut -f8 TF_Information_all_motifs_encode_species.txt | sort | uniq -c
# 
# 
# cut -f8 TF_Information_all_motifs_encode_species.txt | sort | uniq -c
# # 20591 Caenorhabditis_elegans
# #  2475 Drosophila_melanogaster
# # 12361 Homo_sapiens
# # 12724 Mus_musculus
# 
# cd $OUTF/motifs_PWMS/full_database

## comparing with the old version:
cut -f9 ../../../motifs_PWMS/fly/TF_Information_all_motifs.txt | sort | uniq -c 
# 1273 D
#  908 I
#  294 N

cut -f9 ../../../motifs_PWMS/fly/TF_Information_all_motifs_plus.txt | sort | uniq -c 
#  1273 D
# 44256 I
#   294 N

## => the whole database doesn't contain all motifs! So we copy the other folders


cd $OUTF/motifs_PWMS
cp -r ../../motifs_PWMS/worm . 
cp -r ../../motifs_PWMS/fly . 
cp -r ../../motifs_PWMS/mice . 
cp -r ../../motifs_PWMS/human . 

## the code below is the old code how I got these files

# wget http://cisbp.ccbr.utoronto.ca/tmp/Caenorhabditis_elegans_2020_08_05_9:05_am.zip
# wget http://cisbp.ccbr.utoronto.ca/tmp/Drosophila_melanogaster_2020_08_05_9:07_am.zip
# wget http://cisbp.ccbr.utoronto.ca/tmp/Mus_musculus_2020_08_05_9:07_am.zip
# wget http://cisbp.ccbr.utoronto.ca/tmp/Homo_sapiens_2020_08_05_9:09_am.zip
# 
# unzip Caenorhabditis_elegans_2020_08_05_9:05_am.zip -d worm
# unzip Drosophila_melanogaster_2020_08_05_9:07_am.zip -d fly
# unzip Mus_musculus_2020_08_05_9:07_am.zip -d mouse
# unzip Homo_sapiens_2020_08_05_9:09_am.zip -d human
# 
# # need to modify these files as the dollar in the ID make the read.table function fail
# sed -i 's/#/|/g' human/TF_Information.txt
# sed -i 's/#/|/g' human/TF_Information_all_motifs.txt
# 
# R
# 
# library(magrittr)
# library(dplyr)
# 
# 
# get_rank_per_unique_character <- function(v_char){
#   for(char in unique(v_char)){
#     sel = which(v_char == char)
#     v_char[sel] = seq_len(length(sel))
#   }
#   v_char
# }
# 
# get_new_name_by_unique_character_rank <- function(v_char){
# 	u_char = get_rank_per_unique_character(v_char)
# 	v_char = paste(v_char, u_char, sep = '_')
#   v_char = gsub('_1$', '', v_char)
# 	v_char
# }
# 
# ## 3 files are available
# # TF_Information
# # TF_Information_all_motifs
# # TF_Information_all_motifs_plus
# 
# for (specie in c('worm', 'fly', 'mouse', 'human')){
# 
#   filen = paste0(specie, '/TF_Information.txt')
#   meta = read.table(filen, sep = '\t', header = T, stringsAsFactors = F)
#   meta1 = meta %>% filter(Motif_ID != '.')
#   meta1 %<>% filter(!duplicated(Motif_ID))
#   meta1$TF_Name_unique = get_new_name_by_unique_character_rank(meta1$TF_Name)
# 
#   saveRDS(meta1, paste0(specie, '/metadata_filtered.rds'))
# 
#   homer_motifs_file = paste0(specie, '/homer_motifs.txt')
# 
#   sink(homer_motifs_file)
# 
#   for(c1 in 1:nrow(meta1)){
#     TF_Name_unique = meta1$TF_Name_unique[c1]
#     Motif_ID = meta1$Motif_ID[c1]
#     motif = read.table(paste0(specie, '/pwms_all_motifs/', Motif_ID,'.txt'), sep = '\t', header = T, stringsAsFactors = F)
#     motif1 = motif[, -1]
#     log2_odd_score = sum(apply(motif[,-1], 1, function(x) log2(max(x)/0.25)))
#     if(log2_odd_score == 0) next() # some motifs are empty
#     homer_threshold = 0.7 * log2_odd_score
#     consensus_sequence = paste0(colnames(motif1)[apply(motif[,-1], 1, function(x) which.max(x))], collapse = '')
# 
#     cat(paste0('>', consensus_sequence, '\t', TF_Name_unique, '\t', homer_threshold, '\n'))
#     cat(paste0(apply(motif1, 1, paste0, collapse = '\t'), '\n'))
# 
# 
#   }
#   sink()
# 
# }
# 
# # issue: 
# # using TF_Information_all_motifs_plus: there are almost 400 motifs for zfh-2 in worms. 
# 
# # using TF_Information_all_motifs_plus:
# grep zfh-2 worm/homer_motifs.txt | wc -l # 393
# wc -l  worm/homer_motifs.txt # 41397 
# 
# # using TF_Information_all_motifs:
# grep zfh-2 worm/homer_motifs.txt | wc -l # 1
# wc -l  worm/homer_motifs.txt # 32920 
# 
# # still too many motifs for certain TF. i.e. crh-1
# grep crh-1 worm/homer_motifs.txt | wc -l # 34
# 
# # using TF_Information:
# grep zfh-2 worm/homer_motifs.txt | wc -l # 1
# wc -l  worm/homer_motifs.txt # 3444 
# 
# # => let's use this one for now!
# 
# # meta1[1:5, c('TF_Name_unique', 'TF_Name', 'DBID', 'Motif_ID', 'DBID.1', 'TF_Status', 'Motif_Type')]
# # http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html
# 
# 
# awk 'BEGIN { OFS = "\t" } ; { if (substr($1,1,1) == ">") {odd_score = $3 * 100 / 70; $3 = odd_score * 50 / 100} ; print $0 } '  homer_motifs.txt > homer_motifs_threshold_50_percent.txt


##################################################################
## Homer files
##################################################################

# => this should be done with the custom container

## Here one can download organisms, promoters and genomes
# http://homer.ucsd.edu/homer/download.html
# http://homer.ucsd.edu/homer/data/promoters/
# http://homer.ucsd.edu/homer/data/genomes/
# http://homer.ucsd.edu/homer/data/organisms/

# maybe I can even use the galaxy container...

ATAC="/home/jersal/lluis/atac"

URL="http://homer.ucsd.edu/homer/data"

cd $OUTF ; mkdir -p homer_files

cd $OUTF/homer_files

ls ..

cp -r ../../homer_files/* .

## original code

# SPECIES=(worm fly mouse human)
# GENOMES=(ce11 dm6 mm10 hg38)
# 
# for SPECIE in ${SPECIES[@]} ; do wget -nc $URL/promoters/$SPECIE.v5.5.zip ; done
# for SPECIE in ${SPECIES[@]} ; do wget -nc $URL/organisms/$SPECIE.v6.3.zip ; done
# for GENOME in ${GENOMES[@]} ; do wget -nc $URL/genomes/$GENOME.v6.4.zip ; done
# 
# unzip \*.zip
# 
# 
# 
# ATAC=/home/jersal/lluis/atac
# cd $OUTF/homer_files
# wget http://homer.ucsd.edu/homer/custom.motifs
# 
# # # we rename TFs this way:
# # AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer => AP-1_1
# # ZNF652/HepG2-ZNF652.Flag-ChIP-Seq_5626 => ZNF652_2
# 
# awk  'BEGIN { FS = "\t" ; OFS = "\t" ; sum = 0 } { if(substr($1, 1, 1) == ">") { sum += 1 ; gsub("\\(.*", "", $2) ; gsub("\\/.*", "", $2) ; $2 = $2"_"sum } ; print $0 }' custom.motifs > custom_modified.motifs.motifs
# 
# # pb: many TFs are missing (only 416 of them)
# c_homer
# cp /usr/local/share/homer-4.11-2/data/knownTFs/known.motifs .
# # => the right file is there: /usr/local/share/homer-4.11-2/data/knownTFs/known.motifs
# # wc -l /usr/local/share/homer-4.11-2/data/knownTFs/known.motifs # => 12928
# # wc -l custom.motifs # => 5743
# 
# awk  'BEGIN { FS = "\t" ; OFS = "\t" ; sum = 0 } { if(substr($1, 1, 1) == ">") { sum += 1 ; gsub("\\(.*", "", $2) ; gsub("\\/.*", "", $2) ; $2 = $2"_"sum } ; print $0 }' known.motifs > known_modified.motifs





##################################################################
## getting homer_TF_ids_tab_rds (or not)
##################################################################

# This file is used for selecting candidates. That is genes that overlap ATAC and mRNA and that also are Transcription Factors (present in the list used in homer analysis). However, this probably can be removed... Let's see at the end.





##################################################################
## CHIP-Seq files
##################################################################

# # on modEncode there are CHIP-Seq data for only two species: c elegans (209 chip) and d megalonaster (133 CHIP)
# http://data.modencode.org/?Organism=D.%20melanogaster
# 
# https://github.com/kundajelab/ENCODE_downloader
# 
# # filtering chip seq datasets of interest
# https://www.encodeproject.org/search/?type=File&output_type=peaks&assay_term_name=ChIP-seq&file_type=bed+narrowPeak&status!=revoked&target.label!=H3K4me3&target.label!=H3K4me1&target.label!=H3K27me3&target.label!=H3K36me3&target.label!=H3K27ac&target.label!=H3K9me3&target.label!=H3K9ac&target.label!=H3K4me2
# 
# # => doing that it appears that there are mostly Histone CHIP seq in humans.
# 
# # need to find another way to download systematically transcription factors CHIP Seq data
# https://www.factorbook.org/
# 
# # in fact, here we can filter for TF_CHIP_Seq data
# https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&files.file_type=bigBed+narrowPeak
# 
# https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&files.file_type=bigBed+narrowPeak&limit=all

# # The python downloader tool can work with an encode url such as this one:
# $ python encode_downloader.py "https://www.encodeproject.org/search/?type=Experiment&assay_title=ATAC-seq&replicates.library.biosample.life_stage=postnatal&status=released" --ignore-unpublished --file-types fastq "bam:unfiltered alignments"
# 
# # => TO DO: make the query to download all CHIP-Seq files for each species
# 
# ATAC="/home/jersal/lluis/atac"
# 
# cd $OUTF ; mkdir encode_CHIP
# cd $OUTF/encode_CHIP ; mkdir worm fly mouse human
# 
# pip install requests
# cd $ATAC/tools/scripts
# wget https://raw.githubusercontent.com/kundajelab/ENCODE_downloader/master/encode_downloader.py
# chmod +700 encode_downloader.py
# 
# SP=$ATAC/tools/scripts
# 
# cd $OUTF/encode_CHIP/worm
# 
# python $SP/encode_downloader.py "https://www.encodeproject.org/search/?type=File&output_type=peaks&assay_term_name=ChIP-seq&file_type=bed+narrowPeak&status!=revoked&target.label!=H3K4me3&target.label!=H3K4me1&target.label!=H3K27me3&target.label!=H3K36me3&target.label!=H3K27ac&target.label!=H3K9me3&target.label!=H3K9ac&target.label!=H3K4me2&assembly=ce10
# " --ignore-unpublished 
# 
# --file-types bed "bam:unfiltered alignments"
# 
# cd $OUTF/encode_CHIP/mouse
# python $SP/encode_downloader.py  "https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&files.file_type=bigBed+narrowPeak&limit=all&assembly=mm10" --ignore-unpublished  --file-types "bed:narrowPeak"


# new approach!! 
# using: ENCODExplorer
# https://www.bioconductor.org/packages/release/bioc/vignettes/ENCODExplorer/inst/doc/ENCODExplorer.html#obtaining-consensus-peaks-from-chip-seq

cd $OUTF ; mkdir encode_CHIP
cd $OUTF/encode_CHIP #; mkdir worm fly mouse human

cp -r ../../encode_CHIP/* .


R

# BiocManager::install("ENCODExplorer")

library(magrittr)
devtools::load_all('/home/jersal/tools/homemade_R_packages/jeR')


library(ENCODExplorer)
encode_df <- get_encode_df()
# /home/jersal/.cache/AnnotationHub
#   does not exist, create directory? (yes/no): y
# snapshotDate(): 2019-10-29
# downloading 1 resources
# retrieving 1 resource

dim(encode_df) # [1] 454413     63

df_chip = encode_df %>% .[.$assay == 'TF ChIP-seq', ] %T>% pnrow # 141085

df_chip1 = df_chip %>% .[.$file_type == 'bed narrowPeak', ] %T>% pnrow # 30776

df_chip2 = df_chip1 %>% .[.$output_type == 'optimal IDR thresholded peaks', ] %T>% pnrow # 7410

table(df_chip2$assembly)
 ce10   ce11    dm3    dm6 GRCh38   hg19   mm10
  721    801    461    788   1695   2773    171

# paste(grep('chip', tolower(unique(encode_df$assay)), value = T), collapse = '   ')
# # [1] "control chip-seq   histone chip-seq   tf chip-seq   repli-chip   rip-chip"


######################################

df = df_chip2[, c('file_accession', 'assembly', 'target', 'biosample_name', 'biosample_type', 'dataset_biosample_summary', 'dataset_description', 'organism')]


######################################

## creating my own name for the CHIP files
targets = gsub('(eGFP-|3xFLAG-)', '', df$target)

samples = gsub(' ', '_', df$biosample_name)


# # samples %<>% gsub('whole_organism', 'WOGN', .)
# quant10(sapply(targets, nchar))
# # 0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100%
# #  1    3    4    4    5    5    5    6    6    7   15
# quant10(sapply(samples, nchar))
# # 0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100%
# #  2    4    4    4    4    4    4    5    6    7   36
# # => most targets and samples have less than 8 characters so we can cut at 7 charaters

targets %<>% substr(., 1, 7)

# rev(sort(table(targets)))[1:20]
# # CTCF  POLR2A POLR2Ap    REST    EZH2     MYC     YY1   RAD21   elt-3   EP300
# #  410     120      66      55      49      47      45      45      44      43
# # snpc-4    JUND   GABPA  unc-62    TAF1   SIN3A   CEBPB     JUN   FOXA1    NRF1
# #   39      39      37      36      35      35      31      30      30      29


samples %<>% gsub('.*', 'TO_FILL', .) 

# editing the worm samples
is_ce = df$organism == 'Caenorhabditis elegans'
bs = df$dataset_biosample_summary[is_ce]
bs1 = bs %>% gsub('.*whole organism ', '', .) %>% gsub('.*hermaphrodite ', '', .) %>% gsub('(L[1-4]).*', '\\1', .) %>% gsub('young adult.*', 'YA', .) %>% gsub('late embryo.*', 'LE', .) %>% gsub('early embryo.*', 'EE', .) %>% gsub('mixed stage \\(embryo\\).*', 'MSE', .) %>% gsub('dauer.*', 'DA', .) %>% gsub('midembryo.*', 'ME', .)
# sort(table(bs1))
# # DA  ME  EE MSE  L2  L3  LE  YA  L4  L1
# #  6   8  34  98 168 175 217 249 254 313
samples[is_ce] = bs1


# editing the fly samples
is_dm = df$organism == 'Drosophila melanogaster'
bs = df$dataset_biosample_summary[is_dm]
bs1 = bs %>% gsub('.*whole organism male adult.*', 'MA', .) %>% gsub('.*whole organism female adult.*', 'FA', .) %>% gsub('.*whole organism adult.*', 'A', .) %>% gsub('.*whole organism embryo.*', 'E', .) %>% gsub('.*whole organism prepupa.*', 'PP', .) %>% gsub('.*whole organism wandering third.*', 'WT', .) %>% gsub('.*whole organism pupa.*', 'P', .) %>% gsub('.*white prepupa stage.*', 'WP', .) %>% gsub('.*strain Oregon-R S2.*', 'SO', .) %>% gsub('.*ovary female.*', 'OF', .) %>% gsub('.*whole organism larva.*', 'L', .)   
# sort(table(bs1))
# # L    OF    SO     A Kc167    WP    MA    FA     P    WT    PP     E
# # 1     1     6     8    14    18    22    23    23   113   167   851
# # [1] "E"     "PP"    "WT"    "P"     "FA"    "A"     "WP"    "Kc167" "MA"
# # [10] "OF"    "SO"    NA      "L"
samples[is_dm] = bs1


# editing the human samples
is_hs = df$organism == 'Homo sapiens'
bs = df$biosample_name[is_hs]
# rev(sort(table(bs)))
bs1 = sapply(bs, function(x) strsplit(x, ' ')[[1]][1])
# rev(sort(table(bs1)))[1:10]
# # K562   HepG2  HEK293 GM12878   MCF-7      H1    A549 HeLa-S3   liver SK-N-SH
# # 1270     695     462     462     293     185     170     115     104      52
samples[is_hs] = bs1


# editing the mouse samples
is_mm = df$organism == 'Mus musculus'
bs = df$biosample_name[is_mm]
bs1 = sapply(bs, function(x) strsplit(x, ' ')[[1]][1])
# rev(sort(table(bs1)))[1:10]
# # MEL CH12.LX     G1E G1E-ER4   liver    lung   heart  ES-E14    bone  spleen
# #  49      41      10       9       7       5       5       5       4       3
samples[is_mm] = bs1

grep('TO_FILL', samples) # integer(0)

samples %<>% substr(., 1, 7)

targets_samples = paste(targets, samples, sep = '_')

df[, name := targets_samples]
data.table::setcolorder(df, c(ncol(df), seq_len(ncol(df) - 1)))

grep('_', targets) # integer(0)
grep('_', samples) # integer(0)
# => so the sign "_" can be safely used as a separator
grep('^[0-9]+$', samples, value = T) # no samples is composed only of digits

nrow(df) # 7410

saveRDS(df, 'df_metadata_encode_chip.rds')


# here are the most tested conditions
rev(sort(table(df$name)))[1:12]
   # elt-3_L1  POLR2A_K562 CTCF_fibrobl    snpc-4_YA    IRF1_K562   CTCF_MCF-7
   #       28           25           23           20           19           18
   # MYC_K562    Su(var)_E   NR3C1_A549     JUN_K562    E(spl)m_E    ceh-82_L1
   #       17           16           16           15           14           14


R

library(data.table)
library(magrittr)
devtools::load_all('/home/jersal/tools/homemade_R_packages/jeR')

df = readRDS('df_metadata_encode_chip.rds') %T>% pnrow # 7410

# one biosample_name has a / in it (NT2/D1). We remove it has it crashes the file names paths
df[file_accession %in% c('ENCFF597KMH', 'ENCFF226OCL', 'ENCFF577FWL'), 1:5]
df$name %<>% gsub('\\/', '-', .)
df$biosample_name %<>% gsub('\\/', '-', .)
df[file_accession %in% c('ENCFF597KMH', 'ENCFF226OCL', 'ENCFF577FWL'), 1:5]


assemblies = c('ce11', 'dm6', 'mm10', 'GRCh38')

df1 = df[assembly %in% assemblies] %T>% pnrow # 3455
length(which(duplicated(df1$name))) # 1033
# => so one third of the names are duplicated


# here are the most duplicated entries
rev(sort(table(df1$name)))[1:12]
# elt-3_L1    snpc-4_YA    E(spl)m_E    Su(var)_E    IRF1_K562  POLR2A_K562
#       14           10           10            9            9            8
# CTCF_fibrobl    ceh-82_L1   NR3C1_A549    lin-35_YA     efl-1_YA     efl-1_L1
#        7            7            6            6            6            6


new_names = jeR::get_new_name_by_unique_character_rank(df1$name)
new_names %<>% gsub('_1$', '', .)
length(which(duplicated(new_names))) # 0

df1$name = new_names

saveRDS(df1, 'df_metadata_encode_chip_final.rds')
write.csv(df1, 'df_metadata_encode_chip_final.csv')



## downloading the encode bed files

R

library(ENCODExplorer)

df = readRDS('df_metadata_encode_chip_final.rds')[, 1:4]
nr = nrow(df)


assemblies = c('ce11', 'dm6', 'mm10', 'GRCh38')
sapply(assemblies, dir.create, showWarnings = F)

t1 = Sys.time()

for(c1 in 1:nrow(df)){
  dir = df$assembly[c1]
  f_acc = df$file_accession[c1]
  old_bed = paste0(dir, '/', f_acc, '.bed.gz')
  new_bed_unzipped = paste0(dir, '/', df$name[c1], '.bed')
  new_bed = paste0(new_bed_unzipped, '.gz')
  if(!file.exists(new_bed) & !file.exists(new_bed_unzipped)){
    cat(paste0('\n\n ', c1, ' / ', nr, ': ', old_bed, ' -> ', new_bed, '\n\n'))
    downloadEncode(f_acc, dir = dir)
    file.rename(old_bed, new_bed)
  }
}
  
t2 = Sys.time()
t2 - t1

# Time difference of 3.143832 hours

nrow(df) # 3455


## in linux

du -a | cut -d/ -f2 | sort | uniq -c | sort -nr
1696 GRCh38
 802 ce11
 789 dm6
 172 mm10

du -a | grep "bed" | wc -l # 3455 => all good!
 
for assembly in ce11 dm6 mm10 GRCh38; do gunzip $assembly/* ; done

du -h # 2.6G    .



##################################################################
## Trimmed down test datasets 
##################################################################

((mRNA[DataSet Type]) AND ATAC[DataSet Type]) AND C Elegans[Organism] 

((gene expression[DataSet Type]) AND ATAC[DataSet Type]) AND Caenorhabditis elegans[Organism] 
((RNA[Description]) AND ATAC[Description]) AND Caenorhabditis elegans[Organism] 


ATAC-seq 

## worm
# Paper name: Janes et al. 2018, Chromatin accessibility dynamics across C. elegans development and ageing
# Paper link: https://elifesciences.org/articles/37344
# GEO accession number: GSE114494 
# conditions: Emb, L1, l2, L3, L4, D1, D2, D6, D9, D13
# comparisons: L1-Emb, L2-Emb, L3-Emb, L4-Emb, D1-Emb, D2-D1, D6-D1, D9-D1, D13-D1, D6-D13
# comparisons groups: Develop, Adult, All

## fly
# Paper name: Patterns of chromatin accessibility along the anterior-posterior axis in the early Drosophila embryo
# Paper link: https://journals.plos.org/plosgenetics/article?rev=2&id=10.1371/journal.pgen.1007367
# GEO accession number: GSE104957 
# conditions: 
# comparisons: 
# comparisons groups:

## mice
# Paper name: 
# Paper link: 
# conditions: 
# comparisons: 
# comparisons groups:



## human
# Paper name: 
# Paper link: 
# conditions: 
# comparisons: 
# comparisons groups:





