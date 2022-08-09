

###################################################
### Menus
###################################################


rename -v 's/Intro/1_Intro/' *.md
rename -v 's/Overview/2_Overview/' *.md
rename -v 's/Install/3_Install/' *.md
rename -v 's/ID/4_Run/' *.md
rename -v 's/AP/5_AP/' *.md
rename -v 's/MP/6_MP/' *.md
rename -v 's/DA/7_DA/' *.md
rename -v 's/SP/8_SP/' *.md
rename -v 's/Enrich/9_Enrich/' *.md



* Introduction: [Cactus](/README.md), [Quick Start](/docs/Intro__Quick_start.md), 
* Overview: [Graph](/docs/Overview__Graph.md), [Outputs structure](/docs/Overview__Outputs_structure.md)
* Install: [Dependencies](/docs/Install__Dependencies.md), [Containers](/docs/Install__Containers.md), [Data](/docs/Install__Data.md), [Test_datasets](/docs/Install__Test_datasets.md)
* Run: [Input Data](/docs/ID__Input_data.md), [Input Files](/docs/ID__Input_files.md), [Parameters](/docs/ID__Parameters.md)
* Preprocessing: ATAC: [Method](/docs/AP__Method.md), [Figures](/docs/AP__Figures.md), [MultiQC](/docs/MP__MultiQC.md), mRNA: [Method](/docs/MP__Method.md), [MultiQC](/docs/MP__MultiQC.md)
* Differential Abundance: [ATAC](/docs/DA__DiffBind.md), [mRNA](/docs/DA__Sleuth.md), [Figures](/docs/DA__Figures.md), [Tables](/docs/DA__Tables.md)
* Splitting peak sets: [Split](/docs/SP__Splitting_peak_sets.md), [Venn Diagrams](/docs/SP__Venn_diagrams.md)
* Enrichment: [ChIP](/docs/Enrich__CHIP.md), [Chromatin state](/docs/Enrich__Chromatin_states.md), [Motifs](/docs/Enrich__Motifs.md), [Func. Anno.](/docs/Enrich__Functional_annotations.md), [Figures](/docs/Enrich__Figures.md), [Tables](/docs/Enrich__Tables.md)


ls -1 *.md
# 1_Intro__Quick_start.md
# 2_Overview__Graph.md
# 2_Overview__Outputs_structure.md
# 3_Install__Containers.md
# 3_Install__Data.md
# 3_Install__Dependencies.md
# 3_Install__Test_datasets.md
# 4_Run__Input_data.md
# 4_Run__Input_files.md
# 4_Run__Parameters.md
# 5_AP__Figures.md
# 5_AP__Method.md
# 5_AP__MultiQC.md
# 6_MP__Method.md
# 6_MP__MultiQC.md
# 7_DA__DiffBind.md
# 7_DA__Figures.md
# 7_DA__Sleuth.md
# 7_DA__Tables.md
# 8_SP__Splitting_peak_sets.md
# 8_SP__Venn_diagrams.md
# 9_Enrich__CHIP.md
# 9_Enrich__Chromatin_states.md
# 9_Enrich__Figures.md
# 9_Enrich__Functional_annotations.md
# 9_Enrich__Motifs.md
# 9_Enrich__Tables.md





###################################################
### Get test datasets
###################################################

## could not find a bash way to get markdown tables from tsv. Done in R in the end
# # markdown_table.sh comes from here:
# # https://josh.fail/2022/pure-bash-markdown-table-generator/
# 
# input_file=~/workspace/cactus/test_datasets/preprocessing/fly/samplesheet/samples_info_1.tsv
# cat $input_file | paste -d, - - | pandoc -f csv -t markdown_phpextra
# cat $input_file | pandoc --metadata-file -t markdown_phpextra
# cat $input_file | pandoc  -t markdown_phpextra
# 
# pandoc --filter pandoc-csv2table -o test.html input_file


get_markdown_table <- function(specie){
  tmp = read.table(paste0('~/workspace/cactus/test_datasets/preprocessing/', specie, '/samplesheet/samples_info_1.tsv'), header = T)
  knitr::kable(tmp, 'pipe', align = 'c')
}

get_markdown_table('worm')
get_markdown_table('fly')
get_markdown_table('human')
get_markdown_table('mouse')

setwd('~/workspace/cactus/test_datasets')

# specie = 'worm'
# df = data.frame(
#   raw_files = system(paste0('du -h -d0 ', 'preprocessing/', specie, '/fastq')),
#   sampled_files = system(paste0('du -h -d0 ', specie, '/data')),
#   sampled_atac = system(paste0('du -h -d1 ', specie, '/data/atac')),
#   sampled_mrna = system(paste0('du -h -d1 ', specie, '/data/mrna')),
#   stringsAsFactors = F
# )
# 
#  => issue: cannot capture the printed output from du:
# du -h -d1 **/data/ > preprocessing/report/test_datasets_sizes.txt
# du -h -d1 preprocessing/**/fastq >> preprocessing/report/test_datasets_sizes.txt
#  a <- capture.output(a <- system(paste0('du -h -d0 ', 'preprocessing/', specie, '/fastq')))
# system('powershell -noprofile -command "ls -r|measure -s Length"')
# a = system(paste0('du -h -d0 ', 'preprocessing/', specie, '/fastq'), stdout=TRUE, stderr=TRUE)


dir_size <- function(path, recursive = TRUE) {
  stopifnot(is.character(path))
  files <- list.files(path, full.names = T, recursive = recursive)
  vect_size <- sapply(files, function(x) file.size(x))
  size_files <- sum(vect_size)
  order = ifelse(size_files > 10^9, 'GB', 'MB')
  if(order == 'MB') size_files = paste0(format(size_files/10^6, digits = 2), ' MB')
  if(order == 'GB') size_files = paste0(format(size_files/10^9, digits = 2), ' GB')
  size_files
}

get_test_dataset_size <- function(specie){
  df = data.frame(
    specie = specie,
    raw_files = dir_size(paste0('preprocessing/', specie, '/fastq')),
    sampled_files =dir_size(paste0(specie, '/data')),
    sampled_atac = dir_size(paste0(specie, '/data/atac')),
    sampled_mrna = dir_size(paste0(specie, '/data/mrna')),
    stringsAsFactors = F
  )
  return(df)
}

df = do.call(rbind, lapply(c('fly', 'worm', 'mouse', 'human'), get_test_dataset_size))

knitr::kable(df, 'pipe', align = 'c')





###################################################
### Get data
###################################################

cd /home/jersal/workspace/cactus/data/human

tree -d -L 3 --filelimit 5
.
├── bowtie2_indexes_conta
├── CHIP
├── chromatin_states [835 entries exceeds filelimit, not opening dir]
├── genome
│   ├── annotation
│   │   ├── bed_regions
│   │   ├── filtered
│   │   └── R
│   └── sequence
│       └── bowtie2_indexes
└── homer_data
    ├── accession
    ├── genomes
    │   └── hg38
    ├── GO
    └── promoters











