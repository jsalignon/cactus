
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Tutorial](/docs/1_Intro/tutorial.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)

Almost all packages (excepting those with a "/" in their id in the table below) used by Cactus are [Biocontainers](https://biocontainers-edu.readthedocs.io/en/latest/introduction.html). Biocontainer packages have the advantage to be available on conda and Mamba virtual environments, Singularity and Docker containers. As much as possible, single tools Biocontainers were used for easier maintenance. However, in some cases multiple-tools are needed (e.g. analysis in R) in which cases "mulled containers" were created via [BioContainers' multi-package-containers tool](https://github.com/BioContainers/multi-package-containers). The download of all tools happens automatically the first time Cactus is run. Singularity images are hosted on the [Galaxy singularity repository](https://depot.galaxyproject.org/singularity/).   

Below are more details regarding the use of Singularity containers.

Downloading all singularity containers uses 4 Gb of disk space in total. The parameter *singularity_containers_path* set in the *./.cactus.config* file should be set to indicate in which directory the containers should be downloaded.

Here is the detail of the Singularity containers used and their size:


|          name          |  size|id                                                                                            |
|:----------------------:|-----:|:---------------------------------------------------------------------------------------------|
|      pypdf2:2.11.1     |  479M|
|    bowtie2_samtools    |  463M|mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:b6524911af823c7c52518f6c886b86916d062940-0 |
|        figures         |  450M|mulled-v2-e0e17c59e64598cdb16a01c347c673dd021f778a:0dd10d4b12b50eec83ccd1f7b7740a08f2703bdf-0 |
|         homer          |  434M|homer:4.9.1--pl5.22.0_5                                                                       |
| differential_analysis  |  309M|mulled-v2-abcadfb509d8692abe35c0bd02689ab7756d85f8:1b35d287a7c9c53a258be40306bdca167e2e078a-0 |
|      venndiagram       |  306M|pegi3s/r_venn-diagram:1.7.0                                                                   |
|         bbmap          |  293M|bbmap:38.96--h5c4e2a8_0                                                                       |
|        r_basic         |  292M|mulled-v2-b21cd52f0c50bbd777eaed41c0b8228b84cff4bd:b09be1d801d248a5a61257583e629f17052d8181-0 |
|      skewer_pigz       |  291M|mulled-v2-734ede4cc65b3b212388567aac99f6182e023a8f:26fbad413ebdf8aee65d8aa554d52a4f69548508-0 |
|         fastqc         |  261M|fastqc:0.11.7--4                                                                              |
|      bioconductor      |  167M|mulled-v2-0161037c6d8979d1ff5de7e591f5adfb3ffe38b8:2b97ca0a3f4f5409852afe863ef8068a83779815-0 |
|         picard         |  165M|picard:2.26.9--hdfd78af_0                                                                     |
|        diffbind        |  145M|mulled-v2-9ec5efd66a9a09ea4f9ad9bad5485675f031aeb4:cf736786cecad89eca5fea6d119a837e4bad7c08-0 |
|       deeptools        |  107M|deeptools:3.4.3--py_0                                                                         |
| samtools_bedtools_perl |   95M|mulled-v2-95fc59e28f845da0ff950325f8138eff9cedff14:0bc453d1b98bff9aef79c31f643f6b9f93bc7fbd-0 |
|         sleuth         |   55M|r-sleuth:0.30.0--r41hdfd78af_5                                                                |
|         macs2          |   43M|macs2:2.2.7.1--py37hf01694f_1                                                                 |
|        bioperl         |   21M|perl-bioperl-core:1.007002--pl5321hdfd78af_4                                                  |
|        kallisto        |   11M|kallisto:0.46.2--h4f7b962_1                                                                   |
|           -            |     -|-                                                                                             |
|         Total:         | 4081M|-                                                                                             |



Here are more details on the tools in each mulled container:
  - **bowtie2_samtools**: bowtie2=2.4.4, samtools=1.13    
  - **r_basic**: r-base=4.1.3, r-magrittr=2.0.3, r-dplyr=1.0.9, r-purrr=0.3.4, r-ggplot2=3.3.5, r-data.table=1.14.2
  - **samtools_bedtools_perl**: samtools=1.15.1, bedtools=2.30.0, perl=5.32.1
  - **skewer_pigz**: skewer=0.2.2, pigz=2.6
  - **bioconductor**: r-base=4.3, bioconductor-chipseeker=1.36.0, r-magrittr=2.0.3, bioconductor-genomicfeatures=1.52.1, bioconductor-clusterprofiler=4.8.1, bioconductor-annotationdbi=1.62.2, r-purrr=1.0.2, r-ggplot2=3.4.4
  - **diffbind**: bioconductor-diffbind=3.4.0, bioconductor-csaw=1.28.0, bioconductor-edger=3.36.0, r-optparse=1.7.1
  - **figures**: r-base=4.3, r-ggplot2=3.4.4, r-magrittr=2.0.3, r-gridextra=2.3, r-rcolorbrewer=1.1_3, r-data.table=1.14.8
  - **differential_abundance**: r-base=4.1.3, bioconductor-diffbind=3.4.11, r-sleuth=0.30.0, r-ggplot2=3.3.5, r-magrittr=2.0.3, r-openxlsx=4.2.5

<!--

workspace=/home/jersal/workspace
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
singularity_dir=$workspace/singularity_containers/
containers_file=$cactus_dir/conf/containers.config
containers_info_dir=$cactus_dir/docs/util/containers
containers_info=$containers_info_dir/containers_info.tsv
galaxy_info=$containers_info_dir/containers_info_galaxy.tsv
not_galaxy_info=$containers_info_dir/containers_info_not_galaxy.tsv
containers_info_size=$containers_info_dir/containers_info_size.tsv
containers_info_size_1=$containers_info_dir/containers_info_size_1.tsv

grep params $containers_file | grep -v "//" | awk  -v path="$singularity_dir" '{
  name=$1
  gsub("\"", "", $2)
  full_id=id=$2
  gsub("params.", "", name)
  gsub("\\$\\{depot_galaxy\\}\\/", "", id)
  gsub("\\$\\{depot_galaxy\\}\\/", "depot.galaxyproject.org-singularity--", full_id)
  gsub(":", "-", full_id)

  path_full_id= path full_id ".img"
  gsub("params.", "", name)
  print name, id, path_full_id
}' FS=' = ' > $containers_info
grep galaxy $containers_info > $galaxy_info
grep -v galaxy $containers_info > $not_galaxy_info

ls -sh $(grep galaxy $galaxy_info | cut -f3 -d" ") | sed 's/^[[:space:]]*//' | cut -f1 -d" " | paste -d' ' $galaxy_info - | awk '{print $1, $4, $2}' OFS=" " | sort -nrk2 | column -t > $containers_info_size
cat $not_galaxy_info
ls -sh $containers_dir | grep "pegi3s\\|mnuessler" 
echo "venndiagram 306M pegi3s/r_venn-diagram:1.7.0" >> $containers_info_size
echo "pdftk 173M mnuessler/pdftk" >> $containers_info_size
sort -nrk2 $containers_info_size | column -t  > $containers_info_size_1


R
library(data.table)
library(magrittr)
dt = fread('/home/jersal/workspace/cactus/software/docs/util/containers/containers_info_size_1.tsv', sep = ' ', col.names = c('name', 'size', 'id'))
total = gsub('M', '', dt$size) %>% as.integer %>% sum %>% paste0('M')
dt = rbind(dt, data.table(name = '-', size = '-', id = '-'))
dt = rbind(dt, data.table(name = 'Total:', size = total, id = '-'))
knitr::kable(dt, 'pipe', align = c('c', 'r', 'l'))


-->
