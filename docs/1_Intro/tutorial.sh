
# Creating the Cactus configuration file (if needed)
cat >> ${HOME}/.cactus.config << EOL
params.references_dir        = "${HOME}/workspace/cactus/references"
params.singularity_cache_dir = "${HOME}/workspace/singularity_containers"
params.tower_token           = "*"
params.enable_tower          = false
EOL

# Selecting the Nextflow version
export NXF_VER=22.10.8

# Selecting the Cactus version 
cactus_version=0.9.0

# Downloading the worm references and test dataset
nextflow run jsalignon/cactus/scripts/download/download.nf -r $cactus_version \
	--test_datasets --references -profile singularity --species worm

# change folder and check design
cd worm
printf "\n\nATAC-Seq fastq files\n"
cat design/atac_fastq.tsv
printf "\n\nmRNA-Seq fastq files\n"
cat design/mrna_fastq.tsv
printf "\n\nComparisons to be made\n"
cat design/comparisons.tsv
printf "\n\nGroups of comparisons to be plot together in heatmaps\n"
cat design/groups.tsv

# Running Cactus without enrichment analysis
nextflow run jsalignon/cactus -profile singularity -r $cactus_version \
	--res_dir results/tutorial              \
	--species worms                         \
	--chromatin_state iHMM.M1K16.worm_L3    \
	--design__regions_to_remove design/regions_to_remove_empty.tsv \
	--disable_all_enrichments          

# Running Cactus without enrichment analysis and masking the genomic region 
# nearby the RNAi target gene from all ATAC-Seq analyses
nextflow run jsalignon/cactus -profile singularity -r $cactus_version \
	--res_dir results/tutorial              \
	--species worms                         \
	--chromatin_state iHMM.M1K16.worm_L3    \
	--design__regions_to_remove design/regions_to_remove.tsv \
	--disable_all_enrichments      

# Running Cactus with Differential Analysis Subset filters and some enrichment 
# analyses
nextflow run jsalignon/cactus -profile singularity -r $cactus_version \
	--res_dir results/tutorial              \
	--species worms                         \
	--chromatin_state iHMM.M1K16.worm_L3    \
	--design__regions_to_remove design/regions_to_remove.tsv \
	--split__threshold_type     rank        \
	--split__threshold_values [ 200, 1000 ] \
	--do_func_anno_enrichment   false       \
	--do_motif_enrichment       false

# Running the full analysis
nextflow run jsalignon/cactus -profile singularity -r $cactus_version \
	--res_dir results/tutorial              \
	--species worms                         \
	--chromatin_state iHMM.M1K16.worm_L3    \
	--design__regions_to_remove design/regions_to_remove.tsv \
	--split__threshold_type     rank        \
	--split__threshold_values [ 200, 1000 ]

