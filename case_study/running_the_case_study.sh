

#################
## Initialization

## setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus # to change if needed
test_dir=$cactus_dir/test_datasets   # to change if needed
case_study_dir=$test_ds_dir/case_study/data
case_study_dir=$test_ds_dir/application_note

## setting up variables
cpu_nb=30
memory_size='300G'
cactus_version="v0.8.5" ; latest_flag=''
cactus_version=updating_the_bioconductor_container ; latest_flag=''

# nextflow drop jsalignon/cactus # if needed
figshare_version=v4
references_dir=$cactus_dir/references/$figshare_version
padj_breaks="c( 0.2, 1e-5, 1e-50, 1e-100, 1e-200 )"
heatmaps_filters_chip="c( 30, 20, 10, T, 2, 'ward.D', T )"
heatmaps_filters_motifs="c( 30, 20, 10, T, 2, 'ward.D', F )"
heatmaps_params__peaks_self="c( 0.05, T, 'none', T, 50, 'UUDD', 2.5 )"
heatmaps_ggplot__chrom_states="c( 8, 10, 7 )"
heatmaps_ggplot="c( 12, 12, 12 )" 
res_prefix="results/2023_11_26__"

## Initialization
#################


#########################
## Running the case study

# note that only the last 2 argument change between the 3 run commands

# worm run 1: full analysis
species=worm ; run='1' ; cd $case_study_dir/$species
nextflow run jsalignon/cactus \
	-r $cactus_version $latest_flag \
	-profile singularity \
	--references_dir $references_dir \
	--res_dir "${res_prefix}_${species}_run_${run}" \
	--executor_local_cpus $cpu_nb \
	--executor_local_memory $memory_size \
	--common_padj_breaks "$padj_breaks" \
	--heatmaps_filter__CHIP "$heatmaps_filters_chip" \
	--heatmaps_filter__motifs "$heatmaps_filters_motifs" \
	--heatmaps_params__peaks_self  "$heatmaps_params__peaks_self" \
	--common__heatmaps_ggplot "$heatmaps_ggplot" \
	--heatmaps_ggplot__chrom_states "$heatmaps_ggplot__chrom_states" \
	--species $species \
	--chromatin_state 'iHMM.M1K16.worm_L3' \
	--chip_ontology 'all'

# human run 1: full analysis
species=human ; run='1' ; cd $case_study_dir/$species
nextflow run jsalignon/cactus \
	-r $cactus_version $latest_flag \
	-profile singularity \
	--references_dir $references_dir \
	--res_dir "${res_prefix}_${species}_run_${run}" \
	--executor_local_cpus $cpu_nb \
	--executor_local_memory $memory_size \
	--common_padj_breaks "$padj_breaks" \
	--heatmaps_filter__CHIP "$heatmaps_filters_chip" \
	--heatmaps_filter__motifs "$heatmaps_filters_motifs" \
	--heatmaps_params__peaks_self  "$heatmaps_params__peaks_self" \
	--common__heatmaps_ggplot "$heatmaps_ggplot" \
	--heatmaps_ggplot__chrom_states "$heatmaps_ggplot__chrom_states" \
	--species $species \
	--chromatin_state 'ENCFF321DGG' \
	--chip_ontology 'all'

# chromatin state: ENCFF321DGG -> foreskin fibroblast male newborn from donor(s)


# human run 2: changing the chromatin state and the CHIP ontology
species=human ; run='2' ; cd $case_study_dir/$species
nextflow run jsalignon/cactus \
	-r $cactus_version $latest_flag \
	-profile singularity \
	--references_dir $references_dir \
	--res_dir "${res_prefix}_${species}_run_${run}" \
	--executor_local_cpus $cpu_nb \
	--executor_local_memory $memory_size \
	--common_padj_breaks "$padj_breaks" \
	--heatmaps_filter__CHIP "$heatmaps_filters_chip" \
	--heatmaps_filter__motifs "$heatmaps_filters_motifs" \
	--heatmaps_params__peaks_self  "$heatmaps_params__peaks_self" \
	--common__heatmaps_ggplot "$heatmaps_ggplot" \
	--heatmaps_ggplot__chrom_states "$heatmaps_ggplot__chrom_states" \
	--species $species \
	--chromatin_state 'iHMM.M1K16.human_GM' \
	--chip_ontology 'cell_type.fibroblast'

# chromatin states: GM12878 (iHMM.M1K16.human_GM) -> lymphoblastoid cell line 
#	produced from the blood of a female donor with northern and western 
#	European ancestry by EBV transformation

## Running the case study
#########################

