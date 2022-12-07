
# # creating the case study datasets
# homedir=~
# eval homedir=$homedir
# cactus_dir=$homedir/workspace/cactus
# create_test_ds_dir=$cactus_dir/scripts/create_test_datasets
# create_test_datasets_bin_dir=$create_test_ds_dir/bin
# source $create_test_datasets_bin_dir/make_application_note_dirs.sh
# source $create_test_datasets_bin_dir/create_test_datasets_functions.sh


# setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_ds_dir=$cactus_dir/test_datasets
case_study_dir=$test_ds_dir/application_note

# setting up variables
cactus_version=0.9.0
cpu_nb=47
memory_size='250G'



padj_breaks="c( 0.2, 1e-5, 1e-50, 1e-100, 1e-200 )"
heatmaps_filters_chip="c( 30, 20, 10, T, 2, 'ward.D', T )"
heatmaps_filters_motifs="c( 30, 20, 10, T, 2, 'ward.D', F )"

cd $case_study_dir/worm
nextflow run jsalignon/cactus -profile singularity -r $cactus_version -latest -resume --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --species 'worm' --res_dir 'results/2022_12_02' --chromatin_state 'iHMM.M1K16.worm_L3' --chip_ontology 'all' --split__threshold_type 'FDR' --split__threshold_values [1.3] --split__peak_assignment ['all'] --heatmaps_filter__CHIP "$heatmaps_filters_chip" --heatmaps_filter__motifs "$heatmaps_filters_motifs" --padj_breaks__motifs "$padj_breaks" --padj_breaks__CHIP "$padj_breaks" --padj_breaks__chrom_states "$padj_breaks"

cd $case_study_dir/human

# CHIP ontology: 'all'
# chromatin state: ENCFF321DGG -> foreskin fibroblast male newborn from donor(s)
nextflow run jsalignon/cactus -profile singularity -r $cactus_version -latest -resume --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --species 'human' --res_dir 'results/2022_12_02__run_1' --chromatin_state 'ENCFF321DGG' --chip_ontology 'all' --split__threshold_type 'FDR' --split__threshold_values [1.3] --split__peak_assignment ['all'] --heatmaps_filter__CHIP "$heatmaps_filters_chip" --heatmaps_filter__motifs "$heatmaps_filters_motifs" --padj_breaks__motifs "$padj_breaks" --padj_breaks__CHIP "$padj_breaks" --padj_breaks__chrom_states "$padj_breaks"

# CHIP ontology: 'cell_type.fibroblast'
# chromatin states: HiHMM GM12878 -> lymphoblastoid cell line produced from the blood of a female donor with northern and western European ancestry by EBV transformation
nextflow run jsalignon/cactus -profile singularity -r $cactus_version -latest -resume --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --species 'human' --res_dir 'results/2022_12_02__run_2' --chromatin_state 'iHMM.M1K16.human_GM' --chip_ontology 'cell_type.fibroblast' --split__threshold_type 'FDR' --split__threshold_values [1.3] --split__peak_assignment ['all'] --heatmaps_filter__CHIP "$heatmaps_filters_chip" --heatmaps_filter__motifs "$heatmaps_filters_motifs" --padj_breaks__motifs "$padj_breaks" --padj_breaks__CHIP "$padj_breaks" --padj_breaks__chrom_states "$padj_breaks"


