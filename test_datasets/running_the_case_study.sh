
## notes
# chromatin state: ENCFF321DGG -> foreskin fibroblast male newborn from donor(s)
# chromatin states: HiHMM GM12878 -> lymphoblastoid cell line produced from the blood of a female donor with northern and western European ancestry by EBV transformation

## creating the case study datasets
# homedir=~
# eval homedir=$homedir
# cactus_dir=$homedir/workspace/cactus
# create_test_ds_dir=$cactus_dir/scripts/create_test_datasets
# create_test_datasets_bin_dir=$create_test_ds_dir/bin
# source $create_test_datasets_bin_dir/make_application_note_dirs.sh
# source $create_test_datasets_bin_dir/create_test_datasets_functions.sh


## setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_ds_dir=$cactus_dir/test_datasets
case_study_dir=$test_ds_dir/application_note

## setting up variables
# cactus_version=0.9.0 ; latest_flag=""
# cactus_version=main ; latest_flag="-latest"
cactus_version=57a2e7f516d9ac7114f389c0f226ca01ae63db7a ; latest_flag=""
cpu_nb=47
memory_size='250G'
figshare_version=v4
references_dir=$cactus_dir/references/$figshare_version

padj_breaks="c( 0.2, 1e-5, 1e-50, 1e-100, 1e-200 )"
heatmaps_filters_chip="c( 30, 20, 10, T, 2, 'ward.D', T )"
heatmaps_filters_motifs="c( 30, 20, 10, T, 2, 'ward.D', F )"
heatmaps_ggplot="c( 12, 12, 12 )" 


## running Cactus 
# the parameters after the line break are the shared parameters

# worm run 1: full analysis
species=worm ; cd $case_study_dir/$species
nextflow run jsalignon/cactus --res_dir 'results/2022_12_21__run_1' --chromatin_state 'iHMM.M1K16.worm_L3' --chip_ontology 'all' \
--species $species -profile singularity -r $cactus_version $latest_flag -resume --references_dir $references_dir --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --heatmaps_filter__CHIP "$heatmaps_filters_chip" --heatmaps_filter__motifs "$heatmaps_filters_motifs" --common_padj_breaks "$padj_breaks" --common__heatmaps_ggplot "$heatmaps_ggplot"

# worm run 2: peak-self enrichment only; FDR cutoff 10^-4
species=worm ; cd $case_study_dir/$species
nextflow run jsalignon/cactus --res_dir 'results/2022_12_21__run_2' --chromatin_state 'iHMM.M1K16.worm_L3' --split__threshold_values [4] --do_genes_self_enrichment false --do_chip_enrichment false --do_motif_enrichment false --do_chrom_state_enrichment false --do_func_anno_enrichment false \
--species $species -profile singularity -r $cactus_version $latest_flag -resume --references_dir $references_dir --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --heatmaps_filter__CHIP "$heatmaps_filters_chip" --heatmaps_filter__motifs "$heatmaps_filters_motifs" --common_padj_breaks "$padj_breaks" --common__heatmaps_ggplot "$heatmaps_ggplot"

# human run 1: full analysis
species=human ; cd $case_study_dir/$species
nextflow run jsalignon/cactus --res_dir 'results/2022_12_21__run_1' --chromatin_state 'ENCFF321DGG' --chip_ontology 'all' \
--species $species -profile singularity -r $cactus_version $latest_flag -resume --references_dir $references_dir --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --heatmaps_filter__CHIP "$heatmaps_filters_chip" --heatmaps_filter__motifs "$heatmaps_filters_motifs" --common_padj_breaks "$padj_breaks" --common__heatmaps_ggplot "$heatmaps_ggplot"

# human run 2: changing the chromatin state and the CHIP ontology
species=human ; cd $case_study_dir/$species 
nextflow run jsalignon/cactus --res_dir 'results/2022_12_21__run_2' --chromatin_state 'iHMM.M1K16.human_GM' --chip_ontology 'cell_type.fibroblast' \
--species $species -profile singularity -r $cactus_version $latest_flag -resume --references_dir $references_dir --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --heatmaps_filter__CHIP "$heatmaps_filters_chip" --heatmaps_filter__motifs "$heatmaps_filters_motifs" --common_padj_breaks "$padj_breaks" --common__heatmaps_ggplot "$heatmaps_ggplot"

# human run 3: peak-self enrichment only; FDR cutoff 10^-4
species=human ; cd $case_study_dir/$species
nextflow run jsalignon/cactus --res_dir 'results/2022_12_21__run_3' --chromatin_state 'iHMM.M1K16.human_GM' --split__threshold_values [4] --do_genes_self_enrichment false --do_chip_enrichment false --do_motif_enrichment false --do_chrom_state_enrichment false --do_func_anno_enrichment false \
--species $species -profile singularity -r $cactus_version $latest_flag -resume --references_dir $references_dir --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --heatmaps_filter__CHIP "$heatmaps_filters_chip" --heatmaps_filter__motifs "$heatmaps_filters_motifs" --common_padj_breaks "$padj_breaks" --common__heatmaps_ggplot "$heatmaps_ggplot"



