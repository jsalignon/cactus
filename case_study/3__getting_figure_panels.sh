

## setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
# case_study_dir=$cactus_dir/case_study/data
case_study_dir=$cactus_dir/test_datasets/application_note
case_study_copy_dir=$cactus_dir/case_study/figures/panels

cd $cactus_dir/case_study
mkdir -p $case_study_copy_dir


case_study_dir_human_1=$case_study_dir/human/results/2023_11_26___human_run_1
case_study_dir_human_2=$case_study_dir/human/results/2023_11_26___human_run_2
case_study_dir_worm=$case_study_dir/worm/results/2023_11_26___worm_run_1

case_study_dir_worm_fig=$case_study_dir_worm/Figures_Individual
case_study_dir_worm_pro=$case_study_dir_worm/Processed_Data
case_study_dir_worm_fig_2=$case_study_dir_worm_fig/2_Differential_Analysis
case_study_dir_worm_pro_3=$case_study_dir_worm_pro/3_Enrichment_Analysis

case_study_dir_human_1_fig=$case_study_dir_human_1/Figures_Individual
case_study_dir_human_1_pro=$case_study_dir_human_1/Processed_Data
case_study_dir_human_1_tab_3_kegg=$case_study_dir_human_1/Tables_Individual/3_Enrichment_Analysis/func_anno_KEGG
case_study_dir_human_1_fig_1=$case_study_dir_human_1_fig/1_Preprocessing
case_study_dir_human_1_fig_2=$case_study_dir_human_1_fig/2_Differential_Analysis
case_study_dir_human_1_fig_3=$case_study_dir_human_1_fig/3_Enrichment_Analysis
case_study_dir_human_1_pro_3=$case_study_dir_human_1_pro/3_Enrichment_Analysis

case_study_dir_human_2_pro=$case_study_dir_human_2/Processed_Data
case_study_dir_human_2_pro_3=$case_study_dir_human_2_pro/3_Enrichment_Analysis


# Fig 3: examples of saturation curve, volcano, venns, and PA by filter
cp $case_study_dir_human_1_fig_1/ATAC__multiQC.html $case_study_copy_dir
# wget --force-html --input-file $case_study_copy_dir/ATAC__multiQC.html --output-file $case_study_copy_dir/ATAC__multiQC
cp $case_study_dir_human_1_fig_1/ATAC__peaks__saturation_curve/ctl_1__saturation_curve.pdf $case_study_copy_dir
cp $case_study_dir_human_1_fig_2/ATAC__volcano/ssrp1_vs_ctl__ATAC_volcano.pdf $case_study_copy_dir
cp $case_study_dir_human_1_fig_2/Venn_diagrams__four_ways/all__1.3__ssrp1_vs_ctl__venn_up_and_down.pdf $case_study_copy_dir
cp $case_study_dir_human_1_fig_2/Venn_diagrams__two_ways/all__up__1.3__ssrp1_vs_ctl__venn_up_or_down.pdf $case_study_copy_dir
cp $case_study_dir_human_1_fig_2/ATAC__FDR_by_PA/ssrp1_vs_ctl__ATAC_FDR_by_PA.pdf $case_study_copy_dir

# Fig 4: examples of table, barplot and heatmaps
cp $case_study_dir_human_1_tab_3_kegg/ATAC__all__down__1.3__ssrp1_vs_ctl__func_anno_KEGG__enrich.xlsx $case_study_copy_dir
cp $case_study_dir_human_1_fig_3/Barplots__func_anno_KEGG/ATAC__all__down__1.3__ssrp1_vs_ctl__func_anno_KEGG__barplot.pdf $case_study_copy_dir
cp $case_study_dir_human_1_fig_3/Heatmaps__peaks_self/ATAC__all__1.3__all__peaks_self__heatmap.pdf $case_study_copy_dir
cp $case_study_dir_human_1_fig_3/Heatmaps__peaks_self/ATAC__all__1.3__supt16h__peaks_self__heatmap.pdf $case_study_copy_dir

# Fig 5: peak_self heatmaps
cp $case_study_dir_worm_pro_3/Heatmaps__peaks_self/mRNA__Null__1.3__ctl__peaks_self__heatmap.rds $case_study_copy_dir/mRNA__Null__1.3__ctl__peaks_self__heatmap__worm.rds
cp $case_study_dir_human_1_pro_3/Heatmaps__peaks_self/mRNA__Null__1.3__ctl__peaks_self__heatmap.rds $case_study_copy_dir/mRNA__Null__1.3__ctl__peaks_self__heatmap__human.rds
cp $case_study_dir_worm_pro_3/Heatmaps__peaks_self/ATAC__all__1.3__ctl__peaks_self__heatmap.rds $case_study_copy_dir/ATAC__all__1.3__ctl__peaks_self__heatmap__worm.rds
cp $case_study_dir_human_1_pro_3/Heatmaps__peaks_self/ATAC__all__1.3__ctl__peaks_self__heatmap.rds $case_study_copy_dir/ATAC__all__1.3__ctl__peaks_self__heatmap__human.rds

# Fig 6: CHIP and motifs heatmaps
cp $case_study_dir_worm_pro_3/Heatmaps__motifs/ATAC__all__1.3__ctl__motifs__heatmap.rds $case_study_copy_dir/ATAC__all__1.3__ctl__motifs__heatmap__worm.rds
cp $case_study_dir_worm_pro_3/Heatmaps__CHIP/ATAC__all__1.3__ctl__CHIP__heatmap.rds $case_study_copy_dir/ATAC__all__1.3__ctl__CHIP__heatmap__worm.rds
cp $case_study_dir_human_1_pro_3/Heatmaps__motifs/ATAC__all__1.3__ctl__motifs__heatmap.rds $case_study_copy_dir/ATAC__all__1.3__ctl__motifs__heatmap__human.rds
cp $case_study_dir_human_1_pro_3/Heatmaps__CHIP/ATAC__all__1.3__ctl__CHIP__heatmap.rds $case_study_copy_dir/ATAC__all__1.3__ctl__CHIP__heatmap__human.rds

# Fig 7: CHIP both_atac and both_mrna heatmaps
cp $case_study_dir_human_1_pro_3/Heatmaps__CHIP/both_ATAC__all__1.3__ctl__CHIP__heatmap.rds $case_study_copy_dir
cp $case_study_dir_human_1_pro_3/Heatmaps__CHIP/both_mRNA__all__1.3__ctl__CHIP__heatmap.rds $case_study_copy_dir

# Fig 8: chrom states heatmaps
cp $case_study_dir_worm_pro_3/Heatmaps__chrom_states/ATAC__all__1.3__ctl__chrom_states__heatmap.rds $case_study_copy_dir/ATAC__all__1.3__ctl__chrom_states__heatmap__worm.rds
cp $case_study_dir_human_1_pro_3/Heatmaps__chrom_states/ATAC__all__1.3__ctl__chrom_states__heatmap.rds $case_study_copy_dir/ATAC__all__1.3__ctl__chrom_states__heatmap__human_1.rds
cp $case_study_dir_human_2_pro_3/Heatmaps__chrom_states/ATAC__all__1.3__ctl__chrom_states__heatmap.rds $case_study_copy_dir/ATAC__all__1.3__ctl__chrom_states__heatmap__human_2.rds

# Tab S1: reprogramming genes defined in Figure 6L of the Kolundzic et al. study
cp $case_study_dir_human_1/Tables_Merged/2_Differential_Analysis/res_simple.xlsx $case_study_copy_dir

# Fig S2: HA-HE gene tracks
# find "$case_study_dir/worm/results/2024_04_09___worm_run_1/Processed_Data" -regex ".*\(hmg4\|ctl\).*.bw" -exec cp "{}" $case_study_copy_dir \;
find "$case_study_dir/worm/results/2024_04_09___worm_run_1/Processed_Data" -regex ".*\(hmg4\|ctl\).*.bw" -exec bash -c 'cp "$0" "${1}/$(basename "${0%.*}").bigWig"' {} "$case_study_copy_dir" \;
cp "$case_study_dir/worm/results/2024_04_09___worm_run_1/Tables_Individual/2_Differential_Analysis/res_filter/hmg4_vs_ctl__res_filter.xlsx" $case_study_copy_dir
cp "$case_study_dir/worm/results/2024_04_09___worm_run_1/Processed_Data/2_Differential_Analysis/DA_split__bed_regions/both_ATAC__all__up__1.3__hmg4_vs_ctl__regions.bed" $case_study_copy_dir
cp "$case_study_dir/worm/results/2024_04_09___worm_run_1/Processed_Data/2_Differential_Analysis/DA_split__bed_regions/both_ATAC__all__down__1.3__hmg4_vs_ctl__regions.bed" $case_study_copy_dir


# converting the heavy plots to png
pdftoppm -rx 300 -ry 300 -png \
	$case_study_copy_dir/ssrp1_vs_ctl__ATAC_volcano.pdf \
	$case_study_copy_dir/ssrp1_vs_ctl__ATAC_volcano
pdftoppm -rx 300 -ry 300 -png \
	$case_study_copy_dir/ssrp1_vs_ctl__ATAC_FDR_by_PA.pdf \
	$case_study_copy_dir/ssrp1_vs_ctl__ATAC_FDR_by_PA

# note that panel 4a and 5a needs to be created manually by exporting figure
# from multiQC for 4a and by creating a screenshot for 5a

