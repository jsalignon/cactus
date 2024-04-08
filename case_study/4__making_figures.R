
# homedir=~
# eval homedir=$homedir
# cactus_dir=$homedir/workspace/cactus
# cd $cactus_dir/case_study/figures


R


# loading functions and libraries

library(ggplot2)
library(magrittr)
library(data.table)

readRDS1 <- function(x) readRDS(paste0('panels/', x))
pdf_a4   <- function(x, ...) pdf(x, paper = 'a4r', ...)
update_text_size <- function(px) {
	px + theme(
		axis.text.y = element_text(size = 8), 
		axis.text.x = element_text(size = 10), 
		plot.title = element_text(size = 10, hjust = 1)
		)
}

# Fig 5: peak_self heatmaps

p1 = readRDS1('mRNA__Null__1.3__ctl__peaks_self__heatmap__worm.rds') + ggtitle('worm, DEGs') 
p2 = readRDS1('mRNA__Null__1.3__ctl__peaks_self__heatmap__human.rds') + ggtitle('human, DEGs') 
p3 = readRDS1('ATAC__all__1.3__ctl__peaks_self__heatmap__worm.rds') + ggtitle('worm, DARs') 
p4 = readRDS1('ATAC__all__1.3__ctl__peaks_self__heatmap__human.rds') + ggtitle('human, DARs') 

pp = ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = T, 
	legend = 'right', labels = letters[1:4])

pdf('Figure_5.pdf', width = 7.5, height = 7.5)
	print(pp)
dev.off()


# Fig 6: CHIP and motifs heatmaps

p1 = readRDS1('ATAC__all__1.3__ctl__motifs__heatmap__worm.rds') + ggtitle('worm, motifs') 
p2 = readRDS1('ATAC__all__1.3__ctl__CHIP__heatmap__worm.rds') + ggtitle('worm, CHIP') 
p3 = readRDS1('ATAC__all__1.3__ctl__motifs__heatmap__human.rds') + ggtitle('human, motifs') 
p4 = readRDS1('ATAC__all__1.3__ctl__CHIP__heatmap__human.rds') + ggtitle('human, CHIP') 

p1 %<>% update_text_size 
p2 %<>% update_text_size 
p3 %<>% update_text_size 
p4 %<>% update_text_size 

p1 = p1 + theme(
	legend.text = element_text(size = 6), 
	legend.title = element_text(size = 6),
	legend.key.size = unit(0.5, 'cm'))

pp = ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = T, 
	legend = 'right', labels = letters[1:4], widths = c(1.19, 1.3, 1.02, 1.4))

pdf('Figure_6.pdf', width = 7, height = 5.7)
	print(pp)
dev.off()



# Fig 7: CHIP both_atac and both_mrna heatmaps

p1 = readRDS1('both_ATAC__all__1.3__ctl__CHIP__heatmap.rds') + ggtitle('human, DARs, CHIP,\n co-regulated genes')
p2 = readRDS1('both_mRNA__all__1.3__ctl__CHIP__heatmap.rds') + ggtitle('human, promoters, CHIP\n co-regulated genes')
p1 = p1 + theme(plot.title = element_text(hjust = 1))
p2 = p2 + theme(plot.title = element_text(hjust = 1))

pp = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = T, 
	legend = 'right', labels = letters[1:2])

pdf('Figure_7.pdf', width = 7, height = 7)
	print(pp)
dev.off()


# Fig 8: chrom states heatmaps

p1 = readRDS1('ATAC__all__1.3__ctl__chrom_states__heatmap__worm.rds') + ggtitle('HiHMM, worm\nL3 larvae')
p2 = readRDS1('ATAC__all__1.3__ctl__chrom_states__heatmap__human_2.rds') + ggtitle('ChromHMM, human\nfibroblast cell line')
p3 = readRDS1('ATAC__all__1.3__ctl__chrom_states__heatmap__human_1.rds') + ggtitle('HiHMM, human\nlymphoblastoid cell line')

p1 %<>% update_text_size 
p2 %<>% update_text_size 
p3 %<>% update_text_size 

p1 = p1 + theme(
	legend.text = element_text(size = 6), 
	legend.title = element_text(size = 6),
	legend.key.size = unit(0.5, 'cm')
	)

pp = ggpubr::ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = T, 
	legend = 'right', labels = letters[1:3], widths = c(1.2, 1.03, 1.27))

pdf('Figure_8.pdf', width = 7, height = 3.5)
	print(pp)
dev.off()



# Tab S1: reprogramming genes defined in Figure 6L of the Kolundzic et al. study

dt = openxlsx::read.xlsx('panels/res_simple.xlsx') %>% setDT
dt = dt[ET == 'mRNA' & COMP %in% c('ssrp1_vs_ctl', 'supt16h_vs_ctl')]

# reprogramming-promoting genes
dt_up = dt[gene_name %in% c('FGF2', 'CEBPB', 'ESRRB', 'PRDM15', 'BMP2', 'JMJD2C', 'SALL4', 'LIN28B', 'KDM4C')]
dt_up = dt_up[order(-L2FC), c('COMP', 'gene_name', 'gene_id', 'pval', 'padj', 'L2FC', 'FC')]

# reprogramming-preventing genes
dt_dw = dt[gene_name %in% c('SUV39H1', 'SUV39H2', 'PTPN11', 'NR2F1', 'SUMO2', 'CHAF1B', 'RBBP7', 'DRAM1')]
dt_dw = dt_dw[order(-L2FC), c('COMP', 'gene_name', 'gene_id', 'pval', 'padj', 'L2FC', 'FC')]

openxlsx::write.xlsx(list(up = dt_up, down = dt_dw), 'Table_S1.xlsx')



homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
cd $cactus_dir/case_study/figures

# ~/workspace/cactus/docs/examples/png
# /raid/jersal/workspace/cactus/docs/examples/png

library(ggplot2)

png_to_gg <- function(png_file, margin = 20) {
	ggplot2::ggplot() + ggplot2::annotation_custom(
		grid::rasterGrob(png::readPNG(png_file), 
			width = ggplot2::unit(1,"npc"), 
			height = ggplot2::unit(1,"npc")), 
		-Inf, Inf, -Inf, Inf) + theme(plot.margin = 
    margin(t = margin, r = margin, b = margin, l = margin, unit = "pt"))
}

grep1 <- function(x) grep(x, v_png_files, value = T)

v_png_files = c(list.files('../../docs/examples/png', full.names = T), 
				list.files('../../docs/examples/xlsx_png', full.names = T))

v_files_1_QC_reads = grep1('\\/(pca|spearman|ctl_1__(average|reads|insert))')
v_files_1_QC_peaks = grep1('(\\/(ctl_1__peaks|ATAC__peaks)|saturation_curve)')
v_files_1_QC       = c(v_files_1_QC_reads, v_files_1_QC_peaks)
v_files_2_DA_fig   = grep1('\\/hmg4_vs_ctl__') %>% c(grep1('venn'))
# v_files_2_DA_fig   = grep1('\\/hmg4_vs_ctl__') %>%
# 	grep('mRNA_(PCA|volcano)', ., invert = T, value = T) %>% c(grep1('venn'))
v_files_2_DA_tab   = c(grep1('xlsx_png\\/ATAC'), grep1('xlsx_png\\/mRNA'), 
	grep1('xlsx_png\\/res'))
v_files_2_DA       = c(v_files_2_DA_fig, v_files_2_DA_tab)
v_files_3_EN_fig   = c(grep1('__barplot\\.'), grep1('__heatmap\\.'))
v_files_3_EN_tab   = grep1('xlsx_png\\/(CHIP|motifs|func|peaks|genes|chrom)')
v_files_3_EN       = c(v_files_3_EN_fig, v_files_3_EN_tab)

length(v_files_1_QC_reads)
length(v_files_1_QC_peaks)
length(v_files_1_QC) 
length(v_files_2_DA_fig)
length(v_files_2_DA_tab)
length(v_files_3_EN_fig)
length(v_files_3_EN_tab)

plot_nrow_by_ncol <- function(v_files, nrow = 4, ncol = 3, v_labels = 'auto'){
	lp = lapply(v_files, png_to_gg)
	p1 = ggpubr::ggarrange(plotlist = lp, nrow = nrow, ncol = ncol, 
			labels = v_labels)
	return(p1)
}

plot_figures <- function(v_files, nrow = 4, ncol = 3, v_labels = 'auto'){
	plot_nrow_by_ncol(v_files, nrow = nrow, ncol = ncol, v_labels = v_labels)
}

plot_tables <- function(v_files, nrow = 6, ncol = 2, v_labels = v_labels){
	plot_nrow_by_ncol(v_files, nrow = nrow, ncol = ncol, v_labels = v_labels)
}

# plot_1_by_3 <- function(v_files_tab, v_files_fig, nrow = 5, ncol = 1){
# 	len      = length(v_files_2_DA_fig)
# 	v_labels = letters[len:(len + length(v_files_2_DA_tab))]
# 	plot_nrow_by_ncol(v_files_tab, nrow = nrow, ncol = ncol, 
# 		v_labels = v_labels)
# }

p_QC = plot_4_by_3(v_files_1_QC)
pdf('Figure_S3.pdf', width = 7, height = 10)
	print(p_QC)
dev.off()

length(v_files_2_DA_fig) # 17
length(v_files_2_DA_fig) - 12 # 5
length(v_files_2_DA_tab) #  7
p_DA_1   = plot_figures(v_files_2_DA_fig[1:12])
p_DA_2_1 = plot_figures(v_files_2_DA_fig[13:17], nrow = 2,
	v_labels = letters[12 + (1:6)])
p_DA_2_2 = plot_tables(v_files_2_DA_tab, v_labels = letters[12 + (6:13)], 
	nrow = 4)
p_DA_2 = ggpubr::ggarrange(p_DA_2_1, p_DA_2_2, nrow = 2, ncol = 1)
pdf('Figure_S4.pdf', width = 7, height = 10)
	print(p_DA_1)
	print(p_DA_2)
dev.off()


length(v_files_3_EN_fig) # 14
length(v_files_3_EN_fig) - 12 # 2
length(v_files_3_EN_tab) #  7
p_EN_1   = plot_figures(v_files_3_EN_fig[1:12])
p_EN_2_1 = plot_figures(v_files_3_EN_fig[13:14], nrow = 1,
	v_labels = letters[12 + (1:2)])
p_EN_2_2 = plot_tables(v_files_3_EN_tab, 
	v_labels = letters[12 + (3:10)], nrow = 5)
p_EN_2 = ggpubr::ggarrange(p_EN_2_1, p_EN_2_2, nrow = 2, ncol = 1, 
	heights = c(0.25, 0.75) )
pdf('Figure_S5.pdf', width = 7, height = 10)
	print(p_EN_1)
	print(p_EN_2)
dev.off()



v_files_3_EN_fig
v_files_3_EN_tab
plot_4_by_3(v_files_1_QC)



plot_1_by_3 <- function(v_files_tab, v_labels, nrow = 1, ncol = 3){


plot_1_by_3 <- function(v_files_tab, v_labels, nrow = 1, ncol = 3){
	lp = lapply(v_files_tab, png_to_gg)
	p1 = ggpubr::ggarrange(plotlist = lp, nrow = nrow, ncol = ncol, 
		labels = v_labels)
	return(p1)
}

p1 = plot_4_by_3(v_files_1_1_QC)
p2 = plot_4_by_3(v_files_1_1_QC)


v_files_1_2_DA %>% length
v_files_1_3_EN %>% length


p1 = plot_4_by_3(v_files_1_1_QC)
p2 = plot_4_by_3(v_files_1_2_DA)
p3 = plot_4_by_3(v_files_1_3_EN)



list.files('../../docs/examples/png')
list.files('../../docs/examples/')


v_files_1 = c(
	grep('\\/ctl_1__', v_files, value = T),
	grep('without_control_pca', v_files, value = T),
	grep('spearman_correlation_heatmap', v_files, value = T),
	grep('ATAC__peaks__grouped', v_files, value = T),
	grep('hmg4_vs_ctl__ATAC_', v_files, value = T),
	grep('hmg4_vs_ctl__mRNA_', v_files, value = T),
	# grep('mRNA_volcano__', v_files, value = T),
	grep('venn_up', v_files, value = T),
	grep('__barplot', v_files, value = T),
	grep('__heatmap\\.', v_files, value = T)
	)



v_files_1_DA_mrna_ = grep1('\\/hmg4_vs_ctl__mRNA')

plot_4_by_3 <- function(v_files_1){
	lp = lapply(v_files_1, png_to_gg)
	p1 = ggpubr::ggarrange(plotlist = lp, nrow = 4, ncol = 3, 
			labels = 'auto')
	return(p1)
}



plot_by_chunk <- function(lp, nrow, ncol){
	nplots = length(lp)
	npages = ceiling(nplots / (nrow * ncol))
	lp1 = vector('list', length = npages)

	for(page in 1:npages){
		cur_panel = 1 + (page - 1) * 12
		v_panels = cur_panel:min(cur_panel + 11, nplots)
		p1 = ggpubr::ggarrange(plotlist = lp[v_panels], nrow = 4, ncol = 3, 
			labels = v_panels)
		lp1[[page]] = p1
	}

	return(lp1)
}






rep('__barplot', v_files, value = T),
	grep('__heatmap\\.'

v_files_1_DA_mrna_ = grep1('\\/hmg4_vs_ctl__mRNA')

plot_4_by_3 <- function(v_files_1){
	lp = lapply(v_files_1, png_to_gg)
	p1 = ggpubr::ggarrange(plotlist = lp, nrow = 4, ncol = 3, 
			labels = 'auto')
	return(p1)
}



 = ggpubr::ggarrange(plotlist = lp[v_panels], nrow = 4, ncol = 3, 
	labels = v_panels)




grep('mRNA_(PCA|volcano)', grep1('\\/hmg4_vs_ctl__'), invert = T, value = T)


v_files_1_1_QC = c(v_files_1_1_QC_reads, v_files_1_1_QC_peaks)


_1

grep1('__(barplot|heatmap)\\.')


v_files_2 = c(v_files_1_1_QC_reads, v_files_1_1_QC_peaks, v_files_1_2_DA, v_files_1_3_enrich)

v_files %>% .[!. %in% v_files_2]

> v_files %>% .[!. %in% v_files_2]
[1] "../../docs/examples/png/all__1000__hmg4_vs_ctl__venn_up_and_down.png"
[2] "../../docs/examples/png/all__down__1000__hmg4_vs_ctl__venn_up_or_down.png"
[3] "../../docs/examples/png/all__up__1000__hmg4_vs_ctl__venn_up_or_down.png"
[4] "../../docs/examples/png/ctl_1__saturation_curve.png"
[5] "../../docs/examples/png/hmg4_vs_spt16__ATAC_volcano__no_rtr.png"
[6] "../../docs/examples/png/hmg4_vs_spt16__ATAC_volcano.png"
[7] "../../docs/examples/png/mRNA_volcano__no_gtr.png"
[8] "../../docs/examples/png/mRNA_volcano__with_gtr.png"
>


	pca|spearman|ctl_1__reads)')



grep1('(pca|spearman|\\/ctl_1__(average|reads|insert))')

v_files_1_QC_reads = grep('(pca|spearman|\\/ctl_1__(average|reads|insert))', 
	v_files, value = T),


 [2] "../../docs/examples/png/ctl_1__reads_coverage.png"
 [3] "../../docs/examples/png/ctl_1__peaks_coverage.png"
 [4] "../../docs/examples/png/ctl_1__reads_coverage.png"
 [5] "../../docs/examples/png/ctl_1__saturation_curve.png"
 [6] "../../docs/examples/png/pca_top5000_without_control_pca.png"
 [7] "../../docs/examples/png/spearman_correlation_heatmap_without_outliers_without_control_cor.png"
 [8] "../../docs/examples/png/ATAC__peaks__grouped__annotation_barplot.png"
 [9] "../../docs/examples/png/ATAC__peaks__grouped__average_profile.png"
[10] "../../docs/examples/png/ATAC__peaks__grouped__distance_to_TSS.png"
[11] "../../docs/examples/png/hmg4_vs_ctl__ATAC_FDR_by_PA.png"
[12] "../../docs/examples/png/hmg4_vs_ctl__ATAC_other_plots-1.png"
[13] "../../docs/examples/png/hmg4_vs_ctl__ATAC_other_plots-2.png"
[14] "../../docs/examples/png/hmg4_vs_ctl__ATAC_other_plots-3.png"
[15] "../../docs/examples/png/hmg4_vs_ctl__ATAC_PCA_1_2.png"
[16] "../../docs/examples/png/hmg4_vs_ctl__ATAC_PCA_3_4.png"
[17] "../../docs/examples/png/hmg4_vs_ctl__ATAC_volcano.png"
[18] "../../docs/examples/png/hmg4_vs_ctl__mRNA_other_plots-1.png"
[19] "../../docs/examples/png/hmg4_vs_ctl__mRNA_other_plots-2.png"
[20] "../../docs/examples/png/hmg4_vs_ctl__mRNA_PCA_1_2.png"
[21] "../../docs/examples/png/hmg4_vs_ctl__mRNA_PCA_3_4.png"
[22] "../../docs/examples/png/hmg4_vs_ctl__mRNA_volcano.png"
[23] "../../docs/examples/png/all__1000__hmg4_vs_ctl__venn_up_and_down.png"
[24] "../../docs/examples/png/all__down__1000__hmg4_vs_ctl__venn_up_or_down.png"
[25] "../../docs/examples/png/all__up__1000__hmg4_vs_ctl__venn_up_or_down.png"
[26] "../../docs/examples/png/ATAC__all__down__1000__hmg4_vs_ctl__CHIP__barplot.png"
[27] "../../docs/examples/png/ATAC__all__down__1000__hmg4_vs_ctl__chrom_states__barplot.png"
[28] "../../docs/examples/png/ATAC__all__down__1000__hmg4_vs_ctl__func_anno_BP__barplot.png"
[29] "../../docs/examples/png/ATAC__all__down__1000__hmg4_vs_ctl__func_anno_KEGG__barplot.png"
[30] "../../docs/examples/png/ATAC__all__down__1000__hmg4_vs_ctl__genes_self__barplot.png"
[31] "../../docs/examples/png/ATAC__all__down__1000__hmg4_vs_ctl__peaks_self__barplot.png"
[32] "../../docs/examples/png/ATAC__all__down__1000__hmg4_vs_spt16__motifs__barplot.png"
[33] "../../docs/examples/png/ATAC__all__1000__all__CHIP__heatmap.png"
[34] "../../docs/examples/png/ATAC__all__1000__all__chrom_states__heatmap.png"
[35] "../../docs/examples/png/ATAC__all__1000__all__func_anno_BP__heatmap.png"
[36] "../../docs/examples/png/ATAC__all__1000__all__func_anno_KEGG__heatmap.png"
[37] "../../docs/examples/png/ATAC__all__1000__all__genes_self__heatmap.png"
[38] "../../docs/examples/png/ATAC__all__1000__all__motifs__heatmap.png"
[39] "../../docs/examples/png/ATAC__all__1000__all__peaks_self__heatmap.png"
>'



v_files = 

list.files('../../docs/examples/png', full.names = T)
list.files('../../docs/examples/png', full.names = T)



v_files = list.files('../../docs/examples/png', full.names = T)

v_files_QC_reads = c(
	grep('\\/ctl_1__', v_files, value = T),
	grep('without_control_pca', v_files, value = T),
	grep('spearman_correlation_heatmap', v_files, value = T)
)

v_files_1 = c(
	grep('\\/ctl_1__', v_files, value = T),
	grep('without_control_pca', v_files, value = T),
	grep('spearman_correlation_heatmap', v_files, value = T),
	grep('ATAC__peaks__grouped', v_files, value = T),
	grep('hmg4_vs_ctl__ATAC_', v_files, value = T),
	grep('hmg4_vs_ctl__mRNA_', v_files, value = T),
	# grep('mRNA_volcano__', v_files, value = T),
	grep('venn_up', v_files, value = T),
	grep('__barplot', v_files, value = T),
	grep('__heatmap\\.', v_files, value = T)
	)
v_files %>% .[!. %in% v_files_1]
# [1] "../../docs/examples/png/hmg4_vs_spt16__ATAC_volcano__no_rtr.png"
# [2] "../../docs/examples/png/hmg4_vs_spt16__ATAC_volcano.png"
# [3] "../../docs/examples/png/mRNA_volcano__no_gtr.png"
# [4] "../../docs/examples/png/mRNA_volcano__with_gtr.png"

v_files_1 %>% .[which(duplicated(.))] # character(0)

v_files_1 %<>% .[. != '../../docs/examples/png/hmg4_vs_ctl__ATAC_other_plots.png']
v_files_1 %<>% .[. != '../../docs/examples/png/hmg4_vs_ctl__mRNA_other_plots.png']

lp = lapply(v_files_1, png_to_gg)

lp1 = plot_by_chunk(lp, 4, 3)

pdf('Sup_fig_test.pdf', width = 7, height = 10)
	print(lp1)
dev.off()





# lp = lapply(v_files_1, png_to_gg)
# p1 = ggpubr::ggarrange(plotlist = lp, nrow = 4, ncol = 3, 
# 	labels = 'auto')

# pdf('Sup_fig_test.pdf', width = 7, height = 10)
# 	print(p1)
# dev.off()

