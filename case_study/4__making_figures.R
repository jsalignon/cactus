
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





# Supplementary Figures showing all outputs from Cactus

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
v_files_2_DA_tab   = c(grep1('xlsx_png\\/ATAC'), grep1('xlsx_png\\/mRNA'), 
	grep1('xlsx_png\\/res'))
v_files_2_DA       = c(v_files_2_DA_fig, v_files_2_DA_tab)
v_files_3_EN_fig   = c(grep1('__barplot\\.'), grep1('__heatmap\\.'))
v_files_3_EN_tab   = grep1('xlsx_png\\/(CHIP|motifs|func|peaks|genes|chrom)')
v_files_3_EN       = c(v_files_3_EN_fig, v_files_3_EN_tab)

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



