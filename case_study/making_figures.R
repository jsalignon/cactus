
# homedir=~
# eval homedir=$homedir
# cactus_dir=$homedir/workspace/cactus # to change if needed
# cd $cactus_dir/case_study/figures


R


# loading functions and libraries

readRDS1 <- function(x) readRDS(paste0('panels/', x))
pdf_a4   <- function(x, ...) pdf(x, paper = 'a4r', ...)
library(ggplot2)
library(magrittr)

# Fig 6: peak_self heatmaps

p1 = readRDS1('mRNA__Null__1.3__ctl__peaks_self__heatmap__worm.rds') + ggtitle('worm, DEGs') 
p2 = readRDS1('mRNA__Null__1.3__ctl__peaks_self__heatmap__human.rds') + ggtitle('human, DEGs') 
p3 = readRDS1('ATAC__all__1.3__ctl__peaks_self__heatmap__worm.rds') + ggtitle('worm, DARs') 
p4 = readRDS1('ATAC__all__1.3__ctl__peaks_self__heatmap__human.rds') + ggtitle('human, DARs') 

pp = ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = T, 
	legend = 'right', labels = letters[1:4])

pdf_a4('Figure_6.pdf')
	print(pp)
dev.off()


# Fig 7: CHIP and motifs heatmaps

p1 = readRDS1('ATAC__all__1.3__ctl__motifs__heatmap__worm.rds') + ggtitle('worm, motifs') 
p2 = readRDS1('ATAC__all__1.3__ctl__CHIP__heatmap__worm.rds') + ggtitle('worm, CHIP') 
p3 = readRDS1('ATAC__all__1.3__ctl__motifs__heatmap__human.rds') + ggtitle('human, motifs') 
p4 = readRDS1('ATAC__all__1.3__ctl__CHIP__heatmap__human.rds') + ggtitle('human, CHIP') 

pp = ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = T, 
	legend = 'right', labels = letters[1:4])

pdf_a4('Figure_7.pdf')
	print(pp)
dev.off()


# Fig 8: CHIP both_atac and both_mrna heatmaps

p1 = readRDS1('both_ATAC__all__1.3__ctl__CHIP__heatmap.rds') + ggtitle('human, DARs, CHIP,\n co-regulated genes')
p2 = readRDS1('both_mRNA__all__1.3__ctl__CHIP__heatmap.rds') + ggtitle('human, promoters, CHIP\n co-regulated genes')
p1 = p1 + theme(plot.title = element_text(hjust = 1))
p2 = p2 + theme(plot.title = element_text(hjust = 1))

pp = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = T, 
	legend = 'right', labels = letters[1:2])

pdf_a4('Figure_8.pdf')
	print(pp)
dev.off()


# Fig 9: chrom states heatmaps

p1 = readRDS1('ATAC__all__1.3__ctl__chrom_states__heatmap__worm.rds') + ggtitle('HiHMM, worm\nL3 larvae')
p2 = readRDS1('ATAC__all__1.3__ctl__chrom_states__heatmap__human_2.rds') + ggtitle('ChromHMM, human\nfibroblast cell line')
p3 = readRDS1('ATAC__all__1.3__ctl__chrom_states__heatmap__human_1.rds') + ggtitle('HiHMM, human\nlymphoblastoid cell line')

update_text_size <- function(px) {
	px + theme(
		axis.text.y = element_text(size = 8), 
		axis.text.x = element_text(size = 12), 
		plot.title = element_text(size = 10, hjust = 1)
		)
}

p1 %<>% update_text_size 
p2 %<>% update_text_size 
p3 %<>% update_text_size 

p1 = p1 + theme(
	legend.text = element_text(size = 6), 
	legend.title = element_text(size = 6),
	legend.key.size = unit(0.5, 'cm'))

p1$scales$scales[[2]]$labels %<>% gsub('Heterochromatin', 'Heterochrom.', .)
p3$scales$scales[[2]]$labels %<>% gsub('Repressed PolyComb', 'PC repressed', .) %>% 
	gsub('Weak Repressed PolyComb', 'PC repressed wk', .) %>%
	gsub('Downstream', 'Down.', .) %>% gsub('Upstream', 'Up.', .) %>%
	gsub('Bivalent', 'Biv.', .) %>% gsub('repeats', 'rep.', .)

p1 = p1 + coord_fixed(1.3)

pp = ggpubr::ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = T, 
	legend = 'right', labels = letters[1:3], label.y = 0.85)

pdf_a4('Figure_9.pdf')
	print(pp)
dev.off()




