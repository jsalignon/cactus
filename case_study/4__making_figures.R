
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

plot_nrow_by_ncol <- function(v_files, nrow = 4, ncol = 3, v_labels = 'auto',
	margin = 10){
	lp = lapply(v_files, png_to_gg, margin = margin)
	p1 = ggpubr::ggarrange(plotlist = lp, nrow = nrow, ncol = ncol, 
			labels = v_labels)
	return(p1)
}

plot_figures <- function(v_files, nrow = 4, ncol = 3, v_labels = 'auto'){
	plot_nrow_by_ncol(v_files, nrow = nrow, ncol = ncol, v_labels = v_labels)
}

plot_tables <- function(v_files, nrow = 6, ncol = 2, v_labels = v_labels,
	margin = 15){
	plot_nrow_by_ncol(v_files, nrow = nrow, ncol = ncol, v_labels = v_labels,
		margin = margin)
}

v_png_files = c(list.files('../../docs/examples/png', full.names = T), 
				list.files('../../docs/examples/xlsx_png', full.names = T))

v_files_1_QC_reads = grep1('\\/(pca|spearman|ctl_1__(reads|insert))')
v_files_1_QC_peaks = grep1('(\\/(ctl_1__(peaks|average)|ATAC__peaks)|saturation_curve)')
v_files_1_QC       = c(v_files_1_QC_reads, v_files_1_QC_peaks)
p_QC = plot_figures(v_files_1_QC)
pdf('Figure_S2.pdf', width = 7, height = 10)
	print(p_QC)
dev.off()


v_files_2_DA_fig   = grep1('\\/hmg4_vs_ctl__') %>% c(grep1('venn'))
v_files_2_DA_fig %<>% .[. != '../../docs/examples/png/hmg4_vs_ctl__ATAC_other_plots.png']
v_files_2_DA_fig %<>% .[. != '../../docs/examples/png/hmg4_vs_ctl__mRNA_other_plots.png']
v_files_2_DA_fig %<>% sort %>% rev
v_files_2_DA_tab   = c(grep1('xlsx_png\\/ATAC'), grep1('xlsx_png\\/mRNA'), 
	grep1('xlsx_png\\/res'))
length(v_files_2_DA_fig) # 15
length(v_files_2_DA_fig) - 12 # 3
length(v_files_2_DA_tab) #  7
length(v_files_2_DA_tab) + length(v_files_2_DA_fig) - 12 #  10
p_DA_1   = plot_figures(v_files_2_DA_fig[1:12])
p_DA_2_1 = plot_figures(v_files_2_DA_fig[13:15], nrow = 1,
	v_labels = letters[12 + (1:3)])
p_DA_2_2 = plot_tables(v_files_2_DA_tab, v_labels = letters[12 + (4:10)], 
	nrow = 5)
p_DA_2 = ggpubr::ggarrange(p_DA_2_1, p_DA_2_2, nrow = 2, ncol = 1, 
	heights = c(0.25, 0.75) )
pdf('Figure_S3.pdf', width = 7, height = 10)
	print(p_DA_1)
	print(p_DA_2)
dev.off()


v_files_3_EN_fig   = c(grep1('__barplot\\.'), grep1('__heatmap\\.'))
v_files_3_EN_tab   = grep1('xlsx_png\\/(CHIP|motifs|func|peaks|genes|chrom)')
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
pdf('Figure_S4.pdf', width = 7, height = 10)
	print(p_EN_1)
	print(p_EN_2)
dev.off()


# Fig S2: HA-HE gene tracks

# remotes::install_github(repo = "poisonalien/trackplot")
cd ~/install/R/github
git clone https://github.com/CRG-Barcelona/libbeato.git
git clone https://github.com/CRG-Barcelona/bwtool.git
cd libbeato/
./configure
make
sudo make install
cd ../bwtool/
./configure
make
sudo make install

library(trackplot)

install.packages('Signac')


library(Signac)
library(GenomicRanges)

v_bw = list.files('panels', pattern = '*.bw', full.names = T) %>% 
		setNames(., gsub('.*\\/', '', .) %>% gsub('\\.bw', '', .))

inx3_loci = "chr10:5995156-5999243"
OrganismDb 

BigwigTrack(
  region,
  bigwig,



v_bw = R.utils::getAbsolutePath(list.files('../../docs/examples/bw', full.names = T)) %>% 
		setNames(., gsub('.*\\/', '', .) %>% gsub('\\.bw', '', .))


library(EnsDb.Hsapiens.v86)


library(BSgenome.Celegans.UCSC.ce11)

l_bw = list.files('../../docs/examples/bw', full.names = T) %>% 
		setNames(., gsub('.*\\/', '', .) %>% gsub('\\.bw', '', .)) %>% as.list

gr <- GRanges(seqnames = 'chr10',
              ranges = IRanges(start = 5995156, end = 5999243)) # inx3
genome(gr) <- 'ce11'


p1 = BigwigTrack(gr, l_bw)
p2 = AnnotationPlot(atac_small, region = gr)
p3 = CombineTracks(plotlist = list(p1, p2), heights = c(1, 1))

pdf('tmp.pdf', width = 7, height = 10)
	print(p3)
dev.off()

library(Seurat)
set.seed(123) # For reproducibility
# Simulate expression data: 100 genes x 10 cells
data_matrix <- matrix(rpois(1000, lambda = 1), nrow = 100, ncol = 10)

# Assign row and column names to simulate gene names and cell IDs
rownames(data_matrix) <- paste0("Gene", seq_len(nrow(data_matrix)))
colnames(data_matrix) <- paste0("Cell", seq_len(ncol(data_matrix)))
seurat_object <- CreateSeuratObject(counts = data_matrix)
seurat_object$genome <- "ce11"

p2 = AnnotationPlot(seurat_object, region = gr)

gr <- GRanges(seqnames = 'X',
              ranges = IRanges(start = 5995156, end = 5999243))

p1 = BigwigTrack(gr, l_bw)


#genome : "hg19" 
gen<-genome(cpgIslands)
#Chromosme name : "chr7"
chr <- as.character(unique(seqnames(cpgIslands)))
#Ideogram track
itrack <- IdeogramTrack(genome = gen, chromosome = chr)

gtrack <- GenomeAxisTrack()


bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr12",
  start = afrom, end = ato, filter = list(with_ox_refseq_mrna = TRUE), 
  stacking = "dense")


p1 = BigwigTrack(region = 'chr10:5995156-5999243', l_bw)



gr_inx3 <- GRanges(seqnames = 'chrX',
              ranges = IRanges(start = 5995156, end = 5999243))

p1 = BigwigTrack(region = 'chr10:5995156-5999243', l_bw)



gr_inx3 <- GRanges(seqnames = 'chr10',
              ranges = IRanges(start = 5995156, end = 5999243))

p1 = BigwigTrack(region = 'chr10:5995156-5999243', l_bw)




v_bw = list.files('panels', pattern = '*.bw', full.names = T) %>% 
		setNames(., gsub('.*\\/', '', .) %>% gsub('\\.bw', '', .))
gr <- GRanges(seqnames = 'chr10',
              ranges = IRanges(start = 5995156, end = 5999243)) # inx3
genome(gr) <- 'ce11'




dTrack2 <- DataTrack(range = v_bw[1], genome = "hg19",
                     type = "l", chromosome = "chr19", name = "bedGraph")
p1 = plotTracks(dTrack2)

pdf('tmp.pdf', width = 7, height = 10)
	print(p3)
dev.off()


# BiocManager::install("Gviz")

library(Gviz)
library(TxDb.Celegans.UCSC.ce11.ensGene)
txdb = TxDb.Celegans.UCSC.ce11.ensGene
exons = exons(txdb)

test = transcripts(txdb)
test = genes(txdb)

grtrack <- GeneRegionTrack(exons, genome = 'ce11',
                           chromosome = 'chr10', name = "foo", 
                           transcriptAnnotation = "symbol")

grtrack <- GeneRegionTrack(exons, genome = 'ce11',
                           chromosome = 'chr10', name = "foo", 
                           transcriptAnnotation = "symbol", 
                           start = 5995156, end = 5999243)

df_genes = readRDS('../../references/v4/worm/genome/annotation/R/df_genes_metadata.rds')
df_genes = readRDS('../../references/v4/worm/genome/annotation/R/df_genes_transcripts.rds')


library(Gviz)
library(TxDb.Celegans.UCSC.ce11.ensGene)
txdb = TxDb.Celegans.UCSC.ce11.ensGene
exons = exons(txdb)

txTr <- GeneRegionTrack(exons, chromosome = 'chr10', start = 5995156, end = 5999243)
p1 = plotTracks(txTr)

pdf('tmp.pdf', width = 7, height = 10)
	print(grtrack)
dev.off()

genome = 'ce11'
chromosome = 'chrX'
start = 5995156
end   = 5999243

ensGenes <- UcscTrack(genome = genome, chromosome = chromosome, 
                      track="Ensembl Genes", table = "ensGene", 
                      from = start, to = end,
                      trackType = "GeneRegionTrack", 
                      rstarts = "exonStarts", rends = "exonEnds",
                      gene = "name", symbol = "name2", 
                      transcript = "name", strand = "strand", 
                      fill = "#960000", name = "Ensembl Genes")

v_bw = list.files('panels', pattern = '*.bw', full.names = T) %>% 
		setNames(., gsub('.*\\/', '', .) %>% gsub('\\.bw', '', .))
gr <- GRanges(seqnames = chromosome,
              ranges = IRanges(start = start, end = end)) # inx3

dTrack2 <- DataTrack(range = v_bw[1], genome = genome,
                     type = 'l', chromosome = chromosome, name = "bigwig",
                     from = start, to = end)


pdf('tmp.pdf', width = 7, height = 10)
	plotTracks(list(dTrack2, ensGenes), from = start, to = end)
dev.off()




v_bw = list.files('panels', pattern = '*.bigWig', full.names = T) %>% 
		setNames(., gsub('.*\\/', '', .) %>% gsub('\\.bigWig', '', .))
gr <- GRanges(seqnames = chromosome,
              ranges = IRanges(start = start, end = end)) # inx3

dTrack2 <- DataTrack(range = v_bw[1], genome = genome,
                     type = 'l', chromosome = chromosome, name = "bigwig",
                     from = start, to = end)


dTrack2 <- DataTrack(data = v_bw[1], genome = genome,
                     type = 'l', chromosome = chromosome, name = "bigWig",
                     from = start, to = end)


dTrack2 <- DataTrack(data = "panels/ctl_1.bigWig")
dTrack2@data


dTrack2 <- DataTrack(data = "panels/ctl_1.bw", range = gr)
dTrack2@data


dTrack2 <- DataTrack(data = "panels/ctl_1.bigWig", range = gr)
dTrack2@data



dTrack2 <- DataTrack(range = "panels/ctl_1.bw", genome = genome,
                      chromosome = chromosome,
                     from = start, to = end)
dTrack2@data
Track2@range



dTrack2 <- DataTrack(range = "panels/ctl_1.bigWig")
dTrack2@data


# set track for ChIP-seq coverage
chipseqFile=system.file("extdata",
                      "wgEncodeHaibTfbsA549.chr21.bw",
                      package="compGenomRData")

cov.track=DataTrack(chipseqFile,type = "l",
                    name="coverage")

library(rtracklayer)

bw_med1 <- import.bw(v_bw[1], as="GRanges")
bw_med1@seqnames@values %<>% paste0('chr', .) %>% as.factor
bw_med1@seqinfo@seqnames %<>% paste0('chr', .)
bw_med1@seqnames

options(ucscChromosomeNames=F)
dTrack2 <- DataTrack(bw_med1, genome = genome,
                      chromosome = chromosome,
                     from = start, to = end)

pdf('tmp.pdf', width = 7, height = 10)
	plotTracks(list(dTrack2, ensGenes), from = start, to = end)
dev.off()



# BiocManager::install("Gviz")
# BiocManager::install("org.Ce.eg.db")
library(ggbio)
library(org.Ce.eg.db)
orgdb = org.Ce.eg.db
p_txdb <- autoplot(orgdb, which = gr)


library(org.Ce.eg.db)

  region,
  bigwig,

BigwigTrack(
  region,
  bigwig,



v_bw = list.files('../../docs/examples/bw', full.names = T)

bigWigs = read_coldata(bws = v_bw, build = 'ce11')


#
v_files_2_DA_fig %<>% .[. != '../../docs/examples/png/hmg4_vs_ctl__ATAC_other_plots.png']

#Path to bigWig files
bigWigs = c("H1_Oct4.bw", "H1_Nanog.bw", "H1_k4me3.bw", 
            "H1_k4me1.bw", "H1_k27ac.bw", "H1_H2az.bw", "H1_Ctcf.bw")

#Make a table of bigWigs along with ref genome build
bigWigs = read_coldata(bws = bigWigs, build = "hg19")

inx3_loci = "chr10:5995156-5999243"
t = track_extract(colData = bigWigs, loci = inx3_loci)
p1 = track_plot(summary_list = t)





library(ggbio)
library(org.Ce.eg.db)
v_bw = list.files('panels', pattern = '*.bw', full.names = T) %>% 
		setNames(., gsub('.*\\/', '', .) %>% gsub('\\.bw', '', .))

gr <- GRanges(seqnames = 'chrX',
              ranges = IRanges(start = 5995156, end = 5999243)) # inx3
genome(gr) <- 'ce11'

library(TxDb.Celegans.UCSC.ce11.ensGene)
txdb = TxDb.Celegans.UCSC.ce11.ensGene

p1 = autoplot(txdb, which = gr)
p1 %<>% {. + theme_classic() }

pdf('tmp.pdf', width = 7, height = 10)
	print(p1)
dev.off()


install.packages("ggcoverage")

library("rtracklayer")
library("graphics")
library("ggcoverage")
library("ggpattern")

gffread -E  \
	../../references/v4/worm/genome/annotation/filtered/annotation.gff3 \
	-T -o  panels/anno_worm.gtf

dt_meta = data.table(Type = c(1:3, 1:3), 
	Group = c(rep('hmg4', 3), rep('ctl', 3)))
dt_meta[, SampleName := paste(Group, Type, sep = '_')]

track.df = LoadTrackFile(track.folder = 'panels', format = "bw",
                         region = "X:5,995,156-5,999,243", extend = 2000,
                         meta.info = dt_meta)

gtf_file = '../../references/v4/worm/genome/annotation/filtered/annotation.gff3'
gtf.gr = rtracklayer::import.gff(con = gtf_file, format = 'gff3')

gtf_file = 'panels/anno_worm.gtf'
gtf.gr = rtracklayer::import.gff(con = gtf_file, format = 'gtf')



basic.coverage = ggcoverage(data = track.df, color = "auto", plot.type = "joint",)


p1 = ggcoverage(data = track.df, color = "auto", 
	plot.type = "facet", range.position = "in", facet.y.scale = "fixed") + 
  geom_gene(gtf.gr=gtf.gr)



ggcoverage(data = track.df, color = "auto", plot.type = "facet", 
	range.position = "in", facet.y.scale = "fixed")


ggcoverage(data = track.df, color = "auto", plot.type = "facet",
                            mark.region = mark.region, range.position = "in", 
                            facet.y.scale = "fixed")


basic.coverage + 
  geom_gene(gtf.gr=gtf.gr)

track.df = LoadTrackFile(track.folder = 'panels', format = "bw",
                         region = "chrX:5,995,156-5,999,243", extend = 2000,
                         meta.info = dt_meta)

 p1 = ggcoverage(data = track.df, color = "auto",
+ plot.type = "facet", range.position = "in", facet.y.scale = "fixed") +
+   geom_gene(gtf.gr=gtf.gr)
Error in ggplot_add.gene(object, p, objectname) :
  gene_name, gene_type, gene_biotype is not in provided GTF file!
In addition: Warning message:
In geom_coverage(data = data, mapping = mapping, color = color,  :
  Fewer colors provided than there are groups in Group variable, falling back to default colors


 p1 = ggcoverage(data = track.df, color = "auto",
 	plot.type = "facet", range.position = "in", facet.y.scale = "fixed") +
 	geom_gene(gtf.gr=gtf.gr)

BiocManager::install("AnnotationHub")
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c('elegans', 'gtf', '109'))
  #            title
  # AH110645 | Caenorhabditis_elegans.WBcel235.109.abinitio.gtf
  # AH110646 | Caenorhabditis_elegans.WBcel235.109.gtf

gtf <- ah[["AH110646"]]

p1 = ggcoverage(data = track.df, color = "auto",
 	plot.type = "facet", range.position = "in", facet.y.scale = "fixed") +
 	geom_gene(gtf.gr=gtf)

# p1 %<>% {. + theme_classic() }

pdf('tmp.pdf', width = 7, height = 10)
	print(p1)
dev.off()


p1 = ggcoverage(data = track.df, color = "auto",
 	plot.type = "facet", range.position = "in", facet.y.scale = "fixed") +
 	geom_gene(gtf.gr=gtf)

# p1 %<>% {. + theme_classic() }

# p1 = geom_gene(gtf.gr=gtf)
# geom_transcript(gtf.gr=gtf)

p1 = ggcoverage(data = track.df, color = "auto",
 	plot.type = "facet", range.position = "in", facet.y.scale = "fixed") 
p1 = p1 +  geom_ideogram(genome = "ce11",plot.space = 0)

pdf('tmp.pdf', width = 7, height = 10)
	print(p1)
dev.off()


grtrack <- GeneRegionTrack(geneModels, genome = gen,
                           chromosome = chr, name = "Gene Model", 
                           transcriptAnnotation = "symbol",
                           background.panel = "#FFFEDB",
                           background.title = "darkblue")










genome = 'ce11'
chromosome = 'chrX'
start = 5995156
end   = 5999243

ensGenes <- UcscTrack(genome = genome, chromosome = chromosome, 
                      track="Ensembl Genes", table = "ensGene", 
                      from = start, to = end,
                      trackType = "GeneRegionTrack", 
                      rstarts = "exonStarts", rends = "exonEnds",
                      gene = "name", symbol = "name2", 
                      transcript = "name", strand = "strand", 
                      fill = "#960000", name = "Ensembl Genes")

library(rtracklayer)

bw_med1 <- import.bw(v_bw[1], as="GRanges")
bw_med1@seqnames@values %<>% paste0('chr', .) %>% as.factor
bw_med1@seqinfo@seqnames %<>% paste0('chr', .)
bw_med1@seqnames

options(ucscChromosomeNames=F)
dTrack2 <- DataTrack(bw_med1, genome = genome,
                      chromosome = chromosome,
                     from = start, to = end, name = names(v_bw)[1],
                     type = 'l')

pdf('tmp.pdf', width = 7, height = 10)
	plotTracks(list(dTrack2, ensGenes), from = start, to = end)
dev.off()
















# inx-3
genome = 'ce11'
chromosome = 'chrX'
start = 5995156
end   = 5999243



# cav-1
genome = 'ce11'
chromosome = 'chrIV'
start = 9772510
end   = 9773559


# zip-2
genome = 'ce11'
chromosome = 'chrIII'
start = 832982
end   = 835814

ow III:8238001..8239423 (

# tsp-1
genome = 'ce11'
chromosome = 'chrIII'
start = 8238001
end   = 8239423


# https://jnicolaus.com/tutorials-old/gviztutorial.html

library(rtracklayer)
library(purrr)

get_ensGenes <- function(genome, chromosome, start, end){
	ensGenes = UcscTrack(genome = genome, chromosome = chromosome, 
                      track="Ensembl Genes", table = "ensGene", 
                      from = start, to = end,
                      trackType = "GeneRegionTrack", 
                      rstarts = "exonStarts", rends = "exonEnds",
                      gene = "name", symbol = "name2", 
                      transcript = "name", strand = "strand", 
                      fill = "#960000", name = "Ensembl Genes",
                      transcriptAnnotation="gene")
	return(ensGenes)
}

get_bigwig_plot <- function(file, name, genome, chromosome, start, end, 
	ylim = NULL){
	bw_med1 <- import.bw(file, as="GRanges")
	bw_med1@seqnames@values %<>% paste0('chr', .) %>% as.factor
	bw_med1@seqinfo@seqnames %<>% paste0('chr', .)
	bw_med1@seqnames

	options(ucscChromosomeNames = F)
	dTrack2 <- DataTrack(bw_med1, genome = genome,
	                      chromosome = chromosome,
	                     from = start, to = end, name = name,
	                     type = 'l', ylim = ylim, 
		col.histogram  = ifelse(grepl('ctl', name), '#FDE725FF', '#440154FF'),
		fill.histogram = ifelse(grepl('ctl', name), '#FDE725FF', '#440154FF'))
	return(dTrack2)
}


# l_bw = list.files('panels', full.names = T, pattern = '*.bw') %>% 
# 		setNames(., gsub('.*\\/', '', .) %>% gsub('\\.bw', '', .)) %>% as.list
# l_bw_p = imap(l_bw, ~get_bigwig_plot(.x, .y, genome, chromosome, start, end,
# 	ylim = NULL, 
# 		col.histogram  = ifelse(grepl('ctl', .y), '#FDE725FF', '#440154FF'),
# 		fill.histogram = ifelse(grepl('ctl', .y), '#FDE725FF', '#440154FF'))))

# ensGenes = get_ensGenes(genome, chromosome, start, end)

# pdf('tmp.pdf', width = 7, height = 10)
# 	plotTracks(c(l_bw_p, list(ensGenes)), from = start, to = end)
# dev.off()

plot_samples_and_ref <- function(chromosome, start, end, title, ylim = NULL){
	
	genome = 'ce11'
	
	# gtrack <- GenomeAxisTrack()

	l_bw_p = imap(l_bw, ~get_bigwig_plot(.x, .y, genome, chromosome, start, end,
	ylim))

	ensGenes = get_ensGenes(genome, chromosome, start, end)

	plotTracks(c(l_bw_p, list(ensGenes)), from = start, to = end, main = title)

}

pdf('tmp.pdf', width = 7, height = 10)
	plot_samples_and_ref('chrII', 4170264, 4174498, title = 'R03H10.6', ylim = c(0, 500))
dev.off()


pdf('tmp.pdf', width = 7, height = 10)
	plot_samples_and_ref('chrII', 13833095, 13833631, title = 'F08G2.5')
dev.off()


w II:13833095..13833631 (

# tsp-1
genome = 'ce11'
chromosome = 'chrIII'
start = 8237001
end   = 8239423


# F08G2.5
genome = 'ce11'
chromosome = 'chrIII'
start = 8237001
end   = 8239423
dow II:.. (Le


# R03H10.6
genome = 'ce11'
chromosome = 'chrII'
start = 4170264
end   = 4174498

l_bw_p = imap(l_bw, ~get_bigwig_plot(.x, .y, genome, chromosome, start, end,
	ylim = NULL))

ensGenes = get_ensGenes(genome, chromosome, start, end)

pdf('tmp.pdf', width = 7, height = 10)
	plotTracks(c(l_bw_p, list(ensGenes)), from = start, to = end)
dev.off()


