#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
id = args[1]

cur_files = list.files('.')

v_reads <- v_peaks <- c()

c1 = 0
for(percent in seq(10,100,10)){
  c1 = c1+1
  pattern_base = paste0(percent,'_percent')

  reads_file_name = grep(paste0(pattern_base, '*.bed$'), cur_files, value = T)
  reads_file <- read.table(reads_file_name, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

  peaks_file_name = grep(paste0(pattern_base, '.*.narrowPeak$'), cur_files, value = T)
  peaks_file <- tryCatch(read.table(peaks_file_name, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""), error=function(e) NULL)

  v_reads[c1] = nrow(reads_file) / 2

  if(length(peaks_file) != 0) { v_peaks[c1] = nrow(peaks_file) } else { v_peaks[c1] = 0}

}

pdf(file = paste0(id, '__saturation_curve.pdf'))
  plot(v_reads, v_peaks, type = 'b', xlab = 'number of reads', ylab = 'number of peaks', main = id)
dev.off()
