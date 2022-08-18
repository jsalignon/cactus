#!/bin/bash

BW_FILES=`ls *.bw`
LABELS=`echo ${BW_FILES} | sed 's/_raw.bw//g'`

gDNA_PRESENT=${1}
BLAKLISTED=${2}
BINSIZE=${3}

multiBigwigSummary bins \
  --bwfiles $BW_FILES \
  --labels $LABELS \
  --outFileName correlation.npz \
  --blackListFileName ${BLAKLISTED} \
  --binSize ${BINSIZE}

plotPCA \
  --corData correlation.npz \
  --transpose \
  --ntop 1000 \
  --plotHeight 10 \
  --plotWidth 7 \
  --plotTitle pca_raw_top1000_${gDNA_PRESENT} \
  --plotFile pca_raw_top1000_${gDNA_PRESENT}_pca.pdf

plotPCA \
  --corData correlation.npz \
  --transpose \
  --ntop 1000 \
  --log2 \
  --plotHeight 10 \
  --plotWidth 7 \
  --plotTitle pca_log2_top1000_${gDNA_PRESENT} \
  --plotFile pca_log2_top1000_${gDNA_PRESENT}_pca.pdf

plotPCA \
  --corData correlation.npz \
  --transpose \
  --ntop 5000 \
  --log2 \
  --plotHeight 10 \
  --plotWidth 7 \
  --plotTitle pca_log2_top5000_${gDNA_PRESENT} \
  --plotFile pca_log2_top5000_${gDNA_PRESENT}_pca.pdf

plotPCA \
  --corData correlation.npz \
  --transpose \
  --ntop 100 \
  --plotHeight 10 \
  --plotWidth 7 \
  --plotTitle pca_raw_top100_${gDNA_PRESENT} \
  --plotFile pca_raw_top100_${gDNA_PRESENT}_pca.pdf

plotPCA \
  --corData correlation.npz \
  --transpose \
  --log2 \
  --ntop 100 \
  --plotHeight 10 \
  --plotWidth 7 \
  --plotTitle pca_log2_top100_${gDNA_PRESENT} \
  --plotFile pca_log2_top100_${gDNA_PRESENT}_pca.pdf

# heatmap spearman without outliers
plotCorrelation \
  --corData correlation.npz \
  --whatToPlot heatmap \
  --plotNumbers \
  --skipZeros \
  --corMethod spearman \
  --removeOutliers \
  --plotTitle spearman_correlation_heatmap_without_outliers_${gDNA_PRESENT} \
  --plotFile spearman_correlation_heatmap_without_outliers_${gDNA_PRESENT}_cor.pdf

# heatmap spearman with outliers
plotCorrelation \
  --corData correlation.npz \
  --whatToPlot heatmap \
  --plotNumbers \
  --skipZeros \
  --corMethod spearman \
  --plotTitle spearman_correlation_heatmap_with_outliers_${gDNA_PRESENT} \
  --plotFile spearman_correlation_heatmap_with_outliers_${gDNA_PRESENT}_cor.pdf

# heatmap pearson without outliers
plotCorrelation \
  --corData correlation.npz \
  --whatToPlot heatmap \
  --plotNumbers \
  --skipZeros \
  --corMethod pearson \
  --removeOutliers \
  --plotTitle pearson_correlation_heatmap_without_outliers_${gDNA_PRESENT} \
  --plotFile pearson_correlation_heatmap_without_outliers_${gDNA_PRESENT}_cor.pdf

# heatmap pearson with outliers
plotCorrelation \
  --corData correlation.npz \
  --whatToPlot heatmap \
  --plotNumbers \
  --skipZeros \
  --corMethod pearson \
  --plotTitle pearson_correlation_heatmap_with_outliers_${gDNA_PRESENT} \
  --plotFile pearson_correlation_heatmap_with_outliers_${gDNA_PRESENT}_cor.pdf

# # scatterplot spearman with outliers
# plotCorrelation \
#   --corData correlation.npz \
#   --whatToPlot scatterplot \
#   --skipZeros \
#   --corMethod spearman \
#   --plotTitle spearman_correlation_scatterplot_with_outliers \
#   --plotFile spearman_correlation_scatterplot_with_outliers_${gDNA_PRESENT}_cor.pdf

  mv correlation.npz correlation_${gDNA_PRESENT}.npz






