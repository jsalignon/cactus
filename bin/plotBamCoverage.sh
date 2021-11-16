#!/bin/bash

bam=$1
id=$2
binsize_bigwig_creation=$3
nb_threads=$4
blacklisted_regions=$5
effective_genome_size=$6
nb_1bp_site_to_sample_for_coverage=$7

bamCoverage \
  --bam ${bam} \
  --outFileName ${id}_raw.bw \
  --binSize ${binsize_bigwig_creation} \
  --numberOfProcessors ${nb_threads} \
  --blackListFileName ${blacklisted_regions} \
  --effectiveGenomeSize ${effective_genome_size}

bamCoverage \
  --bam ${bam} \
  --outFileName ${id}_RPGC_norm.bw \
  --binSize ${binsize_bigwig_creation} \
  --numberOfProcessors ${nb_threads} \
  --blackListFileName ${blacklisted_regions} \
  --normalizeUsing RPGC \
  --effectiveGenomeSize ${effective_genome_size}

plotCoverage  \
  --bam ${bam} \
  --blackListFileName ${blacklisted_regions} \
  --numberOfProcessors ${nb_threads} \
  --numberOfSamples ${nb_1bp_site_to_sample_for_coverage} \
  --plotTitle ${id}_coverage \
  --plotFile ${id}_coverage.pdf
  