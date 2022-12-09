#!/usr/bin/env nextflow


//// List of all processes in order

// ATAC_reads__merging_reads
// ATAC_reads__trimming_reads
// ATAC_reads__aligning_reads
// ATAC_reads__removing_low_quality_reads
// ATAC_reads__marking_duplicated_reads
// ATAC_reads__removing_duplicated_reads
// ATAC_reads__removing_reads_in_mitochondria_and_small_contigs
// ATAC_reads__converting_bam_to_bed_and_adjusting_for_Tn5

// ATAC_QC_reads__running_fastqc
// ATAC_QC_reads__computing_bigwig_tracks_and_plotting_coverage
// ATAC_QC_reads__computing_and_plotting_bigwig_tracks_correlations
// ATAC_QC_reads__plotting_insert_size_distribution
// ATAC_QC_reads__sampling_aligned_reads
// ATAC_QC_reads__measuring_overlap_with_genomic_regions
// ATAC_QC_reads__estimating_library_complexity
// ATAC_QC_reads__sampling_trimmed_reads
// ATAC_QC_reads__aligning_sampled_reads
// ATAC_QC_reads__gathering_all_stat
// ATAC_QC_reads__gathering_all_samples
// ATAC_QC_reads__splitting_stat_for_multiqc
// ATAC_QC_reads__running_multiQC

// ATAC_peaks__calling_peaks
// ATAC_peaks__splitting_multi_summits_peaks
// ATAC_peaks__removing_blacklisted_regions
// ATAC_peaks__removing_input_control_peaks
// ATAC_peaks__removing_specific_regions

// ATAC_QC_peaks__computing_and_plotting_saturation_curve
// ATAC_QC_peaks__annotating_macs2_peaks
// ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_each_sample
// ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_all_samples_grouped

// MRNA__quantifying_transcripts_abundances

// MRNA_QC__running_fastqc
// MRNA_QC__running_MultiQC

// DA_ATAC__doing_differential_abundance_analysis
// DA_ATAC__annotating_diffbind_peaks
// DA_ATAC__plotting_differential_abundance_results
// DA_ATAC__saving_detailed_results_tables

// DA_mRNA__doing_differential_abundance_analysis
// DA_mRNA__plotting_differential_abundance_results
// DA_mRNA__saving_detailed_results_tables

// DA_split__splitting_differential_abundance_results_in_subsets
// DA_split__plotting_venn_diagrams

// Enrichment__computing_functional_annotations_overlaps
// Enrichment__computing_genes_self_overlaps
// Enrichment__computing_peaks_overlaps
// Enrichment__computing_motifs_overlaps
// Enrichment__reformatting_motifs_results
// Enrichment__computing_enrichment_pvalues

// Figures__making_enrichment_barplots
// Figures__making_enrichment_heatmap
// Figures__merging_pdfs

// Tables__formatting_csv_tables
// Tables__merging_csv_tables
// Tables__saving_excel_tables





//// shortening certain oftenly used parameters

res_dir = params.res_dir
pub_mode = params.pub_mode
out_processed = "${res_dir}/Processed_Data"
out_fig_indiv = "${res_dir}/Figures_Individual"
out_fig_merge = "${res_dir}/Figures_Merged"
out_tab_merge = "${res_dir}/Tables_Merged"
out_tab_indiv = "${res_dir}/Tables_Individual"


//// setting up parameters

save_all_fastq  = params.save_fastq_type == 'all' ? true : false
save_all_bam    = params.save_bam_type   == 'all' ? true : false
save_all_bed    = params.save_bed_type   == 'all' ? true : false
save_last_fastq = params.save_fastq_type in ['last', 'all'] ? true : false
save_last_bam   = params.save_bam_type   in ['last', 'all'] ? true : false
save_last_bed   = params.save_bed_type   in ['last', 'all'] ? true : false

do_atac = params.experiment_types in ['atac', 'both']
do_mRNA = params.experiment_types in ['mRNA', 'both']

params.do_any_enrichment = !params.disable_all_enrichments




//// Creating empty channels

Overlap_tables_channel = Channel.empty()
Formatting_csv_tables_channel = Channel.empty()
Exporting_to_Excel_channel = Channel.empty()
Merging_pdfs_channel = Channel.empty()


//// Creating mRNA-Seq channels

static def returnR2ifExists(r2files) {
  boolean exists = r2files[1].exists();
  boolean diff   = r2files[0] != r2files[1];
  if (exists & diff) { return(r2files)
  } else { return(r2files[0]) }
}

// MRNA_reads_for_merging = Channel.create() // not yet implemented
MRNA_reads_for_kallisto = Channel.create()

fastq_files = file(params.design__mrna_fastq)
Channel
  .from(fastq_files.readLines())
  .map{ it.split() }
  .map{ [ it[0], it[1..-1] ] }
  .transpose()
  .map{ [ it[0], file(it[1]), file(it[1].replace("_R1", "_R2") ) ] }
  .map{ [ it[0], returnR2ifExists(it[1, 2]) ] }
  .dump(tag:'mrna_fastq')
  .tap{ MRNA_reads_for_running_fastqc }
  .map{ it.flatten() }
  .map{ [it[0], it[1..-1] ] }
  .dump(tag:'mrna_kallisto') {"mRNA reads for kallisto: ${it}"}
  .set { MRNA_reads_for_kallisto }


//// Creating ATAC-Seq channels

ATAC_reads_for_merging = Channel.create()
ATAC_reads_for_trimming_1 = Channel.create()

fastq_files = file(params.design__atac_fastq)
Channel
  .from(fastq_files.readLines())
  .map{ it.split() }
  .map{ [ it[0], it[1..-1] ] }
  .transpose()
  .map{ [ it[0], file(it[1]), file(it[1].replace("_R1", "_R2")) ] }
  .tap{ Raw_ATAC_reads_for_running_fastqc }
  .map{ [ it[0], [ it[1], it[2] ] ] }
  .transpose()
  .groupTuple()
  .dump(tag:'atac_raw') {"ATAC raw reads: ${it}"}
  .branch{
       ATAC_reads_for_merging: it[1].size() > 2
    ATAC_reads_for_trimming_1: it[1].size() == 2
  }
  .set{ atac_raw }
  // => if it[1] has 2 elements it means there is just the R1 and R2 files and
  // thus no other file has the same sample_id, so we don't need to merge, if 
  // more it means there are different runs so we merge them, if less it is not 
  // an option since we need paired end reads



//// Let's start!


process ATAC_reads__merging_reads {
  tag "${id}"

  label "fastqc"

  publishDir path: "${out_processed}/1_Preprocessing/ATAC__reads__fastq_merged",
             mode: "${pub_mode}", enabled: save_all_fastq

  when: 
    do_atac
  
  input:
    set id, file(files_R1_R2) from atac_raw.ATAC_reads_for_merging

  output:
    set id, file('*__R1_merged.fastq.gz'), file('*__R2_merged.fastq.gz') \
      into ATAC_reads_for_trimming_2

  script:
    """
    
    cat `ls *R1* | sort` > ${id}__R1_merged.fastq.gz
    cat `ls *R2* | sort` > ${id}__R2_merged.fastq.gz
    
    """
}


atac_raw.ATAC_reads_for_trimming_1
  .map{ it.flatten() }
  .mix( ATAC_reads_for_trimming_2 )
  .dump(tag:'atac_trimming') {"ATAC reads for trimming: ${it}"}
  .set { ATAC_reads_for_trimming_3 }



process ATAC_reads__trimming_reads {
  tag "${id}"

  label "skewer_pigz"
  cpus params.pigz__nb_threads

  publishDir \
    path: "${out_processed}/1_Preprocessing/ATAC__reads__fastq_trimmed",
    mode: "${pub_mode}", pattern: "*.log"

  publishDir \
    path: "${out_processed}/1_Preprocessing/ATAC__reads__fastq_trimmed",
    mode: "${pub_mode}", pattern: "*.fastq.gz", enabled: save_last_fastq

  when: 
    do_atac

  input:
    set id, file(read1), file(read2) from ATAC_reads_for_trimming_3

  output:
    set id, file("*__R1_trimmed.fastq.gz"), file("*__R2_trimmed.fastq.gz") \
             into Trimmed_reads_for_aligning, Trimmed_reads_for_sampling
    set val('trimmed'), id, file("*__R1_trimmed.fastq.gz"),
        file("*__R2_trimmed.fastq.gz") into Trimmed_ATAC_reads_for_running_fastqc
    file("*.log")

  shell:
    '''

    R1=!{read1}
    R2=!{read2}
    id=!{id}
    pigz__nb_threads=!{params.pigz__nb_threads}

    R1TRIM=$(basename ${R1} .fastq.gz)__R1_trimmed.fastq
    R2TRIM=$(basename ${R2} .fastq.gz)__R2_trimmed.fastq

    skewer --quiet -x CTGTCTCTTATA -y CTGTCTCTTATA -m pe ${R1} ${R2} \
           -o ${id} > trimer_verbose.txt

    mv ${id}-trimmed.log ${id}__skewer_trimming.log
    mv ${id}-trimmed-pair1.fastq $R1TRIM
    mv ${id}-trimmed-pair2.fastq $R2TRIM

    pigz -p ${pigz__nb_threads} $R1TRIM
    pigz -p ${pigz__nb_threads} $R2TRIM

    pigz -l $R1TRIM.gz > ${id}__pigz_compression.log
    pigz -l $R2TRIM.gz >> ${id}__pigz_compression.log

    '''
}

// trimming adaptors

// this line crashed the script somehow. I don't really get the grep fast here 
// anyway so I changed it
// pigz -l $R2TRIM.gz | grep fastq - >> !{id}_pigz_compression.log

// cp $R1 $R1TRIM.gz
// cp $R2 $R2TRIM.gz
// touch tmp.log

// R1TRIM=$(basename ${R1} .fastq.gz)_trim.fastq
// R2TRIM=$(basename ${R2} .fastq.gz)_trim.fastq
// 

// skewer potentially useful options
// -m, --mode <str> trimming mode; 
//             1) single-end -- head: 5' end; tail: 3' end; any: anywhere
//             2) paired-end -- pe: paired-end; mp: mate-pair; ap: amplicon (pe)
// -q, --end-quality  <int> Trim 3' end until specified or higher quality 
//      reached; (0). => Maybe one could later set -q 20. But on 
//      "https://informatics.fas.harvard.edu/atac-seq-guidelines.html" they say: 
//      "Other than adapter removal, we do not recommend any trimming of the 
//      reads. Such adjustments can complicate later steps, such as the 
//      identification of PCR duplicates."
// -t, --threads <int>   Number of concurrent threads [1, 32]; (1)
// -A, --masked-output  Write output file(s) for trimmed reads (trimmed bases 
//                      converted to lower case) (no)

// pigz potentially useful options
// -k, --keep           Do not delete original file after processing
// -0 to -9, -11        Compression level (level 11, zopfli, is much slower)
//                       => the default is 6
// --fast, --best       Compression levels 1 and 9 respectively
// -p, --processes n    Allow up to n compression threads 
//              (default is the number of online processors, or 8 if unknown)
// -c, --stdout         Write all processed output to stdout (won't delete)
// -l, --list           List the contents of the compressed input



Raw_ATAC_reads_for_running_fastqc
  .map{ [ 'raw', it[0], it[1], it[2] ] }
  .mix(Trimmed_ATAC_reads_for_running_fastqc)
  .dump(tag:'atac_fastqc') {"ATAC reads for fastqc: ${it}"}
  .set{ All_ATAC_reads_for_running_fastqc}



process ATAC_reads__aligning_reads {
  tag "${id}"

  label "bowtie2_samtools"
  cpus params.botwie2__nb_threads

  publishDir path: "${out_processed}/1_Preprocessing/ATAC__reads__bam",
             mode: "${pub_mode}", pattern: "*.{txt,qc}"
             
  publishDir path: "${out_processed}/1_Preprocessing/ATAC__reads__bam",
             mode: "${pub_mode}", pattern: "*.bam", enabled: save_all_bam

  when: 
    do_atac

  input:
    set id, file(read1), file(read2) from Trimmed_reads_for_aligning

  output:
    set id, file("*.bam") into Bam_for_removing_low_quality_reads,
                               Bam_for_sampling
    file("*.txt")
    file("*.qc")

  script:
    
    def new_bam = id + ".bam"
    
    """

    bowtie2 -p ${params.botwie2__nb_threads} \
      --very-sensitive \
      --end-to-end \
      --no-mixed \
      -X 2000 \
      --met 1 \
      --met-file "${id}_bowtie2_align_metrics.txt" \
      -x "${params.bowtie2_indexes}" \
      -q -1 "${read1}" -2 "${read2}" \
      | samtools view -bS -o ${new_bam} -

    samtools flagstat ${new_bam} > "${id}.qc"

    """
}



process ATAC_reads__removing_low_quality_reads {
  tag "${id}"

  label "bowtie2_samtools"

  publishDir path: "${out_processed}/1_Preprocessing/ATAC__reads__bam_no_lowQ", 
             mode: "${pub_mode}", pattern: "*.qc"

  publishDir path: "${out_processed}/1_Preprocessing/ATAC__reads__bam_no_lowQ", 
             mode: "${pub_mode}", pattern: "*.bam", enabled: save_all_bam

  when: 
    do_atac

  input:
    set id, file(bam) from Bam_for_removing_low_quality_reads

  output:
    set id, file("*.bam") into Bam_for_marking_duplicated_reads
    file("*.qc")

  script:

    def key = id + "__filter_LQ"
    def new_bam = key + ".bam"

    """

    samtools view -F 1804 \
                  -b \
                  -q "${params.sam_MAPQ_threshold}" \
                  ${bam} \
                  | samtools sort - -o ${new_bam}

    samtools flagstat ${new_bam} > "${key}.qc"

    """
}

// => filtering bad quality reads: unmapped, mate unmapped, 
//                                 no primary alignment, low MAPQ



process ATAC_reads__marking_duplicated_reads {
  tag "${id}"

  label "picard"

  when: 
    do_atac

  input:
    set id, file(bam) from Bam_for_marking_duplicated_reads

  output:
    set id, file("*.bam") into Bam_for_removing_duplicated_reads
    file("*.qc")

  script:

    def key = id + "__dup_marked"

    """

    picard -Xmx${params.memory_picard} MarkDuplicates \
      -INPUT "${bam}" \
      -OUTPUT "${key}.bam" \
      -METRICS_FILE "${key}.qc" \
      -VALIDATION_STRINGENCY LENIENT \
      -ASSUME_SORTED true \
      -REMOVE_DUPLICATES false \
      -TMP_DIR "."

    """
}



process ATAC_reads__removing_duplicated_reads {
  tag "${id}"

  label "bowtie2_samtools"

  publishDir \
    path: "${out_processed}/1_Preprocessing/ATAC__reads__bam_no_lowQ_dupli", 
    mode: "${pub_mode}", pattern: "*.qc"

  publishDir \
    path: "${out_processed}/1_Preprocessing/ATAC__reads__bam_no_lowQ_dupli", 
    mode: "${pub_mode}", pattern: "*.bam", enabled: save_all_bam

  when: 
    do_atac

  input:
    set id, file(bam) from Bam_for_removing_duplicated_reads

  output:
    set id, file("*.bam"), file("*.bai") \
          into Bam_for_computing_bigwig_tracks_and_plotting_coverage, 
               Bam_for_removing_mitochondrial_reads
    file("*.qc")

  script:
    def key = id + "__dup_rem"
    def new_bam = key + ".bam"

    """

    samtools view -F 1804 "${bam}" -b -o ${new_bam}
    samtools index -b ${new_bam}
    samtools flagstat ${new_bam} > "${key}.qc"

    """
}

// => removing duplicates, index bam files and generates final stat file





process ATAC_reads__removing_reads_in_mitochondria_and_small_contigs {
  tag "${id}"

  label "bowtie2_samtools"

  publishDir \
   path: "${out_processed}/1_Preprocessing/ATAC__reads__bam_no_lowQ_dupli_mito", 
   mode: "${pub_mode}", pattern: "*.{qc,txt}"

  publishDir \
   path: "${out_processed}/1_Preprocessing/ATAC__reads__bam_no_lowQ_dupli_mito", 
   mode: "${pub_mode}", pattern: "*.bam", enabled: save_last_bam

  when: 
    do_atac

  input:
    set id, file(bam), file(bai) from Bam_for_removing_mitochondrial_reads

  output:
    set id, file("*.bam") \
      into Bam_for_plotting_inserts_distribution, 
           Bam_for_converting_bam_to_bed_and_adjusting_for_Tn5
    file "*_reads_per_chrm_before_removal.txt"
    file "*_reads_per_chrm_after_removal.txt"
    file("*.qc")

  shell:
    '''

    id=!{id}
    bam=!{bam}
    chromosomes_sizes=!{params.chromosomes_sizes}

    key="${id}__no_mito"

    new_bam="${key}.bam"

    samtools view ${bam} | awk " $1 !~ /@/ {print $3}" - | uniq -c \
      > "${id}_reads_per_chrm_before_removal.txt"

    regions_to_keep=$(cut -f1 $chromosomes_sizes | paste -sd " ")

    samtools view -Sb ${bam} $regions_to_keep | tee ${new_bam} \
    | samtools view - | awk '{print $3}' OFS='\t' \
    | uniq -c > "${id}_reads_per_chrm_after_removal.txt"

    samtools index -b ${new_bam}
    samtools flagstat ${new_bam} > "${key}.qc"

    '''
}


// cut -f1 $chromosomes_sizes
// samtools view  ${bam} $regions_to_keep | awk '{print $3}' OFS='\t' - \
//     | sort | uniq -c
// samtools view  ${bam} | awk '{print $3}' OFS='\t' - | sort | uniq -c



process ATAC_reads__converting_bam_to_bed_and_adjusting_for_Tn5 {
  tag "${id}"

  label "samtools_bedtools_perl"

  publishDir \
    path: "${out_processed}/1_Preprocessing/ATAC__reads__bam_asBed_atacShift", 
    mode: "${pub_mode}", pattern: "*.bam", enabled: params.save_1bp_bam

  when: 
    do_atac

  input:
    set id, file(bam) from Bam_for_converting_bam_to_bed_and_adjusting_for_Tn5

  output:
    set id, file("*.bam*") into Reads_in_bam_files_for_diffbind
    set id, file("*.bed") into \
      Reads_in_bed_files_for_gathering_reads_stat, 
      Reads_in_bed_files_for_calling_peaks, 
      Reads_in_bed_files_for_computing_and_plotting_saturation_curve

  script:
    def key = id + "__1bp_shifted_reads"

    """

    bamToBed_and_atacShift.sh ${bam} ${key} ${params.chromosomes_sizes}

    """
}

// Examples of :
// input: hmg4_1__no_mito.bam
// outpus: hmg4_1__1bp_shifted_reads.bam  hmg4_1__1bp_shifted_reads.bam.bai  
//         hmg4_1__1bp_shifted_reads.bed  



// => converting bam to bed, adjusting for the shift of the transposase 
//    (atac seq) and keeping only a bed format compatible with macs2


// input are bam files (aligned reads) that have been filtered out of reads 
//    that are low quality, duplicates, and in mitochondria or small contigs
// output are bed files in which the reads have been converted to 1 base pair 
//    (by keeping only the 5' end) and that have been shifted to take into 
//    account the transposase shift

// Example of a read pair:

// input:
// samtools view n301b170_2_no_mito.bam | grep SRR11607686.9081238 - | cut -f1-9
//       QNAME             FLAG           RNAME    POS    MAPQ    CIGAR  RNEXT    PNEXT   TLEN
// SRR11607686.9081238     163     211000022278033 357     31      37M     =       408     88
// SRR11607686.9081238     83      211000022278033 408     31      37M     =       357     -88

// output:
// grep "SRR11607686.9081238/2\|SRR11607686.9081238/1" n301b170_2_1bp_shifted_reads.bed
// 211000022278033 360     361     SRR11607686.9081238/2   88      +
// 211000022278033 438     439     SRR11607686.9081238/1   -88     -

// reads on the + strand are shifted by + 4 base pair -1 
//    (different coordinate system) (357 + 4 -1 = 360)
// reads on the - strand are shifted by - 5 base pair -1 
//    (different coordinate system) (445 -5 - 1 = 439)

// explanations for the input:
// the read has a size of 37 bp and the insert has a size of 88 bp 
// ranges: 
// - read 1 [357-394 -> size: 37]
// - read 2 [408-445 -> size: 37]
// - insert [357-445 -> size: 88]
// - fragment [345-457 -> size: 112] (not really sure for this. I just added 
//                                     12 for the ATAC-Seq adapter length)
// https://www.frontiersin.org/files/Articles/77572/fgene-05-00005-HTML/image_m/fgene-05-00005-g001.jpg
// https://samtools-help.narkive.com/BpGYAdaT/isize-field-determining-insert-sizes-in-paired-end-data

// Note: here the 1 base pair long (5') transposase-shifted peaks are sent to 
///      DiffBind directly for Differential Accessibility Analysis.

// https://www.nature.com/articles/nmeth.2688
// https://www.nature.com/articles/s42003-020-01403-4
// The Tn5 insertion positions were determined as the start sites of reads 
//    adjusted by the rule of “forward strand +4 bp, negative strand −5 bp”4.
// 


process ATAC_QC_reads__running_fastqc {
  tag "${id}"

  label "fastqc"
  cpus params.fastqc__nb_threads

  publishDir \
    path: "${out_processed}/1_Preprocessing", mode: "${pub_mode}", 
    pattern: "*.html", saveAs: {
           if (reads_type == 'raw')     "ATAC__reads__fastqc_raw/${it}"
      else if (reads_type == 'trimmed') "ATAC__reads__fastqc_trimmed/${it}"
    }

  when: 
    do_atac

  input:
    set val(reads_type), id, file(read1), file(read2) \
        from All_ATAC_reads_for_running_fastqc

  output:
    file("*.{zip, html}") into ATAC_fastqc_reports_for_multiqc

  script:
    """

    fastqc -t ${params.fastqc__nb_threads} ${read1} ${read2}

    """
}



process ATAC_QC_reads__computing_bigwig_tracks_and_plotting_coverage {
  tag "${id}"

  label "deeptools"
  cpus params.deeptools__nb_threads

  publishDir path: "${res_dir}", mode: "${pub_mode}", saveAs: {
    if (it.indexOf(".pdf") > 0) 
        "Figures_Individual/1_Preprocessing/ATAC__reads__coverage/${it}"
    else if (it.indexOf("_raw.bw") > 0) 
        "Processed_Data/1_Preprocessing/ATAC__reads__bigwig_raw/${it}"
  }

  when: 
    do_atac & params.do_bigwig

  input:
    set id, file(bam), file(bai) \
      from Bam_for_computing_bigwig_tracks_and_plotting_coverage

  output:
    set val("ATAC__reads__coverage"), val('1_Preprocessing'), file("*.pdf") \
      into ATAC_reads_coverage_plots_for_merging_pdfs
    file("*.bw") into Bigwigs_for_correlation_1 optional true

  script:
    def key = id + "__reads_coverage"

    """

    bamCoverage \
      --bam ${bam} \
      --binSize             ${params.deeptools__binsize_bigwig_creation} \
      --normalizeUsing      ${params.deeptools__normalization_method} \
      --numberOfProcessors  ${params.deeptools__nb_threads} \
      --blackListFileName   ${params.blacklisted_regions} \
      --effectiveGenomeSize ${params.effective_genome_size} \
      --outFileName ${id}.bw

    plotCoverage \
      --bam ${bam} \
      --numberOfSamples    ${params.deeptools__nb_of_1_bp_samples} \
      --numberOfProcessors ${params.deeptools__nb_threads} \
      --blackListFileName  ${params.blacklisted_regions} \
      --plotTitle ${key} \
      --plotFile ${key}.pdf

    """
}

Merging_pdfs_channel = Merging_pdfs_channel
  .mix(ATAC_reads_coverage_plots_for_merging_pdfs.groupTuple(by: [0, 1]))




Bigwigs_for_correlation_1
  .collect()
  .into{ Bigwigs_with_input_control; Bigwigs_without_input_control_1 }

Bigwigs_without_input_control_1
  .flatten()
  .filter{ !(it =~ /input/) }
  .collect()
  .map{ [ 'without_control', it ] }
  .set{ Bigwigs_without_input_control_1 }

Bigwigs_with_input_control
  .map{ [ 'with_control', it ] }
  .concat(Bigwigs_without_input_control_1)
  .dump(tag:'bigwigs') {"bigwigs for cor and PCA: ${it}"}
  .set{ Bigwigs_for_correlation_2 }


process ATAC_QC_reads__computing_and_plotting_bigwig_tracks_correlations {

  tag "${input_control_present}"

  label "deeptools"

  publishDir \
    path: "${out_fig_indiv}/${out_path}", mode: "${pub_mode}", saveAs: { 
         if      (it.indexOf("_pca.pdf") > 0) "ATAC__reads__PCA/${it}" 
         else if (it.indexOf("_cor.pdf") > 0) "ATAC__reads__correlations/${it}" 
    }

  when: 
    do_atac

  input:
    val out_path from Channel.value('1_Preprocessing') 
    set input_control_present, file("*") from Bigwigs_for_correlation_2

  output:
    set val("ATAC__reads__PCA"),          out_path, file("*_pca.pdf") \
      into ATAC_reads_PCA_plots_for_merging_pdfs
    set val("ATAC__reads__correlations"), out_path, file("*_cor.pdf") \
      into ATAC_reads_correlation_plots_for_merging_pdfs

  script:
    """

    plotPCAandCorMat.sh \
      ${input_control_present} \
      ${params.blacklisted_regions} \
      ${params.deeptools__binsize_bigwig_correlation}

    """
}

Merging_pdfs_channel = Merging_pdfs_channel
  .mix(
    ATAC_reads_PCA_plots_for_merging_pdfs
    .groupTuple(by: [0, 1])
    .map{ it.flatten() }
    .map{ [ it[0], it[1], it[2..-1] ] }
    )

Merging_pdfs_channel = Merging_pdfs_channel
  .mix(
    ATAC_reads_correlation_plots_for_merging_pdfs
    .groupTuple(by: [0, 1])
    .map{ it.flatten() }
    .map{ [ it[0], it[1], it[2..-1] ] }
    )



process ATAC_QC_reads__plotting_insert_size_distribution {
  tag "${id}"

  label "picard"

  publishDir path: "${out_fig_indiv}/${out_path}/ATAC__reads__insert_size", 
             mode: "${pub_mode}"

  when: 
    do_atac

  input:
    val out_path from Channel.value('1_Preprocessing') 
    set id, file(bam) from Bam_for_plotting_inserts_distribution

  output:
    set val("ATAC__reads__insert_size"), out_path, file("*.pdf") \
      into ATAC_reads_insert_size_plots_for_merging_pdfs

  script:
    def key = id + "__insert_size"

    """

    picard -Xmx${params.memory_picard} CollectInsertSizeMetrics \
      -INPUT "${bam}" \
      -OUTPUT "${key}.txt" \
      -METRIC_ACCUMULATION_LEVEL ALL_READS \
      -Histogram_FILE "${key}.pdf" \
      -TMP_DIR .

    """
}

Merging_pdfs_channel = Merging_pdfs_channel
  .mix(
    ATAC_reads_insert_size_plots_for_merging_pdfs.groupTuple(by: [0, 1])
    )


process ATAC_QC_reads__sampling_aligned_reads {
  tag "${id}"

  label "samtools_bedtools_perl"

  when: 
    do_atac

  input:
    set id, file(bam) from Bam_for_sampling

  output:
    set id, file("*.sam") into Bam_for_estimating_library_complexity
    set id, file("*.bam"), file("*.bai") \
        into Bam_for_measuring_overlap_with_genomic_regions, 
             Sampled_aligned_reads_for_gathering_reads_stat
    set id, file("*.NB_ALIGNED_PAIRS*"), file("*.RAW_PAIRS*") \
        into Number_of_aligned_pairs_for_gathering_reads_stat

  script:
    def key = id + "__sampled"
    def new_sam = key + ".sam"
    def new_bam = key + ".bam"

    """

    # saving the header
    samtools view -F 0x904 -H ${bam} > ${new_sam}

    # sampling a certain number of reads
    samtools view -F 0x904 ${bam} \
      | shuf - -n ${params.nb_sampled_aligned_reads} >> ${new_sam}

    # conversion to bam, sorting and indexing of sampled reads
    samtools view -Sb ${new_sam} \
      | samtools sort - -o ${new_bam}
    samtools index ${new_bam}

    NB_ALIGNED_PAIRS=`samtools view -F 0x4 ${bam} | cut -f 1 | sort -T . \
                                                  | uniq | wc -l`
    RAW_PAIRS=`samtools view ${bam}| cut -f 1 | sort -T . | uniq | wc -l`

    touch \$NB_ALIGNED_PAIRS.NB_ALIGNED_PAIRS_${id}
    touch \$RAW_PAIRS.RAW_PAIRS_${id}

    """
}


// The code "-F 0x904" indicates to keep only mapped alignment (flag 0x4) that 
//    are not secondary (flag 0x100, multimappers) or 
//        supplementary (flag 0x800, chimeric entries)


process ATAC_QC_reads__measuring_overlap_with_genomic_regions {
  tag "${id}"

  label "samtools_bedtools_perl"

  when: 
    do_atac

  input:
    set id, file(bam), file(bai) \
      from Bam_for_measuring_overlap_with_genomic_regions

  output:
    set id, file("*.txt") \
      into Overlap_with_genomic_regions_results_for_gathering_reads_stat

  shell:
    '''

    id=!{id}
    bam=!{bam}
    BED_PATH=!{params.bed_regions}

    new_txt="${id}_reads_overlap_with_genomic_regions.txt"

    getTotalReadsMappedToBedFile () { 
      bedtools coverage -a $1 -b $2 | cut -f 4 \
      | awk '{ sum+=$1} END {print sum}'
      }


    PROMOTER=`getTotalReadsMappedToBedFile $BED_PATH/promoters.bed ${bam}`
    EXONS=`getTotalReadsMappedToBedFile $BED_PATH/exons.bed ${bam}`
    INTRONS=`getTotalReadsMappedToBedFile $BED_PATH/introns.bed ${bam}`
    INTERGENIC=`getTotalReadsMappedToBedFile $BED_PATH/intergenic.bed ${bam}`
    GENES=`getTotalReadsMappedToBedFile $BED_PATH/genes.bed ${bam}`
    BAM_NB_READS=`samtools view -c ${bam}`

    echo "PROMOTER EXONS INTRONS INTERGENIC GENES BAM_NB_READS" > ${new_txt}
    echo "$PROMOTER $EXONS $INTRONS $INTERGENIC $GENES $BAM_NB_READS" \
      >> ${new_txt}

    awkDecimalDivision () { 
      awk -v x=$1 -v y=$2 'BEGIN {printf "%.2f\\n", 100 * x / y }'
    }

    P_PROM=$(awkDecimalDivision $PROMOTER $BAM_NB_READS)
    P_EXONS=$(awkDecimalDivision $EXONS $BAM_NB_READS)
    P_INTRONS=$(awkDecimalDivision $INTRONS $BAM_NB_READS)
    P_INTERGENIC=$(awkDecimalDivision $INTERGENIC $BAM_NB_READS)
    P_GENES=$(awkDecimalDivision $GENES $BAM_NB_READS)
    P_READS=$(awkDecimalDivision $BAM_NB_READS $BAM_NB_READS)

    echo "$P_PROM $P_EXONS $P_INTRONS $P_INTERGENIC $P_GENES $P_READS" \
      >> ${new_txt}

    '''

}

// => preparing ATAC-Seq related statistics on aligned reads

// bedtools_coverage_only () { bedtools coverage -a $1 -b $2 ;}
// bedtools_coverage_only $BED_PATH/promoters.bed ${bam} | head
// bedtools_coverage_only $BED_PATH/promoters.bed ${bam} | cut -f 4 | head



process ATAC_QC_reads__estimating_library_complexity {
  tag "${id}"

  label "picard"

  when: 
    do_atac

  input:
    set id, file(sam) from Bam_for_estimating_library_complexity

  output:
    set id, file("*.txt") \
      into Library_complexity_for_gathering_reads_stat

  script:
    """

    picard -Xmx${params.memory_picard} EstimateLibraryComplexity \
      -INPUT ${sam} \
      -OUTPUT ${id}_library_complexity.txt

    """
}



process ATAC_QC_reads__sampling_trimmed_reads {
  tag "${id}"

  label "bbmap"

  when: 
    do_atac

  input:
    set id, file(read1), file(read2) from Trimmed_reads_for_sampling

  output:
    set id, file("*R1.fastq"), file("*R2.fastq") \
      into Sampled_trimmed_reads_for_alignment

  script:
    """

    reformat.sh \
      in1=${read1} \
      in2=${read2} \
      out1=${id}_subsampled_R1.fastq \
      out2=${id}_subsampled_R2.fastq \
      samplereadstarget=${params.nb_sampled_trimmed_reads}

    """
}



process ATAC_QC_reads__aligning_sampled_reads {
  tag "${id}"

  label "bowtie2_samtools"
  cpus params.botwie2__nb_threads

  when: 
    do_atac

  input:
    set id, file(read1), file(read2) from Sampled_trimmed_reads_for_alignment

  output:
    set id, file("*_ref.qc"), file("*_conta.qc") \
      into Flagstat_on_aligned_sampled_trimmed_reads_for_gathering_reads_stat
    set file("*_ref.bam"), file("*_conta.bam")

  script:
    def key_ref = id + "_ref"
    def new_bam_ref = key_ref + ".bam"
    def key_conta = id + "_conta"
    def new_bam_conta = key_conta + ".bam"

    """

    bowtie2 -p ${params.botwie2__nb_threads} \
      --very-sensitive \
      --end-to-end \
      --no-mixed \
      -X 2000 \
      --met 1 \
      -x "${params.bowtie2_indexes}" \
      -q -1 "${read1}" -2 "${read2}" \
      | samtools view -bS -o ${new_bam_ref} -

    samtools flagstat ${new_bam_ref} > "${key_ref}.qc"

    bowtie2 -p ${params.botwie2__nb_threads} \
      --very-sensitive \
      --end-to-end \
      --no-mixed \
      -X 2000 \
      --met 1 \
      -x "${params.bowtie2_indexes_contam}" \
      -q -1 "${read1}" -2 "${read2}" \
      | samtools view -bS -o ${new_bam_conta} -

    samtools flagstat ${new_bam_conta} > "${key_conta}.qc"

    """
}


Overlap_with_genomic_regions_results_for_gathering_reads_stat
    .join(Library_complexity_for_gathering_reads_stat)
    .join(Flagstat_on_aligned_sampled_trimmed_reads_for_gathering_reads_stat)
    .join(Sampled_aligned_reads_for_gathering_reads_stat)
    .join(Number_of_aligned_pairs_for_gathering_reads_stat)
    .join(Reads_in_bed_files_for_gathering_reads_stat)
    .set{ All_stat_for_gathering_reads_stat }



process ATAC_QC_reads__gathering_all_stat {
  tag "${id}"

  label "samtools_bedtools_perl"

  when: 
    do_atac

  input:
    set id, file(overlap_genomic_regions), file(library_complexity),
      file(cel_flagstat), file(op50_flagstat), file(bam), file(bai),
      file(nb_aligned_pairs), file(nb_total_pairs), file(final_bed) \
      from All_stat_for_gathering_reads_stat

  output:
    file("*.txt") into Stats_on_aligned_reads_for_gathering

  shell:
    '''

    id=!{id}
    overlap_genomic_regions=!{overlap_genomic_regions}
    library_complexity=!{library_complexity}
    cel_flagstat=!{cel_flagstat}
    op50_flagstat=!{op50_flagstat}
    bam=!{bam}
    nb_aligned_pairs=!{nb_aligned_pairs}
    nb_total_pairs=!{nb_total_pairs}
    final_bed=!{final_bed}

    # number of paired end reads
    FINAL_PAIRS=`awk 'END{print NR/2}' ${final_bed}`
    ALIGNED_PAIRS=`basename ${nb_aligned_pairs} .NB_ALIGNED_PAIRS_${id}`
    RAW_PAIRS=`basename ${nb_total_pairs} .RAW_PAIRS_${id}`

    # percentage of aligned reads
    PERCENT_ALIGN_REF=`sed '7q;d' ${cel_flagstat} \
    | cut -f 2 -d '(' | cut -f 1 -d '%'`
    PERCENT_ALIGN_CONTA=`sed '7q;d' ${op50_flagstat} \
    | cut -f 2 -d '(' | cut -f 1 -d '%'`

    # percentage of mitochondrial reads
    PERCENT_MITO=$(awk "BEGIN { 
      print 100 * $(samtools view -c ${bam} MtDNA) / $(samtools view -c ${bam}) 
      }" )

    # percentage of reads mapping to various annotated regions
    PERCENT_PROMOTERS=$(sed '3q;d' ${overlap_genomic_regions} | cut -f 1 -d ' ')
    PERCENT_EXONS=$(sed '3q;d' ${overlap_genomic_regions} | cut -f 2 -d ' ')
    PERCENT_INTRONS=$(sed '3q;d' ${overlap_genomic_regions} | cut -f 3 -d ' ')
    PERCENT_INTERGENIC=$(sed '3q;d' ${overlap_genomic_regions} | cut -f 4 -d ' ')
    PERCENT_GENIC=$(sed '3q;d' ${overlap_genomic_regions} | cut -f 5 -d ' ')

    # library size and percentage of duplicated reads 
    LIBRARY_SIZE=0
    RES_LIB_COMP=`sed '8q;d' ${library_complexity}`
    PERCENT_DUPLI=`echo $RES_LIB_COMP | cut -f 9 -d ' '`
    if [ "$PERCENT_DUPLI" = 0 ]; then
    LIBRARY_SIZE=0
    else
    LIBRARY_SIZE=`echo $RES_LIB_COMP | cut -f 10 -d ' '`
    fi
    PERCENT_DUPLI=`awk -v x=$PERCENT_DUPLI 'BEGIN {printf "%.2f\\n", 100 * x }' }`

    # gathering the results
    echo \
    "${id},$PERCENT_MITO,$PERCENT_PROMOTERS,$PERCENT_EXONS,$PERCENT_INTRONS,\
    $PERCENT_INTERGENIC,$PERCENT_GENIC,$PERCENT_ALIGN_REF,$PERCENT_ALIGN_CONTA,\
    $PERCENT_DUPLI,$LIBRARY_SIZE,$RAW_PAIRS,$ALIGNED_PAIRS,$FINAL_PAIRS" \
    > ${id}_bam_stats.txt

    '''
}

// library size = library complexity
// # (=> needs aligned but not filtered reads)



process ATAC_QC_reads__gathering_all_samples {

  label "samtools_bedtools_perl"

  publishDir path: "${out_tab_merge}/1_Preprocessing", mode: "${pub_mode}"
  publishDir path: "${out_tab_indiv}/1_Preprocessing", mode: "${pub_mode}"

  when: 
    do_atac

  input:
    file("*") from Stats_on_aligned_reads_for_gathering.collect()

  output:
    file("*.csv") into Bam_stat_for_splitting

  shell:
    '''

    OUTFILE="ATAC__alignment_statistics.csv"

    echo \
      "LIBRARY_NAME,PERCENT_MITO,PERCENT_PROMOTERS,PERCENT_EXONS,\
PERCENT_INTRONS,PERCENT_INTERGENIC,PERCENT_GENIC,PERCENT_ALIGN_REF,\
PERCENT_ALIGN_CONTA,PERCENT_DUPLI,LIBRARY_SIZE,RAW_PAIRS,ALIGNED_PAIRS,\
FINAL_PAIRS" \
      > ${OUTFILE}

    cat *_bam_stats.txt >> ${OUTFILE}

   '''
}

// => gathering all individual statistics files and making a merged table 



process ATAC_QC_reads__splitting_stat_for_multiqc {

  label "r_basic"

  when: 
    do_atac

  input:
    file(bam_stat_csv) from Bam_stat_for_splitting

  output:
    file("*") into Bam_stat_for_multiqc

  shell:
    '''
    #!/usr/bin/env Rscript

    library(magrittr)
    library(dplyr)

    bam_stat_csv = '!{bam_stat_csv}'


    df = read.csv(bam_stat_csv, stringsAsFactors = F)
    colnames(df) %<>% tolower %>% gsub('percent', 'percentage', .)

    df %<>% rename(
      percentage_mitochondrial       = percentage_mito, 
      percentage_align_reference     = percentage_align_ref, 
      percentage_aligned_contaminant = percentage_align_conta, 
      percentage_duplications        = percentage_dupli
      )

    colnames = colnames(df)[2:ncol(df)]
    for(colname in colnames){
      df1 = df[, c('library_name', colname)]
      write.csv(df1, row.names = F, file = paste0(colname, '_mqc.csv'))
    }

   '''
}


process ATAC_QC_reads__running_multiQC {

  label "multiqc"

  publishDir path: "${out_fig_indiv}/1_Preprocessing", 
             mode: "${pub_mode}", pattern: "*.html"

  publishDir path: "${out_fig_merge}/1_Preprocessing", 
             mode: "${pub_mode}", pattern: "*.html"

  when: 
    do_atac

  input:
    file ('fastqc/*') from ATAC_fastqc_reports_for_multiqc.flatten().toList()
    file(csv_files) from Bam_stat_for_multiqc

  output:
    file "ATAC__multiQC.html"

  script:
    """

    multiqc -f .
    mv multiqc_report.html ATAC__multiQC.html

    """
}

// This kind of commands could maybe make a more pretty multiqc report, 
// but I would need to investigate more on how to do that exactly
// FILENAME=`basename bam_stats.csv .csv`
// cut -d "," -f 1-2 bam_stats.csv > ${FILENAME}_test_mqc.csv
//
// FILENAME=`basename ${bam_stat} .csv`
// cut -d "," -f 1-6 ${bam_stat} > \${FILENAME}_percentage_mqc.csv
// cut -d "," -f 1,7-10 ${bam_stat} > \${FILENAME}_counts_mqc.csv

// the raw multiqc_data folder is probably not useful to publish. 
// Interested people may add it back though by uncommenting these 2 lines.






process ATAC_peaks__calling_peaks {
  tag "${id}"

  label "macs2"

  publishDir \
    path: "${out_processed}/1_Preprocessing/ATAC__peaks__raw", 
    mode: "${pub_mode}", enabled: save_all_bed

  when: 
    do_atac

  input:
    set id, file(bed) from Reads_in_bed_files_for_calling_peaks

  output:
    set id, file("*.narrowPeak") into Peaks_for_sub_peaks_calling

  script:
    """

    export TMPDIR="." PYTHON_EGG_CACHE="."

    macs2 callpeak \
        --treatment "${bed}" \
        --format BED \
        --name "${id}__macs2" \
        --qvalue "${params.macs2__qvalue}" \
        --gsize "${params.effective_genome_size}" \
        --nomodel \
        --shift -75 \
        --extsize 150 \
        --keep-dup all \
        --call-summits

    """
}


//// doc
// https://pypi.org/project/MACS2/
// https://pypi.org/project/MACS3/
// https://github.com/macs3-project/MACS
// https://github.com/macs3-project/MACS
// https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md
// https://macs3-project.github.io/MACSr/articles/MACSr.html

//// other relevant links 
// https://github.com/macs3-project/MACS/discussions/435
// https://github.com/macs3-project/MACS/issues/145
// https://twitter.com/XiChenUoM/status/1336658454866325506
// https://groups.google.com/g/macs-announcement/c/4OCE59gkpKY/m/v9Tnh9jWriUJ
// https://www.biostars.org/p/209592/
// https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/04_peak_calling_macs.html
// https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html



// macs2 manual: https://github.com/jsh58/MACS

//// Here are details for the macs2 parameters (the macs2 website is broken)
// https://www.biostars.org/p/207318/
// --extsize
//    While '--nomodel' is set, MACS uses this parameter to extend reads in 
//    5'->3' direction to fix-sized fragments. For example, if the size of 
//    binding region for your transcription factor is 200 bp, and you want to 
//    bypass the model building by MACS, this parameter can be set as 200. 
//    This option is only valid when --nomodel is set or when MACS fails to 
//    build model and --fix-bimodal is on.
// 
// --shift
//     Note, this is NOT the legacy --shiftsize option which is replaced by 
//    --extsize! You can set an arbitrary shift in bp here. Please Use 
//    discretion while setting it other than default value (0). When --nomodel 
//    is set, MACS will use this value to move cutting ends (5') then apply 
//    --extsize from 5' to 3' direction to extend them to fragments. When this 
//    value is negative, ends will be moved toward 3'->5' direction, otherwise 
//    5'->3' direction. Recommended to keep it as default 0 for ChIP-Seq 
//    datasets, or -1 * half of EXTSIZE together with --extsize option for 
//    detecting enriched cutting loci such as certain DNAseI-Seq datasets. 
//    Note, you can't set values other than 0 if format is BAMPE or BEDPE 
//    for paired-end data. Default is 0.

//// a recommanded setting to use both read pairs is the following:
// https://github.com/macs3-project/MACS/issues/145#issuecomment-742593158
// https://twitter.com/XiChenUoM/status/1336658454866325506
// macs2 callpeak -f BED --shift -100 --extsize 200 ...


//// There is also a justification to use 150 bp instead of 200:
// https://twitter.com/hoondy/status/1337170830997004290
// Is there any reason you chose the size as 200bp? For example, I use 
//    "shift -75 --extsize 150" (approx. size of nucleosome) and wonder if there 
//    were any benefit of making the unit of pileup 200bp.
// We use 200 due to habit. The choice is arbitrary. Not sure how diff it makes. 
//    Intuitively smaller fragLen give better res. but it can’t be too short. 
//    Check Fig1 in this ref (not ATAC though) 
//    https://ncbi.nlm.nih.gov/pmc/articles/PMC2596141/ 75 gives better res. 
//    Maybe that’s the reason ENCODE use it.

//// So to sum up: 
// - the two ends of each read pair are analyzed separately
// - we keep only the 5' end of each read pair. 
//    This is were we want our peak to be centered
// - we shift the reads to take into account the transposase 
//    (+4 bp on + strand, -5 on the - strand)
// - we shift the reads manually by -75 bp (in the 5' direction)
// - we use macs2 with the argument --extsize 150 which indicate to extend the 
//    read by 150 bp in the 3' direction. This way they will be centered in the 
//    original 5' location. 150 bp is the approximate size of nucleosome.
// - the peaks in the .narrowPeak file are 150 base pair long


//// methods from the Daugherty 2017 paper (Brunet lab)
// Chromatin accessibility dynamics reveal novel functional enhancers in 
//    C. elegans
// https://genome.cshlp.org/content/27/12/2096.long
// Supplemental_Material.pdf, section ATAC-seq peak calling
// For every replicate, prior to calling peaks with MACS (Zhang et al. 2008) 
//    (v2.1), single-base 
// reads were shifted 75bp 5’ to mimic read distributions of a 150 bp fragment 
//    of ChIP-seq, thereby
// allowing use of MACS. The following settings were used for MACS: 
//    -g 9e7, -q 5e-2, --nomodel, --extsize 150, -B, --keep-dup all, and 
//    --call-summits.  
// The resulting peaks included a portion that had multiple summits; to 
// maximize our resolution these summits were subsequently separated into 
// individual peaks by treating the midpoint between the two adjacent summits 
// as the 5’ and 3’ end of each individual peak, respectively.
// macs2.1 doc: https://pypi.org/project/MACS2/2.1.0.20151222/


//// Calling macs2 without the --call-summits option
// narrowPeaks file
// I       3609    4118    ctl_2_macs2_peak_1      30998   .       25.663  3103.76 3099.82 325
// I       11019   11724   ctl_2_macs2_peak_2      2539    .       6.18189 255.951 253.917 310
// I       12930   13338   ctl_2_macs2_peak_3      185     .       1.97507 19.9684 18.5887 197
// I       14861   15110   ctl_2_macs2_peak_4      159     .       1.91822 17.3277 15.9723 113
// I       15340   15573   ctl_2_macs2_peak_5      594     .       2.95858 61.078  59.4709 119
// I       16735   17315   ctl_2_macs2_peak_6      6693    .       10.886  671.782 669.305 306
// I       21631   21997   ctl_2_macs2_peak_7      121     .       1.84936 13.4594 12.1469 260
// I       22239   22394   ctl_2_macs2_peak_8      149     .       1.95444 16.3335 14.988  116
// I       24470   24626   ctl_2_macs2_peak_9      230     .       2.21713 24.4981 23.0807 65
// I       26817   27032   ctl_2_macs2_peak_10     1132    .       4.11903 115.028 113.257 85
//
// Summits file
// I       3934    3935    ctl_2_macs2_peak_1      3099.82
// I       11329   11330   ctl_2_macs2_peak_2      253.917
// I       13127   13128   ctl_2_macs2_peak_3      18.5887
// I       14974   14975   ctl_2_macs2_peak_4      15.9723
// I       15459   15460   ctl_2_macs2_peak_5      59.4709
// I       17041   17042   ctl_2_macs2_peak_6      669.305
// I       21891   21892   ctl_2_macs2_peak_7      12.1469
// I       22355   22356   ctl_2_macs2_peak_8      14.988
// I       24535   24536   ctl_2_macs2_peak_9      23.0807
// I       26902   26903   ctl_2_macs2_peak_10     113.257

// 
//
//// Calling macs2 with the --call-summits option
// narrowPeaks file
// I       3609    4118    ctl_2_macs2_peak_1      29989   .       25.0793 3002.86 2998.99 319
// I       11019   11724   ctl_2_macs2_peak_2      2539    .       6.18189 255.951 253.917 310
// I       12930   13338   ctl_2_macs2_peak_3      168     .       1.92215 18.1832 16.8191 210
// I       14861   15110   ctl_2_macs2_peak_4      129     .       1.81726 14.2625 12.9413 127
// I       15340   15573   ctl_2_macs2_peak_5      437     .       2.63535 45.2707 43.7312 108
// I       16735   17315   ctl_2_macs2_peak_6      5854    .       10.0244 587.809 585.402 294
// I       21631   21997   ctl_2_macs2_peak_7a     69      .       1.6287  8.22598 6.995   77
// I       21631   21997   ctl_2_macs2_peak_7b     121     .       1.84936 13.4594 12.1469 261
// I       22239   22394   ctl_2_macs2_peak_8      100     .       1.7653  11.3334 10.0498 85
// I       24470   24626   ctl_2_macs2_peak_9      223     .       2.19612 23.7962 22.3841 72
// Summits file
// I       3928    3929    ctl_2_macs2_peak_1      2998.99
// I       11329   11330   ctl_2_macs2_peak_2      253.917
// I       13140   13141   ctl_2_macs2_peak_3      16.8191
// I       14988   14989   ctl_2_macs2_peak_4      12.9413
// I       15448   15449   ctl_2_macs2_peak_5      43.7312
// I       17029   17030   ctl_2_macs2_peak_6      585.402
// I       21708   21709   ctl_2_macs2_peak_7a     6.995
// I       21892   21893   ctl_2_macs2_peak_7b     12.1469
// I       22324   22325   ctl_2_macs2_peak_8      10.0498
// I       24542   24543   ctl_2_macs2_peak_9      22.3841



// we can see here that when using the --call-summits option, the peak number 7 
//    (summits at coordinates 21891) was split into two different peaks (summits 
//    at coordinates 21708 and 21892)
// however the narrowPeaks file is just duplicated



process ATAC_peaks__splitting_multi_summits_peaks {
  tag "${id}"

  label "samtools_bedtools_perl"

  publishDir path: "${out_processed}/1_Preprocessing/ATAC__peaks__split", 
    mode: "${pub_mode}", enabled: save_all_bed

  when: 
    do_atac

  input:
    set id, file(peaks) from Peaks_for_sub_peaks_calling

  output:
    set id, file("*.narrowPeak") into Peaks_for_removing_blacklisted_regions

  script:
    """

    perl "${projectDir}/bin/splitMACS2SubPeaks.pl" "${peaks}" \
      > "${id}__split_peaks.narrowPeak"

    """
}

// This scripts comes from Anshul Kundaje, and was made for the Daugherty et al 
//    2017 manuscript (Brunet lab)

// Example:

//// raw narrowPeak file
// I       3609    4118    ctl_2_macs2_peak_1      29989   .       25.0793 3002.86 2998.99 319
// I       11019   11724   ctl_2_macs2_peak_2      2539    .       6.18189 255.951 253.917 310
// I       12930   13338   ctl_2_macs2_peak_3      168     .       1.92215 18.1832 16.8191 210
// I       14861   15110   ctl_2_macs2_peak_4      129     .       1.81726 14.2625 12.9413 127
// I       15340   15573   ctl_2_macs2_peak_5      437     .       2.63535 45.2707 43.7312 108
// I       16735   17315   ctl_2_macs2_peak_6      5854    .       10.0244 587.809 585.402 294
// I       21631   21997   ctl_2_macs2_peak_7a     69      .       1.6287  8.22598 6.995   77
// I       21631   21997   ctl_2_macs2_peak_7b     121     .       1.84936 13.4594 12.1469 261
// I       22239   22394   ctl_2_macs2_peak_8      100     .       1.7653  11.3334 10.0498 85
// I       24470   24626   ctl_2_macs2_peak_9      223     .       2.19612 23.7962 22.3841 72
// 
//// splitted narrowPeak file
// I       3609    4118    ctl_2_macs2_peak_1      29989   .       25.0793 3002.86 2998.99 319
// I       11019   11724   ctl_2_macs2_peak_2      2539    .       6.18189 255.951 253.917 310
// I       12930   13338   ctl_2_macs2_peak_3      168     .       1.92215 18.1832 16.8191 210
// I       14861   15110   ctl_2_macs2_peak_4      129     .       1.81726 14.2625 12.9413 127
// I       15340   15573   ctl_2_macs2_peak_5      437     .       2.63535 45.2707 43.7312 108
// I       16735   17315   ctl_2_macs2_peak_6      5854    .       10.0244 587.809 585.402 294
// I       21631   21800   ctl_2_macs2_peak_7a     69      .       1.6287  8.22598 6.995   77
// I       21800   21997   ctl_2_macs2_peak_7b     121     .       1.84936 13.4594 12.1469 92
// I       22239   22394   ctl_2_macs2_peak_8      100     .       1.7653  11.3334 10.0498 85
// I       24470   24626   ctl_2_macs2_peak_9      223     .       2.19612 23.7962 22.3841 72

// here we see that peaks peaks 7a and 7b have the same start (21631) and end 
//    (21997) but different summits after splitting, they are split into two 
//    regions in their middle; so 7a ends at 21800, and 7b starts at 21800


process ATAC_peaks__removing_blacklisted_regions {
  tag "${id}"

  label "samtools_bedtools_perl"

  publishDir path: "${out_processed}/1_Preprocessing/ATAC__peaks__split__no_BL", 
             mode: "${pub_mode}", enabled: save_all_bed

  when: 
    do_atac

  input:
    set id, file(peaks) from Peaks_for_removing_blacklisted_regions

  output:
    file("*.bed")
    set id, file("*__peaks_kept_after_blacklist_removal.bed") \
      into Peaks_without_blacklist_1, 
           Macs2_Peaks_without_blacklist_1_for_annotating_them

  script:
    """

    intersectBed -v -a "${peaks}" -b "${params.blacklisted_regions}" \
      > "${id}__peaks_kept_after_blacklist_removal.bed"

    intersectBed -u -a "${peaks}" -b "${params.blacklisted_regions}" \
      > "${id}__peaks_lost_after_blacklist_removal.bed"

    """
}

// there are two output channels because:
// all peaks including input control are sent for annotation by ChipSeeker 
//    to get the distribution of peaks in various genomic locations
// the peaks are then sent for input control removal before DiffBind analysis


// println "params.use_input_control: ${params.use_input_control}"

Peaks_without_blacklist_1
  .dump(tag:'peaks_wo_bl') {"Peaks without blacklisted regions: ${it}"}
  .branch {
    with_input_control: params.use_input_control
    without_input_control: true
  }
  .set { Peaks_without_blacklist_2 }

// this command redirects the channel items to: 
// Peaks_without_blacklist_2.with_input_control   
//     if params.use_input_control = true
// Peaks_without_blacklist_2.without_input_control 
//    if params.use_input_control = false


Peaks_without_blacklist_2.with_input_control
  .dump(tag:'peaks_wo_bl_w_ic')
  .branch { it ->
    // control: it[0].split("_")[0]== 'input'
    control: it[0] == 'input'
    treatment: true
  }
  .set { Peaks_without_blacklist_3 }


Peaks_without_blacklist_3.treatment
  .combine(Peaks_without_blacklist_3.control)
  .dump(tag:'peaks_input') {"Peaks with input_control controls: ${it}"}
  .map { it[0, 1, 3] }
  .set { Peaks_treatment_with_control }


process ATAC_peaks__removing_input_control_peaks {
  tag "${id}"

  label "samtools_bedtools_perl"

  publishDir \
    path: "${out_processed}/1_Preprocessing/ATAC__peaks__split__no_BL_input", 
    mode: "${pub_mode}", enabled: save_all_bed

  when: 
    do_atac

  input:
    set id, file(peaks), file(input_control_peaks) \
      from Peaks_treatment_with_control

  output:
    file("*.bed")
    set id, file("*__peaks_kept_after_input_control_removal.bed") \
      into Peaks_for_removing_specific_regions_1

  script:
    """

    intersectBed -wa -v -f ${params.input_control_overlap_portion} \
      -a "${peaks}" -b "${input_control_peaks}" \
      > "${id}__peaks_kept_after_input_control_removal.bed"

    intersectBed -wa -u -f ${params.input_control_overlap_portion} \
      -a "${peaks}" -b "${input_control_peaks}" \
      > "${id}__peaks_lost_after_input_control_removal.bed"

    """
}

// https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
// option   Description
// -wa 	  	Write the original entry in A for each overlap.
// -f  	  	Minimum overlap required as a fraction of A. - Default is 1E-9 
//          (i.e., 1bp). - FLOAT (e.g. 0.50)
// -v 	  	Only report those entries in A that have _no overlaps_ with B. 
//          - Similar to ‘grep -v’ (an homage).
// -u 	    Write the original A entry _once_ if _any_ overlaps found in B. 
//          - In other words, just report the fact >=1 hit was found. 
//          - Overlaps restricted by -f and -r.

Peaks_without_blacklist_2.without_input_control
  .concat(Peaks_for_removing_specific_regions_1)
  .set{   Peaks_for_removing_specific_regions_2 }





//// DIFFENRENTIAL BINDING

comparisons_files = file(params.design__comparisons)
Channel
  .from(comparisons_files.readLines())
  .map {
          m = it.split()
          condition1 = m[0]
          condition2 = m[1]
          [ condition1, condition2 ]
       }
  .dump(tag:'comp_file') {"comparison file: ${it}"}
  .into { comparisons_files_for_merging; comparisons_files_for_mRNA_Seq }

Reads_in_bam_files_for_diffbind
  .tap{ Reads_input_control }
  .join( Peaks_for_removing_specific_regions_2, remainder: true )
  .map { [ it[0].split("_")[0], it[0..-1]] }
  .groupTuple()
  .dump(tag:'reads_peaks') {"merged reads and peaks: ${it}"}
  .into { reads_and_peaks_1 ; reads_and_peaks_2 ; reads_and_peaks_3 }

regions_to_remove = Channel.fromPath(params.design__regions_to_remove)

comparisons_files_for_merging
  .combine(reads_and_peaks_1)
  .combine(reads_and_peaks_2)
  .filter { id_comp_1, id_comp_2, 
              id_1, reads_and_peaks_1, 
              id_2, reads_and_peaks_2 -> 
                id_comp_1 == id_1 && id_comp_2 == id_2 }
  .map { 
    [ 
      it[0] + '_vs_' + it[1], 
      it.flatten().findAll { it =~ '\\.bed' }, 
      it.flatten().findAll { it =~ "\\.bam" } 
      ] 
    }
  .tap { Reads_in_bam_files_for_diffbind_1 }
  .map { it[0,1] }
  .combine(regions_to_remove)
  .dump(tag:'clean_peaks') {"peaks for removing regions: ${it}"}
  .set { Peaks_for_removing_specific_regions_3 }



process ATAC_peaks__removing_specific_regions {
  tag "${COMP}"

  label "samtools_bedtools_perl"

  publishDir \
   path: "${out_processed}/1_Preprocessing/ATAC__peaks__split__no_BL_input_RNAi", 
   mode: "${pub_mode}", enabled: save_last_bed

  when: 
    do_atac

  input:
    set COMP, file(bed_files), regions_to_remove \
      from Peaks_for_removing_specific_regions_3

  output:
    file("*.bed")
    set COMP, file("*__peaks_kept_after_specific_regions_removal.bed") \
      into Peaks_for_diffbind

  shell:
    '''

    COMP=!{COMP}
    RTR=!{regions_to_remove}
    BED_FILES="!{bed_files}"

    COMP1=${COMP/_vs_/|}
    ( echo $COMP1 | grep -E -f - $RTR || true ) > rtr_filtered.txt

    if [ ! -s rtr_filtered.txt ]; then
      for FILE in ${BED_FILES}
      do
        id=${FILE/__*/}
        cp $FILE "${id}__peaks_kept_after_specific_regions_removal.bed"
      done
    else
      cat rtr_filtered.txt | sed "s/,//g" | sed "s/.*->//g" \
        | sed "s/[:-]/\\t/g" | sed "s/chr//g" > rtr_filtered_formatted.txt 

      for FILE in ${BED_FILES}
      do
        id=${FILE/__*/}
        intersectBed -v -a $FILE -b rtr_filtered_formatted.txt \
          > "${id}__peaks_kept_after_specific_regions_removal.bed"
        intersectBed -u -a $FILE -b rtr_filtered_formatted.txt \
          > "${id}__peaks_lost_after_specific_regions_removal.bed"
      done
    fi

    '''
}

// really weird nextflow bug: this command makes all commands below ignored 
// when the rtr file is empty:
// echo $COMP1 | grep -E -f - $RTR > rtr_filtered.txt
// https://github.com/nextflow-io/nextflow/issues/1149
// => solution: grep pattern || true
// or: process.shell = ['/bin/bash','-u']


// cat rtr_filtered.txt | sed "s/,//g" | sed "s/.*->//g" | sed "s/[:-]/\\t/g" \
//   | sed "s/chr//g" > rtr_filtered_formatted.txt 
// echo $COMP2 | grep -f - $RTR >> rtr_filtered.txt


// note: the reason why this process is here and not upstread is because we 
// want to remove in all bed files the peaks that are in specific regions 
// (i.e. RNAi) that we want to avoid. This is because, Diffbind will consider 
// all peaks for his analysis, so if we remove one such peak in just one of the 
// two samples to compare, if it is still present in the other sample then it 
// will be included in the analysis and it will likely be found as differential 
// bound during the DBA. I.e. we compare daf-16RNAi vs control. if there is a 
// macs2 peak at daf-16 in the control condition, then even if we remove this 
// peak in the daf-16RNAi condition, it will be included in the final analysis.



Reads_input_control
  // .filter{ id, bam_files -> id.split('_')[0] == 'input'}
  .filter{ id, bam_files -> id == 'input'} 
  // only 1 replicate is allowed for now and it IC should have the name "input"
  .set{ Reads_input_control_1 }

Reads_in_bam_files_for_diffbind_1
  .map { it[0,2] }
  .join(Peaks_for_diffbind)
  .dump(tag:'bam_bai') {"bam and bai files: ${it}"}
  .branch {
    with_input_control: params.use_input_control
    without_input_control: true
  }
  .set { Reads_and_peaks_for_diffbind_1 }

Reads_and_peaks_for_diffbind_1.with_input_control
  .dump(tag:'reads_peaks_w_ic')
  .combine(Reads_input_control_1)
  .map{ 
          comp_id, bam_files, bed_files, input_id, imput_bam -> 
        [ comp_id, bed_files, [bam_files, imput_bam].flatten() ]
      }
  // .map{ it -> [ it[0], it[1], [it[2,4].flatten()]] }
  .set{ Reads_and_peaks_with_input_control }

Reads_and_peaks_for_diffbind_1.without_input_control
  .concat(Reads_and_peaks_with_input_control)
  .dump(tag:'input_diffbind') {"reads and peaks for diffbind: ${it}"}
  .set{ Reads_and_peaks_for_diffbind_2 }






process ATAC_QC_peaks__computing_and_plotting_saturation_curve {
  tag "${id}"

  label "macs2"

  publishDir path: "${out_fig_indiv}/${out_path}/ATAC__peaks__saturation_curve", 
             mode: "${pub_mode}"

  when: 
    do_atac & params.do_saturation_curve

  input:
    val out_path from Channel.value('1_Preprocessing') 
    set id, file(bed) \
      from Reads_in_bed_files_for_computing_and_plotting_saturation_curve

  output:
    set val("ATAC__peaks__saturation_curve"), out_path, file('*.pdf') \
      into ATAC_saturation_curve_plots_for_merging_pdfs

  script:
    """

    export TMPDIR="." PYTHON_EGG_CACHE="."

    PERCENTS=`seq 10 10 100`

    for PERCENT in \${PERCENTS}
    do
      BED_FILE="${id}_sampled_\${PERCENT}_percent.bed"

      macs2 randsample -i ${bed} \
        --seed 38 \
        -o \${BED_FILE} \
        --percentage \${PERCENT} \
        -f BED

      macs2 callpeak \
        --treatment \${BED_FILE} \
        --format BED \
        --name \${BED_FILE}__macs2 \
        --qvalue "${params.macs2__qvalue}" \
        --gsize "${params.effective_genome_size}" \
        --nomodel \
        --shift -75 \
        --extsize 150 \
        --keep-dup all \
        --call-summits
    done

    Rscript "${projectDir}/bin/plot_saturation_curve.R" ${id}

    """
}

Merging_pdfs_channel = Merging_pdfs_channel
  .mix(ATAC_saturation_curve_plots_for_merging_pdfs.groupTuple(by: [0, 1]))




process ATAC_QC_peaks__annotating_macs2_peaks {

  tag "${id}"

  label "bioconductor"

  publishDir path: "${out_processed}/1_Preprocessing/ATAC__peaks__annotated_rds", 
             mode: "${pub_mode}"

  when: 
    do_atac & params.do_raw_peak_annotation

  input:
    set id, file(peaks_bed_file) \
      from Macs2_Peaks_without_blacklist_1_for_annotating_them

  output: 
    set id, file("*.rds") \
      into Annotated_macs2_peaks_for_plotting_each_sample, 
           Annotated_macs2_peaks_for_plotting_all_samples_grouped_1 \
           optional true

  shell:
    '''
    #!/usr/bin/env Rscript

    library(ChIPseeker)
    library(magrittr)
    id             = '!{id}'
    peaks_bed_file = '!{peaks_bed_file}'
    upstream       = !{params.macs2_peaks__promoter_up}
    downstream     = !{params.macs2_peaks__promoter_down}
    tx_db          = AnnotationDbi::loadDb('!{params.txdb}') 

    check_upstream_and_downstream_1 <- function (upstream, downstream){
        if (class(upstream) != class(downstream)) {
            stop("the type of upstream and downstream should be the same...")
        }
        if (!is.numeric(upstream) && !is.null(upstream)) {
            stop("upstream and downstream parameter should be numeric or NULL...")
        }
        if (inherits(upstream, "rel")) {
            if (as.numeric(upstream) < 0 || as.numeric(upstream) >
                1) {
                stop("the value of rel object should be in (0,1)...")
            }
        }
        if (is.numeric(upstream) && !inherits(upstream, "rel")) {
            if (upstream < 1 | downstream < 1) {
                stop("if upstream or downstream is integer, the value of it should be greater than 1...")
            }
        }
    }

    assignInNamespace("check_upstream_and_downstream", check_upstream_and_downstream_1, ns="ChIPseeker")

    nb_of_peaks = 
      system(paste('wc -l', peaks_bed_file), intern = T) %>% 
      gsub(' .*', '', .)
    if(nb_of_peaks == 0) {quit('no')}
    peaks = readPeakFile('!{peaks_bed_file}')

    promoter <- getPromoters(TxDb       = tx_db, 
                             upstream   = upstream, 
                             downstream = downstream)

    tag_matrix = getTagMatrix(peaks, windows = promoter)

    annotated_peaks = annotatePeak(peaks, 
                                  TxDb      = tx_db, 
                                  tssRegion = c(-upstream, downstream), 
                                  level     = 'gene', 
                                  overlap   = 'all')

    genes = as.data.frame(annotated_peaks)$geneId

    lres = list(
      id = id,
      peaks = peaks,
      tag_matrix = tag_matrix,
      annotated_peaks = annotated_peaks,
      genes = genes
    )

    saveRDS(lres, file = paste0(id, '__annotated_peaks.rds'))

    '''
}

// ChIPseeker v1.30.0 the version present in the container, has a wrong internal function.
// For this reason we had to overwrite it with the correct function (from v1.30.3) 
// Here are the details:

// ChIPseeker v1.30.3 
// ChIPseeker:::check_upstream_and_downstream
// function (upstream, downstream)
// {
//     if (class(upstream) != class(downstream)) {
//         stop("the type of upstream and downstream should be the same...")
//     }
// 
// ChIPseeker v1.30.0
// ChIPseeker:::check_upstream_and_downstream
// function (upstream, downstream)
// {
//     if (!identical(upstream, downstream)) {
//         stop("the type of upstream and downstream should be the same...")
//     }


// these are ATAC-Seq peaks from macs

// defaults parameters for the annotation function
// 1 function (peak, tssRegion = c(-3000, 3000), TxDb = NULL, 
// level = "transcript", assignGenomicAnnotation = TRUE, 
// genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", 
// "Downstream", "Intergenic"), annoDb = NULL, addFlankGeneInfo = FALSE, 
// flankDistance = 5000, sameStrand = FALSE, ignoreOverlap = FALSE, 
// ignoreUpstream = FALSE, ignoreDownstream = FALSE, overlap = "TSS", 
// verbose = TRUE)

// the addFlankGeneInfo gives rather confusing output so we ignore it
  // geneId          transcriptId distanceToTSS   flank_txIds    flank_gene

// select(txdb, keys = keys(txdb)[1:1], columns = columns(txdb), 
// keytype = 'GENEID')


Annotated_macs2_peaks_for_plotting_all_samples_grouped_1
  .map { it[1] }
  // .groupTuple () // legacy: before there was the option type that was 
  // equal to either "raw" or "norm", so we would group tupple this way. 
  // Now there is only "raw", so we just collect peaks
  .collect()
  .dump(tag:'anno_list') {"annotated peaks as list: ${it}"}
  .set { Annotated_macs2_peaks_for_plotting_all_samples_grouped_2 }




process ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_each_sample {
  tag "${id}"

  label "bioconductor"

  publishDir \
    path: "${out_fig_indiv}/${out_path}", mode: "${pub_mode}", 
    saveAs: {
           if (it.indexOf("_coverage.pdf") > 0)        
              "ATAC__peaks__coverage/${it}"
      else if (it.indexOf("_average_profile.pdf") > 0) 
              "ATAC__peaks__average_profile/${it}"
  }

  when: 
    do_atac

  input:
    val out_path from Channel.value('1_Preprocessing') 
    set id, file(annotated_peaks_objects_rds) \
      from Annotated_macs2_peaks_for_plotting_each_sample

  output:
    set val("ATAC__peaks__coverage"), out_path, file("*_coverage.pdf") \
      into ATAC_peaks_coverage_for_merging_pdfs
    set val("ATAC__peaks__average_profile"), out_path, 
      file("*_average_profile.pdf") \
      into ATAC_peaks_average_profile_for_merging_pdfs

  shell:
    '''
    #!/usr/bin/env Rscript

    library(ChIPseeker)
    library(ggplot2)

    id = '!{id}'
    upstream = !{params.macs2_peaks__promoter_up}
    downstream = !{params.macs2_peaks__promoter_down}
    lres = readRDS('!{annotated_peaks_objects_rds}')

    pdf(paste0(id, '__peaks_coverage.pdf'))
      covplot(lres$peaks, weightCol="V5") + ggtitle(id) + 
        theme(plot.title = element_text(hjust = 0.5))
    dev.off()

    pdf(paste0(id, '__average_profile.pdf'))
      plotAvgProf(lres$tag_matrix, xlim=c(-upstream, downstream), 
        xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency") + 
        ggtitle(id) + theme(plot.title = element_text(hjust = 0.5))
    dev.off()

    '''
}

Merging_pdfs_channel = Merging_pdfs_channel
  .mix(ATAC_peaks_coverage_for_merging_pdfs.groupTuple(by: [0, 1]))

Merging_pdfs_channel = Merging_pdfs_channel
  .mix(ATAC_peaks_average_profile_for_merging_pdfs.groupTuple(by: [0, 1]))




process ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_all_samples_grouped {

  label "bioconductor"

  publishDir path: "${out_fig_indiv}/${out_path}/ATAC__peaks__grouped_plots", 
             mode: "${pub_mode}"
  publishDir path: "${out_fig_merge}/${out_path}", mode: "${pub_mode}"

  when: 
    do_atac

  input:
    val out_path from Channel.value('1_Preprocessing') 
    file ("*") from Annotated_macs2_peaks_for_plotting_all_samples_grouped_2

  output:
    file("*.pdf")

  shell:
    '''
    #!/usr/bin/env Rscript

    upstream = !{params.macs2_peaks__promoter_up}
    downstream = !{params.macs2_peaks__promoter_down}

    library(ChIPseeker)
    library(ggplot2)
    library(purrr)

    rds_files = list.files(pattern = '*.rds')
    lres = lapply(rds_files, readRDS)

    names0 = map_chr(lres, 'id')
    tag_matrix_list = map(lres, 'tag_matrix') %>% setNames(., names0)
    annotated_peaks_list = map(lres, 'annotated_peaks') %>% setNames(., names0)

    size_facet = ifelse(length(tag_matrix_list) > 12, 2.8, 5.5)

    p1 = plotAvgProf(tag_matrix_list, c(-upstream, downstream), facet='row')
    p1 = p1 + theme(axis.text.y = element_text(size = 3.5), 
          strip.text.y = element_text(size = size_facet))

    pdf('ATAC__peaks__grouped__average_profile.pdf', font = 'mono')
      print(p1)
    dev.off()

    pdf('ATAC__peaks__grouped___annotation_barplot.pdf')
      plotAnnoBar(annotated_peaks_list)
    dev.off()

    pdf('ATAC__peaks__grouped___distance_to_TSS.pdf')
      plotDistToTSS(annotated_peaks_list)
    dev.off()

    '''
}




process MRNA__quantifying_transcripts_abundances {
  tag "${id}"

  label "kallisto"
  cpus params.kallisto__nb_threads
  
  publishDir path: "${out_processed}/1_Preprocessing/mRNA__kallisto_output", 
             mode: "${pub_mode}"

  when: 
    do_mRNA

  input:
    set id, file(reads) from MRNA_reads_for_kallisto

  output:
    set id, file("kallisto_${id}") into Kallisto_results_for_sleuth_1

  script:

    def single = reads instanceof Path

    if( single ) {
      """
        mkdir kallisto_${id}
        kallisto quant \
          --single \
          -l ${params.kallisto__fragment_len} \
          -s ${params.kallisto__fragment_sd} \
          -b ${params.kallisto__bootstrap} \
          -t ${params.kallisto__nb_threads} \
          -i ${params.kallisto_transcriptome} \
          -o kallisto_${id} \
          ${reads}
      """
    }
    else {
      """
        mkdir kallisto_${id}
        kallisto quant \
          -b ${params.kallisto__bootstrap} \
          -t ${params.kallisto__nb_threads} \
          -i ${params.kallisto_transcriptome} \
          -o kallisto_${id} \
          ${reads}
      """
    }
}

// note: I adapted a basic mRNA-Seq pipeline from 
//                    https:////github.com/cbcrg/kallisto-nf



process MRNA_QC__running_fastqc {
  tag "${id}"

  label "fastqc"
  cpus params.fastqc__nb_threads

  publishDir path: "${out_processed}/1_Preprocessing/mRNA__fastqc", 
             mode: "${pub_mode}", pattern: "*.html" 

  when: 
    do_mRNA

  input:
    set id, file(reads) from MRNA_reads_for_running_fastqc

  output:
    file("*.{zip, html}") into MRNA_fastqc_reports_for_multiqc

  script:
    """

    fastqc -t ${params.fastqc__nb_threads} ${reads}

    """
}



process MRNA_QC__running_MultiQC {

  label "multiqc"

  publishDir path: "${out_fig_indiv}/1_Preprocessing", 
             mode: "${pub_mode}", pattern: "*.html"

  publishDir path: "${out_fig_merge}/1_Preprocessing", 
             mode: "${pub_mode}", pattern: "*.html"

  when: 
    do_mRNA

  input:
    file ('fastqc/*') from MRNA_fastqc_reports_for_multiqc.flatten().toList()

  output:
    set "mRNA__multiQC.html", "*multiqc_data" optional true

  script:
    """

    multiqc -f .
    mv multiqc_report.html mRNA__multiQC.html

    """
}



// making the groups for differential gene expression

Kallisto_results_for_sleuth_1
  .map{ [ it[0].split('_')[0], it[1] ] }
  .groupTuple()
  .into{ Kallisto_results_for_sleuth_2; Kallisto_results_for_sleuth_3 }

Kallisto_results_for_sleuth_2
  .combine(Kallisto_results_for_sleuth_3)
  .map { it[0,2,1,3]}
  .join(comparisons_files_for_mRNA_Seq, by: [0,1])
  .dump(tag:'kalisto_sleuth') {"${it}"}
  .map{
    def list = []
    list.add( it[0,1].join('_vs_') )
    list.addAll( it[0..3] )
    return(list)
    }
  .set { Kallisto_results_for_sleuth_4 }






process DA_ATAC__doing_differential_abundance_analysis {
  tag "${COMP}"

  label "diffbind"
  // label "differential_abundance"
  // Error in loadNamespace(x) : there is no package called ‘edgeR’

  publishDir path: "${out_processed}/2_Differential_Abundance", 
    mode: "${pub_mode}", saveAs: {
     if (it.indexOf("__all_peaks.bed") > 0) "ATAC__all_peaks__bed/${it}"
     else if (it.indexOf("__diffbind_peaks_dbo.rds") > 0) 
                                            "ATAC__all_peaks__DiffBind/${it}"
     else if (it.indexOf("__diffbind_peaks_gr.rds") > 0) 
                                            "ATAC__all_peaks__gRange/${it}"
  }

  when: 
    do_atac

  input:
    set COMP, file(bed), file(bam) from Reads_and_peaks_for_diffbind_2

  output:
    set COMP, file('*__diffbind_peaks_gr.rds'), 
      file('*__diffbind_peaks_dbo.rds') \
      into Diffbind_peaks_for_annotating_them optional true
    set COMP, file('*__all_peaks.bed') \
      into All_detected_diffbind_peaks_for_background optional true

  shell:
    '''
    #!/usr/bin/env Rscript

    ##### loading data and libraries

    library(DiffBind)
    library(magrittr)
    library(parallel)

    source('!{projectDir}/bin/export_df_to_bed.R')
    cur_seqinfo = readRDS('!{params.cur_seqinfo}')

    COMP              = '!{COMP}'
    use_input_control = '!{params.use_input_control}'
    make_grey_list    = !{params.diffbind__make_grey_list}
    min_overlap       = !{params.diffbind__min_overlap}
    min_count         = !{params.diffbind__min_count}
    analysis_method   = !{params.diffbind__analysis_method}
    normalization     = !{params.diffbind__normalization}
    summits           = !{params.diffbind__summits}
    design            = !{params.diffbind__design}
    edger_tagwise     = !{params.diffbind__edger_tagwise}
    
    conditions = strsplit(COMP, '_vs_')[[1]]
    cond1 = conditions[1]
    cond2 = conditions[2]
    use_input_control %<>% toupper %>% as.logical

    if(make_grey_list & !use_input_control) stop('Grey lists cannot be created without an input control')
    

    ##### Preparing the metadata table
    bed_bam_files = list.files(pattern = '*.bam$|*_removal.bed$')
    cur_files = grep('diffbind_peaks', bed_bam_files, value = T, invert = T)
    df = data.frame(path = cur_files, stringsAsFactors=F)
    cursplit = sapply(cur_files, strsplit, '_')
    df$condition = sapply(cursplit, '[[', 1)
    df$replicate = sapply(cursplit, '[[', 2)
    df$id = paste0(df$condition, '_', df$replicate)
    df$id[df$condition == 'input'] = 'input'
    df$type = sapply(df$path, function(c1) {
                ifelse(
                  length(grep('reads', c1)) == 1, 'reads', ifelse(
                    length(grep('peaks', c1)) == 1, 'peaks', '')
                )
              })

    names_df1 = c('SampleID', 'Condition', 'Replicate', 'bamReads', 
                  'ControlID', 'bamControl', 'Peaks','PeakCaller')
    all_id = unique(df$id) %>% .[. != 'input']
    df1 = data.frame(matrix(nrow = length(all_id), 
                     ncol = length(names_df1)), stringsAsFactors=F)
    names(df1) = names_df1

    for(c1 in 1:length(all_id)){
      cur_id = all_id[c1]
      sel_reads = which(df$id == cur_id & df$type == 'reads')
      sel_peaks = which(df$id == cur_id & df$type == 'peaks')

      df1$SampleID[c1]   = cur_id
      df1$Condition[c1]  = df$condition[sel_reads]
      df1$Replicate[c1]  = df$replicate[sel_reads]
      df1$bamReads[c1]   = df$path[sel_reads]
      df1$Peaks[c1]      = df$path[sel_peaks]
      df1$PeakCaller[c1] = 'bed'

      if(use_input_control){
        sel_input_control_reads = 
            which(df$condition == 'input' & df$type == 'reads')
        df1$ControlID[c1] = df$id[sel_input_control_reads]
        df1$bamControl[c1] = df$path[sel_input_control_reads]
      } else {
        df1$ControlID <- NULL
        df1$bamControl <- NULL
      }
    }

    empty_peaks = which(sapply(df1$Peaks, file.size) == 0)
    if(length(empty_peaks) > 0) df1$Peaks[empty_peaks] = ''
    if(length(empty_peaks) == nrow(df1)) { quit('no') }



    ##### Running DiffBind

    dbo = tryCatch(
      expr = {

        config = data.frame(AnalysisMethod = analysis_method, RunParallel = F, 
                            singleEnd = F)

        dbo <- dba(sampleSheet = df1, minOverlap = min_overlap, config = config)

        dbo$config$edgeR$bTagwise = edger_tagwise

        if(make_grey_list) dbo <- dba.blacklist(dbo, blacklist = F, 
          greylist = cur_seqinfo)

        dbo <- dba.count(dbo, bParallel = F, bRemoveDuplicates = F, 
          fragmentSize = 1, minOverlap = min_overlap, 
          score = DBA_SCORE_NORMALIZED, bUseSummarizeOverlaps = F, 
          bSubControl = F, minCount = min_count, summits = summits)

        dbo <- dba.normalize(dbo, normalize = normalization, 
          library = DBA_LIBSIZE_BACKGROUND,  background = TRUE)

        dbo <- dba.contrast(dbo, categories = DBA_CONDITION, minMembers = 2, reorderMeta = list(Condition = cond2), design = design)

        dbo <- dba.analyze(dbo, bBlacklist = F, bGreylist = F,   bReduceObjects = F)

      }, error = function(e) {
        print(e$message)
        quit('no')
      }
    )

    saveRDS(dbo, paste0(COMP, '__diffbind_peaks_dbo.rds'))


    ##### Exporting all peaks as a data frame

    # extracting all peaks (note: th is the fdr threshold, so with th = 1 
    # we keep all peaks)
    all_peaks_gr = suppressWarnings(dba.report(dbo, th = 1))

    # recomputing the FDR to have more precise values (DiffBind round them 
    # at 3 digits)
    all_peaks_gr$FDR <- NULL
    all_peaks_gr$padj = p.adjust(data.frame(all_peaks_gr)$p.value, 
                                 method = 'BH')

    # adding the raw reads counts of each replicate
    cp_raw = dba.peakset(dbo, bRetrieve = TRUE, score = DBA_SCORE_READS)
    mcols(cp_raw) = apply(mcols(cp_raw), 2, round, 2)
    m <- findOverlaps(all_peaks_gr, cp_raw)
    names_subject = names(mcols(cp_raw))
    mcols(all_peaks_gr)[queryHits(m), names_subject] = 
        mcols(cp_raw)[subjectHits(m), names_subject]
    saveRDS(all_peaks_gr, paste0(COMP, '__diffbind_peaks_gr.rds'))


    ##### Exporting all peaks as a bed file

    all_peaks_df = as.data.frame(all_peaks_gr)
    all_peaks_df %<>% dplyr::mutate(
      name = rownames(.), 
      score = round(-log10(p.value), 2))
    all_peaks_df %<>% dplyr::rename(chr = seqnames)
    all_peaks_df %<>% dplyr::select(chr, start, end, name, score, strand)
    export_df_to_bed(all_peaks_df, paste0(COMP, '__all_peaks.bed'))

    '''
}


// not sure these lines of codes will be needed
// string_as_object <- function(x) eval(parse(text = x))
// min_overlap     = string_as_object('!params.dba_min_overlap')
// analysis_method = string_as_object('!params.dba_analysis_method')
// min_count       = string_as_object('!params.dba_min_count')
// score           = string_as_object('!params.dba_min_count')


// => generating the set of all peaks found in all replicates, and a set of 
// differentially abundant/accessible peaks (can also be called differentially 
// bound regions)


////// justification for the parameters used

//// dba: minOverlap = 1
// could be set to 1 (any peak), or 2 peaks present in at least 2 peak set. The latter is much more stringent in situations where one has only 2 replicates per condition, which is why I decided to go with 1 for now.
// https://support.bioconductor.org/p/57809/

//// grey list is applied only if we have an input control as the greylist is made based on the input control sample. It is applied at the very beginning to filter out the grey regions from the start.

//// dbo$config$AnalysisMethod = DBA_DESEQ2 
// https://support.bioconductor.org/p/98736/
// For normalization, edgeR uses the TMM method which relies on a core of sites that don't systematically change their binding affinities. While this assumption of usually true for RNA-seq, there are many ChIP-seq experiments where one sample group has a completely different binding profile than the other sample group. Often, in one condition there is very little binding, but in the other condition much additional binding has been induced. ...  This is the reason we changed the default analysis method in DiffBind from edgeR to DESeq2.

//// dba.count: summits = 75
// DiffBind Vignette
// When forming the global binding matrix consensus peaksets, DiﬀBind ﬁrst identiﬁes all unique
// peaks amongst the relevant peaksets. As part of this process, it merges overlapping peaks,
// replacing them with a single peak representing the narrowest region that covers all peaks
// that overlap by at least one base. There are at least two consequences of this that are worth
// noting.
// First, as more peaksets are included in analysis, the average peak width tends to become
// longer as more overlapping peaks are detected and the start/end points are adjusted outward
// to account for them. Secondly, peak counts may not appear to add up as you may expect
// due to merging. For example, if one peakset contains two small peaks near to each other,
// while a second peakset includes a single peak that overlaps both of these by at least one
// base, these will all be replaced in the merged matrix with a single peak. As more peaksets
// are added, multiple peaks from multiple peaksets may be merged together to form a single,
// wider peak. Use of the "summits" parameter is recommended to control for this widening
// eﬀect.
// ?dba.count
// summits: unless set to ‘FALSE’, summit heights (read pileup) and
// locations will be calculated for each peak.  
// If the value of ‘summits’ is ‘TRUE’ (or ‘0’), the summits
// will be calculated but the peaksets will be unaffected.  If
// the value is greater than zero, all consensus peaks will be
// re-centered around a consensus summit, with the value of
// ‘summits’ indicating how many base pairs to include upstream
// and downstream of the summit (so all consensus peaks will be
// of the same width, namely ‘2 * summits + 1’).
// https://support.bioconductor.org/p/9138959/
// "Generally, ATAC-seq benefits from using the summits parameter at a relatively low value (50 or 100) to focus on clear regions of open chromatin with less background signal. Note that often there are two "peaks" in ATAC fragment lengths."
// => by setting summits to 75, the summits of the consensus peaks will be extende by 75 bp on both side. 

//// dba.count: fragmentSize = 1, bUseSummarizeOverlaps = F
//// dbo$config$singleEnd = T
// ?dba.count
// bUseSummarizeOverlaps: logical indicating that ‘summarizeOverlaps’
//           should be used for counting instead of the built-in counting
//           code.  This option is slower but uses the more standard
//           counting function. If ‘TRUE’, all read files must be BAM
//           (.bam extension), with associated index files (.bam.bai
//           extension).  The ‘fragmentSize’ parameter must absent.
// fragmentSize: This value will be used as the length of the reads.  Each
//           read will be extended from its endpoint along the appropriate
//           strand by this many bases.  If set to zero, the read size
//           indicated in the BAM/BED file will be used. ‘fragmentSize’
//           may also be a vector of values, one for each ChIP sample plus
//           one for each unique Control library.
// https://support.bioconductor.org/p/9138959/
// The fragmentSize parameter is only used if you have single-end data (otherwise the paired-end insert size is used). In that case, it is a good idea to specify a value for this parameter to "extend" the single-end reads by the mean fragment size.
// bUseSummarizeOverlaps=TRUE is generally the preferred option. If it is set to FALSE, it will use some internal code that, among other things, does not match paired-end reads, resulting in counts that may be up to double (each end counted separately). That option also gives you more fine-grained control over how the counting is done via $config parameters.
// => we do want each pair to be analyzed separately. And we want to keep reads at 1 bp as this is the most precise signal we have (transposase-shifted 5' end of reads)



//// dbo <- dba.normalize(dbo, normalize = DBA_NORM_LIB, )
// https://support.bioconductor.org/p/p133577/

//// dba.normalize: 
//// dba.normalize: normalize = DBA_NORM_RLE, library = DBA_LIBSIZE_BACKGROUND,  background = TRUE
// => background
// help(dba.normalize)
// background: This parameter controls the option to use "background"
//           bins, which should not have differential enrichment between
//           samples, as the basis for normalizing (instead of using reads
//           counts overlapping consensus peaks).  When enabled, the
//           chromosomes for which there are peaks in the consensus
//           peakset are tiled into large bins and reads overlapping these
//           bins are counted.
//
// Vignette:
// An assumption in RNA-seq analysis, that the read count matrix reﬂects an unbiased repre-
// sentation of the experimental data, may be violated when using a narrow set of consensus
// peaks that are chosen speciﬁcally based on their rates of enrichment. It is not clear that
// using normalization methods developed for RNA-seq count matrices on the consensus reads
// will not alter biological signals; it is probably a good idea to avoid using the consensus count matrix (Reads in Peaks) for normalizing unless there is a good prior reason to expect balanced changes in binding.
// https://support.bioconductor.org/p/p133577/
// If you are getting different results, it is likely that there is a systematic shift in enriched regions in one of the conditions. Unless you have a specific reason to believe this is due to technical issues, this is likely a biological signal you want to retain. In this case, you should NOT use RLE on the binding matrix.
// DBA_NORM_LIB is a quick way to perform a minimal normalization. I would try:
// normalize <- dba.normalize(CTCF_narrowpeakscount, method = DBA_DESEQ2, 
//                            normalize = DBA_NORM_RLE,  background=TRUE)
// to run RLE against background bins rather than the binding matrix.
// This should probably be the default, but computing the background bins is somewhat compute intensive.
// => normalize
// Vignette:
// # About normalization
// DiﬀBind relies on three underlying core methods for normalization. These include the "na-
// tive" normalization methods supplied with DESeq2 and edgeR , as well as a simple library-
// based method. The normalization method is speciﬁed with the normalize parameter to
// dba.normalize
// The native DESeq2 normalization method is based on calculating the geometric mean for
// each gene across samples[6], and is referred to "RLE" or DBA_NORM_RLE in DiﬀBind .
// The native edgeR normalization method is based on the trimmed mean of M-values approach[7],
// and is referred to as "TMM" or DBA_NORM_TMM in DiﬀBind .
// A third method is also provided that avoids making any assumptions regarding the distribution
// of reads between binding sites in diﬀerent samples that may be speciﬁc to RNA-seq analysis
// and inappropriate for ChIP-seq analysis. This method ( "lib" or DBA_NORM_LIB ) is based on
// the diﬀerent library sizes for each sample, computing normalization factors to weight each
// sample so as to appear to have the same library size. For DESeq2 , this is accomplished
// by dividing the number of reads in each library by the mean library size. For edgeR , the
// normalization factors are all set to 1.0 , indicating that only library-size normalization should occur.
// Note that any of these normalization methods can be used with either the DESeq2 or
// edgeR analysis methods, as DiﬀBind converts normalization factors to work appropriately
// in either DESeq2 or edgeR .


//// Other references
// https://support.bioconductor.org/p/86594/
// https://support.bioconductor.org/p/107679/


// note: it may be better to use a fragment size of 150 together with the slopBed 75 bp shifted reads. Need to check that in more details some day. For now "fragmentSize = 1" should be good enough.




// cut -f1,1 ctl_2_peaks_kept_after_blacklist_removal_filtered.bed | sort | uniq // X


// I get this error with the fly test dataset:
// Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
//   every gene contains at least one zero, cannot compute log geometric means
// https://www.biostars.org/p/440379/

// df_tmp = data.frame(dba.peakset(dbo, bRetrieve = TRUE, score = DBA_SCORE_TMM_READS_FULL_CPM))
// data.frame(df_tmp)[, 6:ncol(df_tmp)]
// which(apply(data.frame(df_tmp)[, 6:ncol(df_tmp)], 1, function(x) all(x == 0)))

// here is the description of the changes: https://bioconductor.org/packages/release/bioc/news/DiffBind/NEWS

// # cur_greylist <- dba.blacklist(dbo, Retrieve=DBA_GREYLIST)


// ?dba.count
// bSubControl: logical indicating whether Control read counts are
//           subtracted for each site in each sample. If there are more
//           overlapping control reads than ChIP reads, the count will be
//           set to the ‘minCount’ value specified when ‘dba.count’ was
//           called, or zero if no value is specified.
// 
//           If ‘bSubControl’ is not explicitly specified, it will be set
//           to ‘TRUE’ unless a greylist has been applied (see
//           ‘dba.blacklist’).



// adding bUseSummarizeOverlaps = F to the dba.count call. This way we count as before, with 1 bp fragment length as it should be.
// https://rdrr.io/bioc/DiffBind/man/dba.count.html

// => removing the creation of a dbo1 object; the score DBA_SCORE_READS is obtained with the subsequent call to dba.peakset
// >         dbo1 <- dba.count(dbo, bParallel = F, bRemoveDuplicates = FALSE, fragmentSize = 1, minOverlap = 0, score = DBA_SCORE_READS)
// Warning message:
// In dba.count(dbo, bParallel = F, bRemoveDuplicates = FALSE, fragmentSize = 1,  :
//   No action taken, returning passed object...



// # note: I also disable bScaleControl since it is used only to normalize controls reads prior to the substraction with bSubControl. But since we don't do that anymore and use greylist it doesn't matter so we can remove this parameter
// https://support.bioconductor.org/p/69924/

// # About normalization
// DiﬀBind relies on three underlying core methods for normalization. These include the "na-
// tive" normalization methods supplied with DESeq2 and edgeR , as well as a simple library-
// based method. The normalization method is speciﬁed with the normalize parameter to
// dba.normalize
// The native DESeq2 normalization method is based on calculating the geometric mean for
// each gene across samples[6], and is referred to "RLE" or DBA_NORM_RLE in DiﬀBind .
// The native edgeR normalization method is based on the trimmed mean of M-values approach[7],
// and is referred to as "TMM" or DBA_NORM_TMM in DiﬀBind .
// A third method is also provided that avoids making any assumptions regarding the distribution
// of reads between binding sites in diﬀerent samples that may be speciﬁc to RNA-seq analysis
// and inappropriate for ChIP-seq analysis. This method ( "lib" or DBA_NORM_LIB ) is based on
// the diﬀerent library sizes for each sample, computing normalization factors to weight each
// sample so as to appear to have the same library size. For DESeq2 , this is accomplished
// by dividing the number of reads in each library by the mean library size. For edgeR , the
// normalization factors are all set to 1.0 , indicating that only library-size normalization should occur.
// Note that any of these normalization methods can be used with either the DESeq2 or
// edgeR analysis methods, as DiﬀBind converts normalization factors to work appropriately
// in either DESeq2 or edgeR .

// DBA_NORM_NATIVE: equal DBA_NORM_TMM for edgeR and DBA_NORM_RLE for DESeq2


// => it is in fact recommended to use greylist regions instead of substracting controls! We update the code accordingly
// https://support.bioconductor.org/p/97717/
// https://github.com/crazyhottommy/ChIP-seq-analysis/blob/master/part3_Differential_binding_by_DESeq2.md
// https://support.bioconductor.org/p/72098/#72173

// From the Diffbind manual it also says:
// Another way of thinking about greylists is that they are one way of using the information in
// the control tracks to improve the reliability of the analysis. Prior to version 3.0, the default
// in DiﬀBind has been to simply subtract control reads from ChIP reads in order to dampen
// the magnitude of enrichment in anomalous regions. Greylists represent a more principled way
// of accomplishing this. If a greylist has been applied, the current default in DiﬀBind is to not
// subtract control reads.


// we should keep the bFullLibrarySize parameter to TRUE
// The issue is how the data are normalized. Setting bFullLibrarySize=FALSE uses a standard RNA-seq normalization method (based on the number of reads in consensus peaks), which assumes that most of the "genes" don't change expression, and those that do are divided roughly evenly in both directions. Using the default bFullLibrarySize=TRUE avoids these assumptions and does a very simple normalization based on the total number of reads for each library.

// THe bTagwise argument should not be needed but we keep it just in case and to be clearer
// Next edgeR::estimateGLMTrendedDisp is called with the DGEList and a design matrix derived
// from the design formula. By default, edgeR::estimateGLMTagwiseDisp is called next; this
// can be bypassed by setting DBA$config$edgeR$bTagwise=FALSE . 15

// the diffbind package changed a lot since the version I was using. The contrast command should be changed, and the dba.analyze command too. I will need to read in more details the method
// dbo <- dba.contrast(dbo, ~Condition, minMembers = 2)
// Error in dba.analyze(dbo, bTagwise = FALSE, bFullLibrarySize = TRUE, bSubControl = TRUE,  :
//   unused arguments (bTagwise = FALSE, bFullLibrarySize = TRUE, bSubControl = TRUE)


// here is the description of the changes: https://bioconductor.org/packages/release/bioc/news/DiffBind/NEWS
  // Changes in version 3.0     

  // - Moved edgeR bTagwise parameter to $config option
  // - Remove normalization options bSubControl, bFullLibrarySize, filter, and filterFun from dba.analyze(), now set in dba.normalize().
  // - dba.normalize(): is a new function; - bSubControl default depend on presence of Greylist
  // - Default for bUseSummarizeOverlaps in dba.count is now TRUE
  // 
  // 	The previous methods for modelling are maintained for backward
  // 	compatibility, however they are not the default.  To repeat
  // 	earlier analyses, dba.contrast() must be called explicitly with
  // 	design=FALSE. See ?DiffBind3 for more information.
  // 
  //   The default mode for dba.count() is now to center around
  //   summits (resulting in 401bp intervals).  To to avoid
  //   recentering around summits as was the previous default, set
  //   summits=FALSE (or to a more appropriate value).
  // 
  //   Normalization options have been moved from dba.analyze() to the
  //   new interface function dba.normalize(). Any non-default
  //   normalization options must be specified using dba.normalize().

// from the manuel:
//   2. bFullLibrarySize: This is now part of the library parameter for dba.normalize . li
// brary=DBA_LIBSIZE_FULL is equivalent to bFullLibrarySize=TRUE , and library=DBA_LIBSIZE_PEAKREADS
// is equivalent to bFullLibrarySize=FALSE .

// from the R help:
// DBA_LIBSIZE_FULL: Full library size (all reads in library)
// DBA_LIBSIZE_PEAKREADS: Library size is Reads in Peaks


// rtracklayer::export(promoters, 'promoters.bed')

// note: diffbind_peaks: means the peaks from diffbind, not that these peaks are 
// diffbound (differentially bound). This set is in fact all the peaks that 
// diffbind found in all replicates. The corresponding bed file will be used 
// as a control for downstream enrichment tasks (CHIP, motifs, chromatin states).

//// Scores are not for anything besides downstream analysis (i.e. PCA, plots)
// https://support.bioconductor.org/p/69924/
// 2. Peak scores. Peak scores can be exported, but are mostly used for plotting 
// data that doesn't have a differential analysis run via dba.analyze(). For 
// example, after calling dba.count(), the heatmaps and PCA plots will use the 
// peak scores. You can see the peak scores using dba.peakset() with 
// bRetrieve=TRUE. The documentation for dba.count() describes the  different 
// scoring methods available (currently 16). These range from using the raw read 
// counts (with or without control reads subtracted), to variations of RPKM 
// normalized read counts and TMM normalized counts (the normalization method 
//   used by edgeR) and some scores based on summits (we do recommend using the 
//     summits option in dba.count). When you run an analysis using dba.analyze(), 
//     the normalization method used by the underlying differential expression 
//     package (edgeR or DESeq2) will be used for plots and reports (ie, when 
//       bCounts=TRUE in a call to dba.report).


process DA_ATAC__annotating_diffbind_peaks {
  tag "${COMP}"

  label "bioconductor"

  publishDir path: "${out_processed}/2_Differential_Abundance", 
    mode: "${pub_mode}", saveAs: {
     if (it.indexOf("_df.rds") > 0) "ATAC__all_peaks__dataframe/${it}"
     else if (it.indexOf("_cs.rds") > 0) "ATAC__all_peaks__ChIPseeker/${it}"
   }

  when: 
    do_atac

  input:
    set COMP, file(diffbind_peaks_gr), file(diffbind_peaks_dbo) \
      from Diffbind_peaks_for_annotating_them

  output:
    file('*.rds')
    set COMP, file("*_df.rds") into Annotated_diffbind_peaks_for_saving_tables
    set COMP, file("*_df.rds"), file(diffbind_peaks_dbo) \
      into Annotated_diffbind_peaks_for_plotting

  shell:
    '''
    #!/usr/bin/env Rscript

    ##### loading data, libraries and parameters
    library(GenomicFeatures)
    library(ChIPseeker)
    library(magrittr)

    COMP = '!{COMP}'
    tx_db <- loadDb('!{params.txdb}') 
    upstream = !{params.diffbind_peaks__promoter_up}
    downstream = !{params.diffbind_peaks__promoter_down}
    diffbind_peaks_gr = readRDS('!{diffbind_peaks_gr}')


    ##### annotating peaks
    anno_peak_cs = annotatePeak(diffbind_peaks_gr, 
                                TxDb = tx_db, 
                                tssRegion = c(-upstream, downstream), 
                                level = 'gene', 
                                overlap = 'all')

    ##### creating the data frame object
    anno_peak_gr = anno_peak_cs@anno
    df = as.data.frame(anno_peak_gr)
    anno_peak_df = cbind(peak_id = rownames(df), df, 
      anno_peak_cs@detailGenomicAnnotation)
    anno_peak_df$peak_id %<>% as.character %>% as.integer

    ##### saving results
    name0 = paste0(COMP, '__diffb_anno_peaks')
    saveRDS(anno_peak_df, paste0(name0, '_df.rds'))
    saveRDS(anno_peak_cs, paste0(name0, '_cs.rds'))

    '''
}

// cs: for ChIPseeker
// rtracklayer::export(anno_peak_df, paste0(name0, '.bed'))
// these are peaks of differential chromatin accessibility called by DiffBind

// overlap: one of 'TSS' or 'all', if overlap="all", then gene overlap with peak 
// will be reported as nearest gene, no matter the overlap is at TSS region or not.


// note that in df_annotated_peaks: the geneStart and geneEnd field reflect actually the transcript start and transcript end. This is if the level = 'transcript' option is used. If the option level = 'gene' is used then the coordinates correspond to gene. (same goes for distanceToTSS, genLength...)



process DA_ATAC__plotting_differential_abundance_results {
  tag "${COMP}"

  label "diffbind"
  // label "differential_abundance"

  publishDir path: "${out_fig_indiv}/${out_path}", mode: "${pub_mode}", 
    saveAs: {
      if (it.indexOf("_volcano.pdf") > 0) "ATAC__volcano/${it}"
      else if (it.indexOf("_PCA_1_2.pdf") > 0) "ATAC__PCA_1_2/${it}"
      else if (it.indexOf("_PCA_3_4.pdf") > 0) "ATAC__PCA_3_4/${it}"
      else if (it.indexOf("_other_plots.pdf") > 0) "ATAC__other_plots/${it}"
    }

  publishDir path: "${out_processed}/${out_path}/ATAC__non_annotated_peaks", 
             pattern: "*__ATAC_non_annotated_peaks.txt", mode: "${pub_mode}"

  input:
    val out_path from Channel.value('2_Differential_Abundance')
    set COMP, file(annotated_peaks), file(diffbind_object_rds) \
      from Annotated_diffbind_peaks_for_plotting

  output:
    set val("ATAC__volcano"), out_path, file('*__ATAC_volcano.pdf') \
      into ATAC_Volcano_for_merging_pdfs
    set val("ATAC__PCA_1_2"), out_path, file('*__ATAC_PCA_1_2.pdf') \
      into ATAC_PCA_1_2_for_merging_pdfs
    set val("ATAC__PCA_3_4"), out_path, file('*__ATAC_PCA_3_4.pdf') \
      into ATAC_PCA_3_4_for_merging_pdfs
    set val("ATAC__other_plots"), out_path, file('*__ATAC_other_plots.pdf') \
      into ATAC_Other_plot_for_merging_pdfs
    file("*__ATAC_non_annotated_peaks.txt")

  shell:
    '''
    #!/usr/bin/env Rscript


    ##### loading data and libraries

    library(ggplot2)
    library(magrittr)
    library(grid)
    library(DiffBind)

    source('!{projectDir}/bin/functions_plot_volcano_PCA.R')

    COMP = '!{COMP}'

    dbo = readRDS('!{diffbind_object_rds}')
    fdr_threshold = !{params.diffbind_plots__fdr_threshold}
    top_n_labels = !{params.diffbind_plots__top_n_labels}
    dbo$config$th = fdr_threshold
    df_genes_metadata = readRDS('!{params.df_genes_metadata}')
    df_annotated_peaks = readRDS('!{annotated_peaks}')


    ##### volcano plots

    res = df_annotated_peaks %>% dplyr::rename(L2FC = Fold, gene_id = geneId)
    df_genes_metadata_1 = df_genes_metadata[, c('gene_id', 'gene_name')]
    res %<>% dplyr::inner_join(df_genes_metadata_1, by = 'gene_id')

    pdf(paste0(COMP, '__ATAC_volcano.pdf'))
      plot_volcano_custom(res, sig_level = fdr_threshold, 
          label_column = 'gene_name', title = paste(COMP, 'ATAC'),
          top_n_labels = top_n_labels)
    dev.off()


    ##### PCA plots

    if(nrow(dbo$binding) == nrow(res)){
      
      # Handling unnanotated peaks
      v_gene_names = res$gene_name[match(1:nrow(dbo$binding), res$peak_id)]
      sel = which(is.na(v_gene_names))
      v_gene_names[is.na(v_gene_names)] = paste0('no_gene_', 1:length(sel))
      sink(paste0(COMP, '__ATAC_non_annotated_peaks.txt'))
      print(dbo$peaks[[1]][sel,])
      sink()
      
      # doing and plotting the PCA
      prcomp1 <- DiffBind__pv_pcmask__custom(dbo, nrow(dbo$binding), 
      cor = F, bLog = T)$pc
      rownames(prcomp1$x) = v_gene_names
      
      lp_1_2 = get_lp(prcomp1, 1, 2, paste(COMP, ' ', 'ATAC'))
      lp_3_4 = get_lp(prcomp1, 3, 4, paste(COMP, ' ', 'ATAC'))
      
      pdf(paste0(COMP, '__ATAC_PCA_1_2.pdf'))
      make_4_plots(lp_1_2)
      dev.off()
      
      pdf(paste0(COMP, '__ATAC_PCA_3_4.pdf'))
      make_4_plots(lp_3_4)
      dev.off()
      
    }


    ##### other plots

    pdf(paste0(COMP, '__ATAC_other_plots.pdf'))
        dba.plotMA(dbo, bNormalized = T)
        dba.plotHeatmap(dbo, main = 'all reads')
        first_2_replicates = sort(c(which(dbo$masks$Replicate.1), 
          which(dbo$masks$Replicate.2)))
        dba.plotVenn(dbo, mask = first_2_replicates, main = 'all reads')
    dev.off()

    '''
}

Merging_pdfs_channel = Merging_pdfs_channel
  .mix(ATAC_Volcano_for_merging_pdfs.groupTuple(by: [0, 1]))
Merging_pdfs_channel = Merging_pdfs_channel
  .mix(ATAC_PCA_1_2_for_merging_pdfs.groupTuple(by: [0, 1]))
Merging_pdfs_channel = Merging_pdfs_channel
  .mix(ATAC_PCA_3_4_for_merging_pdfs.groupTuple(by: [0, 1]))
Merging_pdfs_channel = Merging_pdfs_channel
  .mix(ATAC_Other_plot_for_merging_pdfs.groupTuple(by: [0, 1]))

// DiffBind__pv_pcmask__custom was adapted from DiffBind::pv.pcmask. 
// This little hacking is necessary because DiffBind functions plots PCA but 
// do not return the object.
// alternatively, the PCA could be made from scratch as done here: 
// https://support.bioconductor.org/p/76346/



process DA_ATAC__saving_detailed_results_tables {
  tag "${COMP}"

  label "r_basic"

  when: 
    do_atac

  input:
    set COMP, file(annotated_peaks) \
      from Annotated_diffbind_peaks_for_saving_tables

  output:
    set val('ATAC_detailed'), val('2_Differential_Abundance'), file('*.rds') \
      into ATAC_detailed_tables_for_formatting_table
    set COMP, file("*.rds") into ATAC_detailed_tables_for_splitting_in_subsets

  shell:
    '''
    #!/usr/bin/env Rscript

    ##### loading data and libraries
    library(magrittr)
    library(purrr)

    COMP = '!{COMP}'

    df_annotated_peaks = readRDS('!{annotated_peaks}')
    df_genes_metadata = readRDS('!{params.df_genes_metadata}')


    # setting up parameters
    conditions = tolower(strsplit(COMP, '_vs_')[[1]])
    cond1 = conditions[1]
    cond2 = conditions[2]

    # creating the res_detailed table
    res_detailed = df_annotated_peaks
    # adding gene metadata in a more readable format than what provided by 
    # default by ChIPseeker
    df_genes_metadata1 = dplyr::rename(df_genes_metadata, gene_chr = chr, 
      gene_start = start, gene_end = end, gene_width = width, 
      gene_strand = strand)
    res_detailed %<>% dplyr::select(-geneChr, -geneStart, -geneEnd, 
      -geneLength, -geneStrand)
    colnames(res_detailed) %<>% tolower
    res_detailed %<>% dplyr::rename(gene_id = geneid)
    res_detailed %<>% dplyr::inner_join(df_genes_metadata1, by = 'gene_id')

    # renaming columns
    res_detailed %<>% dplyr::rename(chr = seqnames, L2FC = fold, 
      pval = p.value, distance_to_tss = distancetotss, five_UTR = fiveutr, 
      three_UTR = threeutr)
    res_detailed$COMP = COMP

    # collapsing replicates values and renaming these columns as well
    colnames(res_detailed)[grep('conc_', colnames(res_detailed))] = 
      c('conc_cond1', 'conc_cond2')
    colns = colnames(res_detailed)
    colns1 = grep(paste0('^', cond1, '_'), colns)
    colns2 = grep(paste0('^', cond2, '_'), colns)
    res_detailed$counts_cond1 = 
      apply(res_detailed[, colns1], 1, paste, collapse = '|')
    res_detailed$counts_cond2 = 
      apply(res_detailed[, colns2], 1, paste, collapse = '|')
    res_detailed[, c(colns1, colns2)] <- NULL

    # reordering columns
    res_detailed$strand <- NULL  
      # strand information is not available for ATAC-Seq
    res_detailed %<>% dplyr::select(COMP, peak_id, chr, start, end, width, 
      gene_name, gene_id, pval, padj, L2FC, distance_to_tss, annotation, conc, 
      conc_cond1, conc_cond2, counts_cond1, counts_cond2, dplyr::everything())

    # adding the filtering columns
    res_detailed %<>% dplyr::mutate(
      FC_up = L2FC > 0,
      FC_down = L2FC < 0,

      PA_8kb = abs(distance_to_tss) < 8000,
      PA_3kb = abs(distance_to_tss) < 3000,
      PA_2u1d = distance_to_tss > -2000 & distance_to_tss < 1000,
      PA_TSS = distance_to_tss == 0,
      PA_genProm = genic | promoter,
      PA_genic = genic,
      PA_prom = promoter,
      PA_distNC = distal_intergenic | ( intron & !promoter & !five_UTR  & 
        !three_UTR  & !exon)
    )

    # saving table
    saveRDS(res_detailed, paste0(COMP, '__res_detailed_atac.rds'))

    '''
}


Formatting_csv_tables_channel = Formatting_csv_tables_channel
  .mix(ATAC_detailed_tables_for_formatting_table)



params.design__genes_to_remove_1 = workflow.launchDir.toString() + '/' +
                                   params.design__genes_to_remove

process DA_mRNA__doing_differential_abundance_analysis {
  tag "${COMP}"

  label "sleuth"

  publishDir path: "${out_processed}/2_Differential_Abundance", 
    mode: "${pub_mode}", saveAs: {
       if (it.indexOf("__mRNA_DEG_rsleuth.rds") > 0) 
          "mRNA__all_genes__rsleuth/${it}"
       else if (it.indexOf("__all_genes_prom.bed") > 0) 
          "mRNA__all_genes__bed_promoters/${it}"
  }

  when: 
    do_mRNA

  input:
    set COMP, cond1, cond2, file(kallisto_cond1), file(kallisto_cond2) \
      from Kallisto_results_for_sleuth_4

  output:
    set COMP, file('*__mRNA_DEG_rsleuth.rds') into Sleuth_results_for_plotting
    set COMP, file('*__mRNA_DEG_df.rds') into Sleuth_results_for_saving_tables
    set COMP, file('*__all_genes_prom.bed') \
      into All_detected_sleuth_promoters_for_background

  shell:
    '''
    #!/usr/bin/env Rscript

    library(sleuth)
    library(ggplot2)
    library(magrittr)

    df_genes_transcripts = readRDS('!{params.df_genes_transcripts}')
    df_genes_metadata    = readRDS('!{params.df_genes_metadata}')
    
    df_genes_to_remove = tryCatch(
      read.table('!{params.design__genes_to_remove_1}'), 
      error = function(e) NULL
    )

    COMP  = '!{COMP}'
    cond1 = '!{cond1}'
    cond2 = '!{cond2}'

    promoters_df = readRDS('!{params.promoters_df}')
    source('!{projectDir}/bin/export_df_to_bed.R')
    source('!{projectDir}/bin/get_prom_bed_df_table.R')

    if(is.null(df_genes_to_remove)){
      df_genes_transcripts_1 = df_genes_transcripts
    } else {
      df_genes_to_remove_1 = df_genes_to_remove %>% 
                             .[.$V1 %in% c(cond1, cond2), ]
      v_genes_to_remove = unique(df_genes_to_remove_1$V2)
      v_genes_id_to_remove = df_genes_metadata %>% 
                             .[.$gene_name %in% v_genes_to_remove,] %>% 
                             .$gene_id
      df_genes_transcripts_1 = df_genes_transcripts %>% 
                               .[!.$GENEID %in% v_genes_id_to_remove, ]
    }

    s2c = data.frame(path = dir(pattern = paste('kallisto', '*')), 
      stringsAsFactors = F)
    s2c$sample = sapply(s2c$path, function(x) strsplit(x, 'kallisto_')[[1]][2])
    s2c$condition = sapply(s2c$path, function(x) strsplit(x, '_')[[1]][2])
    levels(s2c$condition) = c(cond1, cond2)

		test_cond = paste0('condition', cond2)
		cond <- factor(s2c$condition)
		cond <- relevel(cond, ref = cond2)
		md <- model.matrix(~cond, s2c)
		colnames(md)[2] <- test_cond

    t2g <- dplyr::rename(df_genes_transcripts_1, target_id = TXNAME, 
        gene_id = GENEID)
    t2g$target_id %<>% paste0('transcript:', .)

    # Load the kallisto data, normalize counts and filter genes
	  sleo <- sleuth_prep(sample_to_covariates = s2c, full_model = md, 
      target_mapping = t2g, aggregation_column = 'gene_id', 
      transform_fun_counts = function(x) log2(x + 0.5), gene_mode = T)

    # Estimate parameters for the sleuth response error measurement (full) model
    sleo <- sleuth_fit(sleo)

    # Performing test and saving the sleuth object
    sleo <- sleuth_wt(sleo, test_cond)
    saveRDS(sleo, file = paste0(COMP, '__mRNA_DEG_rsleuth.rds'))

    # saving as dataframe and recomputing the FDR
    res <- sleuth_results(sleo, test_cond, test_type = 'wt', 
      pval_aggregate  = F)
    res %<>% .[!is.na(.$b), ]
    res$padj = p.adjust(res$pval, method = 'BH')
    res$qval <- NULL
    res %<>% dplyr::rename(gene_id = target_id)
    res1 = dplyr::inner_join(df_genes_metadata, res, by = 'gene_id')
    saveRDS(res1, file = paste0(COMP, '__mRNA_DEG_df.rds'))

    # exporting promoters of all detected genes
    promoters_df1 = promoters_df
    promoters_df1 %<>% .[.$gene_id %in% res1$gene_id, ]
    prom_bed_df = get_prom_bed_df_table(promoters_df1, res1)
    export_df_to_bed(prom_bed_df, paste0(COMP, '__all_genes_prom.bed'))

    '''
}

// => doing differential gene expression analysis


// with rsleuth v0.30
// length(which(is.na(res$b))) # 17437

// with rsleuth v0.29
// dim(res) # 2754   11 => NA values are automatically filtered out
// dim(df_genes_metadata) # 20191     8
// 20191 - 2754 # 17437

// sleuth_prep message
// 3252 targets passed the filter
// 2754 genes passed the filter

// the filtering function is the following:
// https://www.rdocumentation.org/packages/sleuth/versions/0.29.0/topics/basic_filter
// basic_filter(row, min_reads = 5, min_prop = 0.47)
// row        this is a vector of numerics that will be passedin
// min_reads  the minimum mean number of reads
// min_prop   the minimum proportion of reads to pass this filter

// note: we compute the FDR for both ATAC and mRNA seq to be sure to have consistent values generated by the same FDR method




process DA_mRNA__plotting_differential_abundance_results {
  tag "${COMP}"

  publishDir path: "${out_fig_indiv}/${out_path}", mode: "${pub_mode}", 
    saveAs: { 
      if (it.indexOf("__mRNA_volcano.pdf") > 0) "mRNA__volcano/${it}"
      else if (it.indexOf("__mRNA_PCA_1_2.pdf") > 0) "mRNA__PCA_1_2/${it}"
      else if (it.indexOf("__mRNA_PCA_3_4.pdf") > 0) "mRNA__PCA_3_4/${it}"
      else if (it.indexOf("__mRNA_other_plots.pdf") > 0) "mRNA__other_plots/${it}"
    }


  label "sleuth"

  input:
    val out_path from Channel.value('2_Differential_Abundance')
    set COMP, file(mRNA_DEG_rsleuth_rds) from Sleuth_results_for_plotting

  output:
    set val("mRNA__volcano"), out_path, file('*__mRNA_volcano.pdf') \
      into MRNA_Volcano_for_merging_pdfs
    set val("mRNA__PCA_1_2"), out_path, file('*__mRNA_PCA_1_2.pdf') \
      into MRNA_PCA_1_2_for_merging_pdfs
    set val("mRNA__PCA_3_4"), out_path, file('*__mRNA_PCA_3_4.pdf') \
      into MRNA_PCA_3_4_for_merging_pdfs
    set val("mRNA__other_plots"), out_path, file('*__mRNA_other_plots.pdf') \
      into MRNA_Other_plot_for_merging_pdfs

  shell:
    '''
    #!/usr/bin/env Rscript

    ##### Loading libraries and data

    library(sleuth)
    library(ggplot2)
    library(magrittr)
    library(grid)

    source('!{projectDir}/bin/functions_plot_volcano_PCA.R')

    sleo = readRDS('!{mRNA_DEG_rsleuth_rds}')
    df_genes_metadata = readRDS('!{params.df_genes_metadata}')

    COMP = '!{COMP}'
    conditions = strsplit(COMP, '_vs_')[[1]]
    cond1 = conditions[1]
    cond2 = conditions[2]
    test_cond = paste0('condition', cond2)

    fdr_threshold = !{params.sleuth_plots__fdr_threshold}
    top_n_labels  = !{params.sleuth_plots__top_n_labels}


    ##### volcano plots

    res_volcano <- sleuth_results(sleo, test_cond)
    res_volcano %<>% dplyr::rename(gene_id = target_id, L2FC = b)
    df_genes_metadata_1 = df_genes_metadata[, c('gene_id', 'gene_name')]
    res_volcano %<>% dplyr::inner_join(df_genes_metadata_1, by = 'gene_id')
    res_volcano %<>% dplyr::mutate(padj = p.adjust(pval, method = 'BH'))

    pdf(paste0(COMP, '__mRNA_volcano.pdf'))
      plot_volcano_custom(res_volcano, sig_level = fdr_threshold, 
        label_column = 'gene_name', title = paste(COMP, 'mRNA'),
        top_n_labels = top_n_labels
        )
    dev.off()


    ##### PCA plots

    # the pca is computed using the default parameters in the sleuth functions 
    # sleuth::plot_pca
    mat = sleuth:::spread_abundance_by(sleo$obs_norm_filt, 
      'scaled_reads_per_base')
    prcomp1 <- prcomp(mat)
    v_gene_id_name = df_genes_metadata_1$gene_name %>% 
      setNames(., df_genes_metadata_1$gene_id)
    rownames(prcomp1$x) %<>% v_gene_id_name[.]

    lp_1_2 = get_lp(prcomp1, 1, 2, paste(COMP, ' ', 'mRNA'))
    lp_3_4 = get_lp(prcomp1, 3, 4, paste(COMP, ' ', 'mRNA'))

    pdf(paste0(COMP, '__mRNA_PCA_1_2.pdf'))
      make_4_plots(lp_1_2)
    dev.off()

    pdf(paste0(COMP, '__mRNA_PCA_3_4.pdf'))
      make_4_plots(lp_3_4)
    dev.off()


    ##### other plots

    pdf(paste0(COMP, '__mRNA_other_plots.pdf'))
      plot_ma(sleo, test = test_cond, sig_level = fdr_threshold) + 
        ggtitle(paste('MA:', COMP))
      plot_group_density(sleo, use_filtered = TRUE, 
        units = "scaled_reads_per_base", trans = "log", 
        grouping = setdiff(colnames(sleo$sample_to_covariates), "sample"), 
        offset = 1) + ggtitle(paste('Estimated counts density:', COMP))
      # plot_scatter(sleo) + ggtitle(paste('Scatter:', COMP))
      # plot_fld(sleo, 1) + ggtitle(paste('Fragment Length Distribution:', 
      # COMP))
    dev.off()

    '''
}

Merging_pdfs_channel = Merging_pdfs_channel
  .mix(MRNA_Volcano_for_merging_pdfs.groupTuple(by: [0, 1]))
Merging_pdfs_channel = Merging_pdfs_channel
  .mix(MRNA_PCA_1_2_for_merging_pdfs.groupTuple(by: [0, 1]))
Merging_pdfs_channel = Merging_pdfs_channel
  .mix(MRNA_PCA_3_4_for_merging_pdfs.groupTuple(by: [0, 1]))
Merging_pdfs_channel = Merging_pdfs_channel
  .mix(MRNA_Other_plot_for_merging_pdfs.groupTuple(by: [0, 1]))




process DA_mRNA__saving_detailed_results_tables {
  tag "${COMP}"

  label "r_basic"

  publishDir path: "${out_tab_indiv}/2_Differential_Abundance/mRNA", 
    mode: pub_mode, enabled: params.tables__save_csv

  when: 
    do_mRNA

  input:
    set COMP, file(mRNA_DEG_df) from Sleuth_results_for_saving_tables

  output:
    set val('mRNA_detailed'), val('2_Differential_Abundance'), file('*.rds') \
      into MRNA_detailed_tables_for_formatting_table
    set COMP, file('*.rds') into MRNA_detailed_tables_for_splitting_in_subsets

  shell:
    '''
    #!/usr/bin/env Rscript

    ## Loading libraries and data

    library(magrittr)

    COMP = '!{COMP}'

    mRNA_DEG_df = readRDS('!{mRNA_DEG_df}')


    ## Saving a detailed results table

    res_detailed = mRNA_DEG_df
    res_detailed %<>% dplyr::rename(L2FC = b)
    res_detailed$COMP = COMP
    res_detailed %<>% dplyr::select(COMP, chr, start, end, width, strand, 
      gene_name, gene_id, entrez_id, pval, padj, L2FC, dplyr::everything())
    saveRDS(res_detailed, paste0(COMP, '__res_detailed_mRNA.rds'))

    '''
}


Formatting_csv_tables_channel = Formatting_csv_tables_channel
  .mix(MRNA_detailed_tables_for_formatting_table)






if(!do_atac & do_mRNA){
  MRNA_detailed_tables_for_splitting_in_subsets
  .set{ Differational_abundance_results_for_splitting_in_subsets }
}

if(do_atac & !do_mRNA){
  ATAC_detailed_tables_for_splitting_in_subsets
  .set{ Differational_abundance_results_for_splitting_in_subsets }
}

if(do_atac & do_mRNA){
  ATAC_detailed_tables_for_splitting_in_subsets
  .mix(MRNA_detailed_tables_for_splitting_in_subsets)
  // format: COMP, res_detailed
  .groupTuple()
  // format: COMP, [ res_detailed_atac, res_detailed_mRNA ]
  .filter{ comp, res_files -> 
    ( do_atac && !do_mRNA  && res_files.size() == 1) ||
    (!do_atac &&  do_mRNA  && res_files.size() == 1) ||
    ( do_atac &&  do_mRNA  && res_files.size() == 2)
  }
  .dump(tag: 'both_data')
  .set{ Differational_abundance_results_for_splitting_in_subsets }
}

do_mRNA_lgl = do_mRNA.toString().toUpperCase()
do_atac_lgl = do_atac.toString().toUpperCase()


process DA_split__splitting_differential_abundance_results_in_subsets {
  tag "${COMP}"

  label "r_basic"

  publishDir path: "${out_processed}/2_Differential_Abundance", 
  mode: "${pub_mode}", saveAs: {
    if (it.indexOf("__genes.rds") > 0) "DA_split__genes_rds/${it}"
    else if   (it.indexOf(".bed") > 0) "DA_split__bed_regions/${it}"
  }

  input:
    set COMP, file(res_detailed) \
      from Differational_abundance_results_for_splitting_in_subsets

  output:
    set val('res_simple'), val('2_Differential_Abundance'), 
      file("*__res_simple.rds") \
      into Res_simple_table_for_formatting_table optional true
    set val('res_filter'), val('2_Differential_Abundance'), 
      file("*__res_filter.rds") \
      into Res_filter_table_for_formatting_table optional true
    set COMP, file("*__genes.rds") \
      into DA_genes_split_for_doing_enrichment_analysis, 
           DA_genes_for_plotting_venn_diagrams optional true
    set COMP, file("*__regions.bed") \
      into DA_regions_split_for_doing_enrichment_analysis optional true

  shell:
    '''
    #!/usr/bin/env Rscript


    ################################
    ## loading data and libraries
    library(magrittr)
    library(purrr)
    library(data.table)

    source('!{projectDir}/bin/splitting_DAR_in_subsets_functions.R')
    source('!{projectDir}/bin/export_df_to_bed.R')
    source('!{projectDir}/bin/get_prom_bed_df_table.R')
    source('!{projectDir}/bin/read_from_nextflow.R')
    source('!{projectDir}/bin/grepl_filter.R')
    
    COMP = '!{COMP}'

    do_mRNA = !{do_mRNA_lgl}
    do_atac = !{do_atac_lgl}

    lf = list.files()

    promoters_df = readRDS('!{params.promoters_df}')

    TT       = '!{params.split__threshold_type}'
    TV_split = read_from_nextflow(
      '!{params.split__threshold_values}') %>% as.numeric
    FC_split = read_from_nextflow(
      '!{params.split__fold_changes}')
    PA_split = read_from_nextflow(
      '!{params.split__peak_assignment}')


    ################################
    ## creating the res_simple table

    # combining atac and mRNA results
    if(do_atac) {
      res_detailed_atac = readRDS(grep('atac.rds', lf, value = T))

      # adding the aggregated PA filter column
      res_detailed_atac$PA_all = T
      PA_columns_all = grep('PA_', colnames(res_detailed_atac), value = T)
      res_detailed_atac$PA = 
        get_merged_columns(res_detailed_atac, paste0('PA_', PA_split), 'PA')

      res_simple_atac = res_detailed_atac %>% 
        dplyr::mutate(transcript_id = NA, ET = 'ATAC') %>% 
        dplyr::select(COMP, peak_id, chr, gene_name, gene_id, pval, padj, 
          L2FC, PA, ET)

      # adding the rank column
      res_simple_atac %<>% dplyr::arrange(padj, desc(abs(L2FC)))
      dt = data.table(res_simple_atac)
      dt[, FC := ifelse(L2FC > 0, 'up', 'down')]
      dt[L2FC == 0, FC := 'NA']
      dt$rank = 0
      dt[(!duplicated(dt[, c('gene_id', 'FC')])), rank := 1:.N]
      dt[, rank := rank[1], .(cumsum(rank != 0))]
      dt[, FC := NULL]
      res_simple_atac = dt %>% copy %>% setDF

      res_simple = res_simple_atac
    }

    if(do_mRNA) {
      res_detailed_mRNA = readRDS(grep('mRNA.rds', lf, value = T))
      res_simple_mRNA = res_detailed_mRNA %>% dplyr::select(COMP, chr, 
        gene_name, gene_id, pval, padj, L2FC)
      res_simple_mRNA = cbind(peak_id = 'Null', res_simple_mRNA, ET = 'mRNA', 
        PA = 'Null', stringsAsFactors = F)
      res_simple_mRNA %<>% dplyr::select(COMP, peak_id, chr, gene_name, 
        gene_id, pval, padj, L2FC, PA, ET)

      # adding the rank column
      res_simple_mRNA %<>% dplyr::arrange(padj, desc(abs(L2FC)))
      res_simple_mRNA$rank = 1:nrow(res_simple_mRNA)

      res_simple = res_simple_mRNA
    }

    if(do_mRNA & do_atac) res_simple = rbind(res_simple_atac, res_simple_mRNA)

    # adding the aggregated FC filter column
    res_simple %<>% dplyr::mutate(FC_all = T, FC_up = L2FC > 0, 
      FC_down = L2FC < 0)
    res_simple$FC = get_merged_columns(res_simple, paste0('FC_', FC_split), 
      'FC')

    # adding the aggregated TV filter column
    for(TV in TV_split) res_simple[[paste0('TV_', TV)]] = 
      filter_entries_by_threshold(res_simple, TT, TV)
    res_simple$TV = 
      get_merged_columns(res_simple, paste0('TV_', TV_split), 'TV')
    res_simple$TV %<>% gsub('^$', 'NS', .)  # NS: None Significant
    
    # filtering and reordering columns
    res_simple %<>% dplyr::select(ET, PA, FC, TV, COMP, peak_id, chr, 
      gene_name, gene_id, pval, padj, L2FC)

    saveRDS(res_simple, paste0(COMP, '__res_simple.rds'))


    ################################
    ## creating the res_filter table

    df_split = expand.grid(TV = TV_split, FC = FC_split, PA = PA_split, 
      stringsAsFactors = F)
    lres_filter = list()

    for(c1 in 1:nrow(df_split)){
      TV1 = df_split$TV[c1]
      FC1 = df_split$FC[c1]
      PA1 = df_split$PA[c1]
            
      res_filter_atac = res_simple %>% 
        dplyr::filter(grepl_filter(TV1, TV) & 
                      grepl_filter(FC1, FC) & grepl_filter(PA1, PA) & ET == 'ATAC')
      res_filter_mRNA = res_simple %>% 
        dplyr::filter(grepl_filter(TV1, TV) & grepl_filter(FC1, FC) & ET == 'mRNA')

      # adding the Experiment Type "both"
      ATAC_genes = res_filter_atac$gene_id
      mRNA_genes = res_filter_mRNA$gene_id
      both_genes = res_filter_atac$gene_id %>% .[. %in% mRNA_genes]

      res_filter_atac_both = res_filter_atac %>% 
        dplyr::filter(gene_id %in% both_genes) %>% 
        dplyr::mutate(ET = 'both_ATAC')
      res_filter_mRNA_both = res_filter_mRNA %>% 
        dplyr::filter(gene_id %in% both_genes) %>% 
        dplyr::mutate(ET = 'both_mRNA')

      res_filter_tmp = rbind(res_filter_atac_both, res_filter_mRNA_both, 
        res_filter_atac, res_filter_mRNA)
      res_filter_tmp %<>% dplyr::mutate(TV = TV1, FC = FC1, PA = PA1)
      res_filter_tmp$PA[res_filter_tmp$ET == 'mRNA'] = 'Null'

      lres_filter[[c1]] = res_filter_tmp

    }

    res_filter = do.call(rbind, lres_filter)
    res_filter %<>% .[!duplicated(.), ]

    saveRDS(res_filter, paste0(COMP, '__res_filter.rds'))


    ################################
    ## splitting results in subsets and exporting as bed and gene lists

    if(do_atac) {
      atac_bed_df = res_detailed_atac
      atac_bed_df %<>% dplyr::mutate(score = round(-log10(padj), 2), 
                    name = paste(gene_name, peak_id, sep = '_'), strand = '*')
      atac_bed_df %<>% dplyr::select(chr, start, end, name, score, strand, 
        gene_id)
    }

    if(do_mRNA) {
      promoters_df1 = promoters_df
      prom_bed_df = get_prom_bed_df_table(promoters_df1, res_simple_mRNA)
    }

    res_filter$gene_peak = 
      apply(res_filter[, c('gene_name', 'peak_id')], 1, paste, collapse = '_')

    df_split1 = res_filter %>% dplyr::select(ET:TV) %>% .[!duplicated(.),]

    for(c1 in 1:nrow(df_split1)){
      ET1  = df_split1$ET[c1]
      PA1  = df_split1$PA[c1]
      FC1  = df_split1$FC[c1]
      TV1 = df_split1$TV[c1]

      df = subset(res_filter, ET == ET1 & PA == PA1 & FC == FC1 & TV == TV1)
      DA_genes = unique(df$gene_id)
      NDA_genes = subset(res_simple, ET == ET1 & !gene_id %in% DA_genes, 
        'gene_id')$gene_id
      lgenes = list(DA = DA_genes, NDA = NDA_genes)

      key = paste(ET1, PA1, FC1, TV1, COMP, sep = '__')

      # exporting bed files
      if(do_mRNA && ET1 %in% c('mRNA', 'both_mRNA')) {
        cur_bed = prom_bed_df %>% .[.$gene_id %in% DA_genes, ]
      }

      if(do_atac && ET1 %in% c('ATAC', 'both_ATAC')) {
        cur_bed = atac_bed_df %>% .[.$name %in% df$gene_peak, ]
      }

      bed_name = paste0(key, '__regions.bed')
      export_df_to_bed(cur_bed, bed_name)

      # exporting gene list
      key %<>% gsub('both_....', 'both', .) # renaming both_{ATAC,mRNA} as both
      saveRDS(lgenes, paste0(key, '__genes.rds'))

    }

    '''
}

// bed_name = paste0(key, '__diff_expr_genes_prom.bed')
// bed_name = paste0(key, '__diff_ab_peaks.bed')

// note that the rank column for atac increase for each peak that is assigned to 
// a new gene. This way as many genes are assigned in ATAC-Seq and mRNA-Seq

// commands to check venn diagrams sizes:
// dt[rank < 1001][L2FC > 0]$gene_id %>% unique %>% length
// dt[rank < 1001][L2FC < 0]$gene_id %>% unique %>% length


// a=      res_simple %>% .[.$ET == 'ATAC', ] %>% 
//     .[.$rank < 1001 & .$FC > 0, ] %>% .$gene_id %>% unique
// c = a %>% .[!. %in% b]
// dt[gene_id == c[1]]


// res_simple[, c('ET', 'FC', 'TV_1000')] %>% table
// to run after this line:    for(TV in TV_split) res_simple[[paste0('TV_', ...

// res_simple_atac %>% .[.$rank < 1000] %>%  .[.$L2FC < 0 ,]  %>% 
// .[!duplicated(.$gene_id),] %>% nrow
// res_simple_atac %>%  .[.$L2FC > 0 ,]  %>% .[!duplicated(.$gene_id),] %>% nrow


Formatting_csv_tables_channel = Formatting_csv_tables_channel
  .mix(Res_simple_table_for_formatting_table)
Formatting_csv_tables_channel = Formatting_csv_tables_channel
  .mix(Res_filter_table_for_formatting_table)




// Note: in the res_filter table, each entry is multiplied at maximum by the 
// number of rows in the df_split table times 2 (for the ET column)

// why we remove 1bp at the start for bed export from a grange object
// https://www.biostars.org/p/89341/
// The only trick is remembering the BED uses 0-based coordinates.
// note that the df that are input to the function export_df_to_bed all come 
// from a grange object (output of )






//// Adding keys to genes sets

// note: there are slightly more items in the peaks channels, since the both 
// sets are duplicated (both_ATAC and both_mRNA)


DA_genes_split_for_doing_enrichment_analysis
  // COMP, [ multiple_rds_files ] 
  //  (files format: path/ET__PA__FC__TV__COMP__genes.rds)
  .map{ it[1] }.flatten().toList()
  // [ all_rds_files ]
  .into{ 
    DA_genes_list_for_computing_genes_self_overlaps ;
    DA_genes_for_computing_genes_self_overlaps_1
  }

DA_genes_for_computing_genes_self_overlaps_1
  .flatten()
  // one rds_file per line
  .map{ [ it.name.replaceFirst(~/__genes.*/, ''), it ] }
  // key (ET__PA__FC__TV__COMP), rds_file
  .into{ 
    DA_genes_for_computing_genes_self_overlaps_2 ; 
    DA_genes_for_computing_functional_annotations_overlaps_1
  }


DA_genes_for_computing_genes_self_overlaps_2
  // format: key (ET__PA__FC__TV__COMP), rds_file (lgenes = list(DA, NDA))
  .combine(DA_genes_list_for_computing_genes_self_overlaps)
  // format: key, rds_file, rds_files
  .map{ [ "${it[0]}__genes_self", 'genes_self', it[1], it[2..-1] ] }
  // format: key (ET__PA__FC__TV__COMP__DT), DT, rds_file, [ rds_files ]
  .dump(tag: 'genes_self')
  .set{ DA_genes_for_computing_genes_self_overlaps_3 }



//// Adding backgrounds to peak sets

// A background is needed for downstream analysis to compare DEG or DBP to all 
// genes or all peaks. This background is all_peaks found by DiffBind for the 
// ATAC and both_ATAC entries. And it is the promoters of genes detected by 
// sleuth for mRNA and both_mRNA entries.

DA_regions_split_for_doing_enrichment_analysis
  // COMP, multiple_bed_files 
  //   (files format: path/ET__PA__FC__TV__COMP__peaks.bed)
  .map{ it[1] }.flatten()
  // bed_file
  .map{ [ it.name.replaceFirst(~/__regions.bed/, ''), it ] }
  // key (ET__PA__FC__TV__COMP), bed_file
  .dump(tag:'DA_regions')
  .into{ DA_regions_split_ATAC ; DA_regions_split_mRNA }

DA_regions_split_ATAC
  // key, DA_regions
  .filter{ it[0] =~ /.*ATAC.*/ }
  .combine(All_detected_diffbind_peaks_for_background)
  // key, DA_regions, COMP, all_peaks
  .dump(tag:'DA_regions_ATAC')
  .set{ DA_regions_split_ATAC1 }

DA_regions_split_mRNA
  .filter{ it[0] =~ /.*mRNA.*/ }
  .combine(All_detected_sleuth_promoters_for_background)
  // key, diff_expr_prom_bed, COMP, all_genes_prom_bed
  .dump(tag:'DA_regions_mRNA')
  .set{ DA_regions_split_mRNA1 }

DA_regions_split_ATAC1
  .mix(DA_regions_split_mRNA1)
  // format: key, DA_regions, all_regions
  .filter{ it[0].split('__')[4] ==~ it[2] }
  // keeping only entries for which the COMP match with the key
  .map{ it[0, 1, 3] }
  // key, DA_regions, all_regions
  .filter{ it[1].readLines().size() >= params.min_entries_DA_bed }
  // keeping only entries with more than N dif bound/expr regions
  .dump(tag:'DA_regions_with_bg')
  .into{ 
    DA_regions_with_bg_for_computing_peaks_overlaps_1
    DA_regions_with_bg_for_computing_motifs_overlaps_1
  }



process DA_split__plotting_venn_diagrams {
  tag "${COMP}"

  label "venndiagram"

  publishDir path: "${out_fig_indiv}/2_Differential_Abundance", 
    mode: "${pub_mode}", saveAs: {
      if (it.indexOf("__venn_up_or_down.pdf") > 0) 
        "Venn_diagrams__two_ways/${it}"
      else if (it.indexOf("__venn_up_and_down.pdf") > 0) 
        "Venn_diagrams__four_ways/${it}"
  }

  input:
    set COMP, file('*') from DA_genes_for_plotting_venn_diagrams

  output:
    set val("Venn_diagrams__two_ways"), val("2_Differential_Abundance"), 
      file('*__venn_up_or_down.pdf') \
      into Venn_up_or_down_for_merging_pdfs optional true
    set val("Venn_diagrams__four_ways"), val("2_Differential_Abundance"), 
      file('*__venn_up_and_down.pdf') \
      into Venn_up_and_down_for_merging_pdfs optional true

  shell:
    '''
    #!/usr/bin/env Rscript

    library(VennDiagram)

    source('!{projectDir}/bin/plot_venn_diagrams.R')


    all_files = list.files(pattern = '*.rds')

    df = data.frame(do.call(rbind, lapply(all_files, 
      function(x) strsplit(x, '__')[[1]][-6])), stringsAsFactors = F)
    colnames(df) = c('ET', 'PA', 'FC', 'TV', 'COMP')

    PAs = unique(df$PA)
    PAs = PAs[PAs != 'Null']
    TVs = unique(df$TV)
    COMP = df$COMP[1]

    get_file_name <- function(ET, PA, FC, TV){
      paste(ET, PA, FC, TV, COMP, 'genes.rds', sep = '__')
    }

    for(PA1 in PAs){
      for(TV1 in TVs){

        atac_up = get_file_name('ATAC', PA1, 'up', TV1)
        atac_down = get_file_name('ATAC', PA1, 'down', TV1)
        mrna_up = get_file_name('mRNA', 'Null', 'up', TV1)
        mrna_down = get_file_name('mRNA', 'Null', 'down', TV1)

        a_u = file.exists(atac_up)
        a_d = file.exists(atac_down)
        m_u = file.exists(mrna_up)
        m_d = file.exists(mrna_down)

        if(a_u & m_u) {
          lgenes = list(atac_up = readRDS(atac_up)$DA, mrna_up = 
            readRDS(mrna_up)$DA)
          prefix = paste(PA1, 'up', TV1, COMP, sep = '__')
          plot_venn_diagrams(lgenes, prefix)
        }

        if(a_d & m_d) {
          lgenes = list(atac_down = readRDS(atac_down)$DA, mrna_down = 
            readRDS(mrna_down)$DA)
          prefix = paste(PA1, 'down', TV1, COMP, sep = '__')
          plot_venn_diagrams(lgenes, prefix)
        }

        if(a_u & a_d & m_u & m_d) {
          lgenes = list(atac_up = readRDS(atac_up)$DA, mrna_up = 
            readRDS(mrna_up)$DA, atac_down = readRDS(atac_down)$DA, 
            mrna_down = readRDS(mrna_down)$DA)
          prefix = paste(PA1, TV1, COMP, sep = '__')
          plot_venn_diagrams(lgenes, prefix)
        }

      }
    }

  '''
}

Merging_pdfs_channel = Merging_pdfs_channel.mix(Venn_up_or_down_for_merging_pdfs
  .groupTuple(by: [0, 1])
  .map{ it.flatten() }
  .map{ [ it[0], it[1], it[2..-1] ] })

Merging_pdfs_channel = Merging_pdfs_channel.mix(Venn_up_and_down_for_merging_pdfs
  .groupTuple(by: [0, 1])
  .map{ it.flatten() }
  .map{ [ it[0], it[1], it[2..-1] ] })




//// importing groups of comparisons to plot together on the overlap matrix and 
// the heatmaps

comparisons_grouped = file(params.design__groups)
Channel
  .from( comparisons_grouped.readLines() )
  .map { it.split() }
  .map { [ it[0], it[1..-1].join('|') ] }
  // format: GRP, [ comp_order ]
  .dump(tag:'comp_group') {"comparisons grouped: ${it}"}
  // .into { comparisons_grouped_for_overlap_matrix ; 
  // comparisons_grouped_for_heatmap }
  .set { comparisons_grouped_for_heatmap }





DA_genes_for_computing_functional_annotations_overlaps_1
  // key (ET__PA__FC__TV__COMP), rds_file
  .combine(params.func_anno_databases)
  // key (ET__PA__FC__TV__COMP), rds_file, func_anno
  .map{ key, rds_file, func_anno -> 
    data_type = "func_anno_" + func_anno
    new_key = key + "__" + data_type
    [ new_key, data_type, func_anno, rds_file ]
  }
  // key (ET__PA__FC__TV__COMP__DT), data_type, func_anno, rds_file 
  .dump(tag: 'func_anno')
  .set{ DA_genes_for_computing_functional_annotations_overlaps_2 }



process Enrichment__computing_functional_annotations_overlaps {
  tag "${key}"

  label "bioconductor"

  when: 
    params.do_func_anno_enrichment && params.do_any_enrichment

  input:
    set key, data_type, func_anno, file(gene_set_rds) \
      from DA_genes_for_computing_functional_annotations_overlaps_2

  output:
    set key, data_type, file('*__counts.csv') \
      into Functional_annotations_overlaps_for_computing_pvalues optional true

  shell:
    '''
    #!/usr/bin/env Rscript

    library(clusterProfiler)
    library(magrittr)

    key = '!{key}'
    lgenes = readRDS('!{gene_set_rds}')
    func_anno = '!{func_anno}'
    org_db = AnnotationDbi::loadDb('!{params.org_db}')
    df_genes_metadata = readRDS('!{params.df_genes_metadata}')
    kegg_environment = readRDS('!{params.kegg_environment}')
    use_nda_as_bg_for_func_anno = !{params.use_nda_as_bg_for_func_anno}
    simplify_cutoff = !{params.simplify_cutoff}



    gene_set_enrich <- function(lgenes, type){
      DA_genes = lgenes$DA
      universe = NULL
      if(use_nda_as_bg_for_func_anno) universe = lgenes$NDA

      if(type == 'KEGG') {

        vec = df_genes_metadata$entrez_id %>% 
          setNames(., df_genes_metadata$gene_id)
        gene_set_entrez = vec[DA_genes]

        res <- clusterProfiler:::enricher_internal(gene_set_entrez, 
          pvalueCutoff  = 1, qvalueCutoff  = 1, pAdjustMethod = 'BH', 
          universe = universe, USER_DATA = kegg_environment)

        if(is.null(res)) error('NULL output')
        res
      }
        else {
          simplify( 
            enrichGO( 
              gene = DA_genes, OrgDb = org_db, keyType = 'ENSEMBL', 
              ont = type, pAdjustMethod = 'BH', pvalueCutoff = 1, 
              qvalueCutoff  = 1, universe = universe), 
            cutoff = simplify_cutoff, by = 'p.adjust', select_fun = min)
        }
    }

    robust_gene_set_enrich <- function(lgenes, type){
      tryCatch(gene_set_enrich(lgenes, type), error = 
        function(e) {print(paste('No enrichment found'))})
    }

    res = robust_gene_set_enrich(lgenes, func_anno)

    if(class(res) == 'enrichResult') {
      df = res@result

      if(nrow(df) != 0) {

          # converting IDs from entrez to Ensembl
          if(func_anno == 'KEGG'){
            vec = df_genes_metadata$gene_id %>% 
              setNames(., df_genes_metadata$entrez_id)
            ensembl_ids = purrr::map_chr(strsplit(df$geneID, '/'), 
                            ~paste(unname(vec[.x]), collapse = '/'))
            df$geneID = ensembl_ids
          }

          # extracting count columns, reformatting and saving results
          df %<>% tidyr::separate(GeneRatio, c('ov_da', 'tot_da'), sep = '/')
          df %<>% tidyr::separate(BgRatio, c('ov_nda', 'tot_nda'), sep = '/')
          df[, c('ov_da', 'tot_da', 'ov_nda', 'tot_nda')] %<>% 
            apply(2, as.integer)

          df %<>% dplyr::rename(tgt_id = ID, tgt = Description, 
            genes_id = geneID)
          df %<>% dplyr::select(tgt, tot_da, ov_da, tot_nda, ov_nda, tgt_id, 
            genes_id)

          write.csv(df, paste0(key, '__counts.csv'), row.names = F)

      }
    }

    '''
}

Overlap_tables_channel = Overlap_tables_channel
  .mix(Functional_annotations_overlaps_for_computing_pvalues)

// universe = ifelse(use_nda_as_bg_for_func_anno, lgenes$NDA, NULL) # => this fails: "error replacement has length zero"

// https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
// hacking the enrichKEGG function to get it to work offline and inside the containers; otherwise I get the error 
// In utils::download.file(url, quiet = TRUE, method = method, ...) :
//   URL 'https://rest.kegg.jp/link/cel/pathway': status was 'Failure when receiving data from the peer'
// source code: https://rdrr.io/bioc/clusterProfiler/src/R/enrichKEGG.R
// KEGG_DATA <- prepare_KEGG(species, "KEGG", keyType)
// search_kegg_organism('cel', by='kegg_code')
// my_kegg_data = clusterProfiler:::prepare_KEGG('cel', "KEGG")
// 
// gene_set_entrez = c("176003", "178531", "175223", "177619", "180999", "3565479",
// "174059", "180862", "189649", "173657", "185570", "181697", "188978",
// "180929", "36805029", "186080", "179225", "179595", "179712",
// "179425", "181325", "171788", "178203", "173301", "175423", "178130",
// "183903", "179990", "187830", "176278", "176155", "3565939",
// "177051", "180901", "174198", "177074", "178191", "188458", "175469",
// "182932")
// organism_key = 'cel'
// universe = NULL
// res = enrichKEGG( gene = gene_set_entrez, pvalueCutoff = 1, qvalueCutoff  = 1, pAdjustMethod = 'BH', organism = organism_key, keyType = 'ncbi-geneid', universe = universe)
// res = enrichKEGG( gene = gene_set_entrez, pvalueCutoff = 1, qvalueCutoff  = 1, pAdjustMethod = 'BH', organism = organism_key, keyType = 'ncbi-geneid', universe = universe, use_internal_data = T)
// 

// For enrichGO, the go terms are stored into a package and accessed liked that:  goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
// so everything is run in local within the clusterProfiler:::get_GO_data function

// for enrichKEGG, the KEGG database is fetched online. We use a preparsed one here instead


// res <- clusterProfiler:::enricher_internal(gene_set_entrez,
//                          pvalueCutoff  = 1,
//                          pAdjustMethod = 'BH',
//                          universe      = universe,
//                          qvalueCutoff  = 1,
//                          USER_DATA = my_kegg_data1)
// 
//  R.utils::setOption("clusterProfiler.download.method",'auto')
// 
// specie_long = 'caenorhabditis_elegans'
// 
// 
// 
// 
// options("clusterProfiler.download.method")
// 
// my_kegg_data1 = clusterProfiler:::prepare_KEGG('cel', "KEGG", 'kegg')
// 
// 
// my_kegg_data1 = clusterProfiler:::prepare_KEGG('cel', "KEGG", 'kegg')
// my_kegg_data1$EXTID2PATHID %>% head
// my_kegg_data$EXTID2PATHID %>% head
// 
// my_kegg_data1 = clusterProfiler:::prepare_KEGG('cel', "KEGG")
//  my_kegg_data1 = clusterProfiler:::prepare_KEGG('cel', 'KEGG', "ENTREZID")
// 
// res <- clusterProfiler:::enricher_internal(gene_set_entrez,
//                           pvalueCutoff  = 1,
//                           pAdjustMethod = 'BH',
//                           universe      = universe,
//                           qvalueCutoff  = 1,
//                           USER_DATA = my_kegg_data1)
// detach("package:clusterProfiler", unload=TRUE)
// library(clusterProfiler)



process Enrichment__computing_genes_self_overlaps {
  tag "${key}"

  label "r_basic"

  when:
    params.do_genes_self_enrichment && params.do_any_enrichment

  input:
    set key, data_type, file(lgenes), file('all_gene_sets/*') \
      from DA_genes_for_computing_genes_self_overlaps_3

  output:
    set key, data_type, file("*__counts.csv") \
      into Genes_self_overlaps_for_computing_pvalues

  shell:
    '''
    #!/usr/bin/env Rscript

    library(magrittr)

    key = '!{key}'
    lgenes = readRDS('!{lgenes}')


    DA = lgenes$DA
    NDA = lgenes$NDA

    all_gene_sets = list.files('all_gene_sets')
    all_gene_sets %<>% setNames(., gsub('__genes.rds', '', .))
    all_DA = purrr::map(all_gene_sets, 
      ~readRDS(paste0('all_gene_sets/', .x))$DA)

    df = purrr::imap_dfr(all_DA, function(genes_tgt, target){
      tgt = target
      tot_tgt = length(genes_tgt)
      tot_da = length(DA)
      ov_da = length(intersect(DA, genes_tgt))
      tot_nda = length(NDA)
      ov_nda = length(intersect(NDA, genes_tgt))
      data.frame(tgt, tot_tgt, tot_da, ov_da, tot_nda, ov_nda)
    })

    write.csv(df, paste0(key, '__genes_self__counts.csv'), row.names = F)

    '''
}

Overlap_tables_channel = Overlap_tables_channel
  .mix(Genes_self_overlaps_for_computing_pvalues)





DA_regions_with_bg_for_computing_peaks_overlaps_1
  .into{ 
    DA_regions_with_bg_for_computing_peaks_overlaps_2
    DA_regions_with_bg_for_computing_peaks_overlaps_3
  }
DA_regions_channel = DA_regions_with_bg_for_computing_peaks_overlaps_2
  .map{ it[1] }
  .toList()
  .map{ [ 'peaks_self'  , it ] }
  .dump(tag:'peaks_self')

// filtering the selected chromatin state file with the parameter 
// "params.chromatin_state_1"
Chrom_states_channel = 
  Channel.fromPath( "${params.chromatin_state_1}/*bed" )
  .toList()
  .map{ [ 'chrom_states', it ] }
  .dump(tag:'chrom_states')

// filtering the selected CHIP ontology group with the parameter 
// "params.chip_ontology"
CHIP_channel_all = Channel.fromPath( "${params.encode_chip_files}/*bed" )

Channel
  .fromPath( params.chip_ontology_groups )
  .splitCsv( header: false, sep: '\t' )
  .filter{ group, chip_bed -> group == params.chip_ontology }
  .map{ [it[1], it[0] ] }
  .map{ [it[0] ] }
  .dump(tag:'chip_ontology')
  .set{ chip_files_to_keep }

CHIP_channel_all
  .map{ [ it.name, it ] }
  .join(chip_files_to_keep)
  .map{ it[1] }
  .toList()
  .map{ [ 'CHIP', it ] }
  .set{ CHIP_channel_filtered }


Bed_regions_to_overlap_with = Channel.empty()
if(params.do_any_enrichment){
  
  if( params.do_peaks_self_enrichment ) Bed_regions_to_overlap_with = Bed_regions_to_overlap_with.mix(DA_regions_channel)
  if( params.do_chrom_state_enrichment ) Bed_regions_to_overlap_with = Bed_regions_to_overlap_with.mix(Chrom_states_channel)
  if( params.do_chip_enrichment ) Bed_regions_to_overlap_with = Bed_regions_to_overlap_with.mix(CHIP_channel_filtered)
  
}


DA_regions_with_bg_for_computing_peaks_overlaps_3
  // format: key (ET__PA__FC__TV__COMP), DA_regions, all_regions
  .combine(Bed_regions_to_overlap_with)
  // format: key, DA_regions, all_regions, data_type, bed_files
  .map{ [ it[0,3].join('__'), it[3], it[1], it[2], it[4] ] }
  .dump(tag:'bed_overlap')
  // format: key (ET__PA__FC__TV__COMP__DT), data_type, DA_regions, 
  //        all_regions, bed_files
  .set{ DA_regions_with_bg_and_bed_for_computing_peaks_overlaps }


process Enrichment__computing_peaks_overlaps {
  tag "${key}"

  label "samtools_bedtools_perl"

  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
    set key, data_type, file(DA_regions), file(all_regions), file("BED_FILES/*") \
      from DA_regions_with_bg_and_bed_for_computing_peaks_overlaps

  output:
    set key, data_type, file("*__counts.csv") \
      into Peaks_self_overlaps_for_computing_pvalues

  shell:
    '''

    DB=!{DA_regions}
    ALL=!{all_regions}
    KEY=!{key}


    intersectBed -v -a ${ALL} -b ${DB} > not_diffbound.bed
    NDB="not_diffbound.bed"

    tot_da=`wc -l < ${DB}`
    tot_nda=`wc -l < ${NDB}`

    OUTPUT_FILE="${KEY}__counts.csv"

    echo "tgt, tot_tgt, tot_da, ov_da, tot_nda, ov_nda" > $OUTPUT_FILE

    BEDS=($(ls BED_FILES))

    for BED1 in ${BEDS[@]}
      do
      	BED=BED_FILES/$BED1
        tgt=`basename ${BED} .bed`
        tgt=`basename ${tgt} __regions`
        tot_tgt=`wc -l < ${BED}`
        intersectBed -u -a ${DB} -b ${BED} > overlap_DB.tmp
        ov_da=`wc -l < overlap_DB.tmp`
        intersectBed -u -a ${NDB} -b ${BED} > overlap_NDB.tmp
        ov_nda=`wc -l < overlap_NDB.tmp`
        echo \
          "${tgt}, ${tot_tgt}, ${tot_da}, ${ov_da}, ${tot_nda}, ${ov_nda}" \
          >> $OUTPUT_FILE
      done

    '''
}

Overlap_tables_channel = Overlap_tables_channel
  .mix(Peaks_self_overlaps_for_computing_pvalues)


// note: we need to save tmp files (overlap_DB.tmp and overlap_NDB.tmp, 
// otherwise it crashes when there is zero overlap)

// note: here we use the option intersectBed -u instead of -wa to just indicate 
// if at least one chip peak overlap with a given atac seq peak. Thus the number 
// of overlap cannot be higher than the number of atac peaks.





DA_regions_with_bg_for_computing_motifs_overlaps_1
  // key (ET__PA__FC__TV__COMP), DA_regions, all_regions
  .map{ [ "${it[0]}__motifs", "motifs", it[1], it[2] ] }
  .dump(tag:'peaks_for_homer')
  // format: key (ET__PA__FC__TV__COMP__DT), data_type, DA_regions, all_regions
  .set{ DA_regions_with_bg_for_computing_motifs_overlaps_2 }



process Enrichment__computing_motifs_overlaps {
  tag "${key}"

  label "homer"
  cpus params.homer__nb_threads

  publishDir path: "${out_processed}/3_Enrichment/${data_type}/${key}", 
             mode: "${pub_mode}"

  when: 
    params.do_motif_enrichment && params.do_any_enrichment

  input:
    set key, data_type, file(DA_regions), file(all_regions) \
      from DA_regions_with_bg_for_computing_motifs_overlaps_2

  output:
    file("**")
    set key, data_type, file("*__homer_results.txt") optional true \
      into Motifs_overlaps_for_reformatting_results

  shell:
    '''

    findMotifsGenome.pl !{DA_regions} !{params.homer_genome} "." \
      -size given \
      -p !{params.homer__nb_threads} \
      -bg !{all_regions} \
      -mknown !{params.pwms_motifs} \
      -nomotif

    FILE="knownResults.txt"
    if [ -f $FILE ]; then
       mv $FILE "!{key}__homer_results.txt"
    fi

    '''
}

// note: homer automatically removes overlapping peaks between input and bg, 
// so it doesn't matter that we don't separtate them.



process Enrichment__reformatting_motifs_results {
  tag "${key}"

  label "r_basic"

  input:
    set key, data_type, file(motifs_results) \
      from Motifs_overlaps_for_reformatting_results

  output:
    set key, data_type, file("*__counts.csv") \
      into Formatted_motifs_results_for_computing_pvalues

  shell:
    '''
    #!/usr/bin/env Rscript

    library(magrittr)

    key = '!{key}'
    filename = '!{motifs_results}'


    df = read.csv(file = filename, sep = '\t', stringsAsFactors = F)

    total = purrr::map_chr(strsplit(names(df), 'Motif.of.')[c(6,8)], 2) %>% 
              gsub('.', '', ., fixed = T) %>% as.integer

    names(df) = c('tgt', 'consensus', 'pvalue', 'log_pval', 'qval', 'ov_da', 
                  'pt_da', 'ov_nda', 'pt_nda')
    df$tot_da  = total[1]
    df$tot_nda = total[2]
    df %<>% dplyr::select(tgt, tot_da, ov_da, tot_nda, ov_nda, consensus)

    ov_da_too_high = which(df$ov_da > df$tot_da)
    if(length(ov_da_too_high) > 0) df$ov_da[ov_da_too_high] = 
      df$tot_da[ov_da_too_high]
    ov_nda_too_high = which(df$ov_nda > df$tot_nda)
    if(length(ov_nda_too_high) > 0) df$ov_nda[ov_nda_too_high] = 
      df$tot_nda[ov_nda_too_high]

    write.csv(df, paste0(key, '__motifs__counts.csv'), row.names = F)

    '''
}

Overlap_tables_channel = Overlap_tables_channel
  .mix(Formatted_motifs_results_for_computing_pvalues)

// Overlap_tables_channel = Overlap_tables_channel.view()

// note : the "wrong_entries" variable is used because ov_nda is sometimes 
// slightly higher (by a decimal) than tot_nda; i.e.: tot_nda = 4374 and 
// ov_nda = 4374.6. This makes the Fischer test crash later on. This change 
// is minimal so we just fix it like that.

// This one liner works and is cleaner but it fails in nextflow due to the 
// double escape string
// total = as.integer(stringr::str_extract_all(paste0(names(df), collapse = ' '),
//          "\\(?[0-9]+\\)?")[[1]])





Overlap_tables_channel = Overlap_tables_channel.dump(tag: 'overlap_tables')

process Enrichment__computing_enrichment_pvalues {
  tag "${key}"

  label "r_basic"
  
  publishDir path: "${out_processed}/3_Enrichment/Enrichment__${data_type}", 
             mode: "${pub_mode}", enabled: params.enrichment__save_df

  input:
    set key, data_type, file(df_count_csv) from Overlap_tables_channel

  output:
    file("*.rds") into Enrichment_results_for_plotting optional true
    set data_type, val("3_Enrichment"), file("*.rds") \
      into Enrichment_results_for_formatting_table optional true

  shell:
    '''
    #!/usr/bin/env Rscript

    library(magrittr)
    source('!{projectDir}/bin/get_chrom_states_names_vec.R')

    key = '!{key}'
    data_type = '!{data_type}'
    df1 = read.csv('!{df_count_csv}', stringsAsFactors = F)
    motifs_test_type = '!{params.motifs_test_type}'


    # computing pvalue and L2OR
    df = df1
    for(c1 in 1:nrow(df)){
      ov_da   =   df$ov_da[c1]
      tot_da  =  df$tot_da[c1]
      ov_nda  =  df$ov_nda[c1]
      tot_nda = df$tot_nda[c1]
      mat = rbind(c(ov_da, ov_nda), c(tot_da - ov_da, tot_nda - ov_nda))
      fisher_test = fisher.test(mat, alternative = 'two.sided')
      df$pval[c1] = fisher_test$p.value
      df$L2OR[c1] = log2(fisher_test$estimate)
      if(data_type == 'motifs' & motifs_test_type == 'binomial') {
        df$pval[c1] = binom.test(ov_da, tot_da, ov_nda / tot_nda, 
          alternative = 'two.sided')$p.value
      }
    }

    # adding padj and percentage of overlap
    df$padj = p.adjust(df$pval, method = 'BH')
    df %<>% dplyr::mutate(pt_da = ov_da  / tot_da  )
    df %<>% dplyr::mutate(pt_nda = ov_nda / tot_nda )
    df$pt_da %<>% {. * 100 } %>% round(2) %>% paste0(., '%')
    df$pt_nda %<>% {. * 100 } %>% round(2) %>% paste0(., '%')

    # renaming chromatin states
    if(data_type == 'chrom_states'){
      vec = get_chrom_states_names_vec(df$tgt)
      df$tgt %<>% vec[.]
    }

    # reordering columns
    df %<>% dplyr::select(tgt, pval, padj, L2OR, pt_da, ov_da, tot_da, 
      pt_nda, ov_nda, tot_nda, dplyr::everything())

    # adding a Gene Enrichment Type column for func_anno
    if(grepl('func_anno', data_type)) {
      GE = gsub('func_anno_', '', data_type)
      data_type1 = 'func_anno'
    } else { data_type1 = data_type }

    # sorting by padj and then overlap counts
    df %<>% dplyr::arrange(padj, desc(ov_da))

    # adding the key and saving for plots
    key_split = strsplit(key, '__')[[1]]
    key_df = as.data.frame(t(key_split[-length(key_split)]), 
      stringsAsFactors = F)
    cln = c('ET', 'PA', 'FC', 'TV', 'COMP')
    key_df %<>% set_colnames(cln)
    if(data_type1 == 'func_anno') key_df %<>% cbind(GE = GE, .)
    df2 = cbind(key_df, df)
    saveRDS(df2, paste0(key, '__enrich.rds'))

    '''
}


// => the input should be a df with these columns: tgt tot_da ov_da tot_nda ov_nda
// other extra columns are facultatory. These are: tot_tgt for bed_overlap, 
// consensus for motifs, geneIds for ontologies/pathways

// data_types can be either of: func_anno(BP|CC|MF|KEGG), genes_self, 
// peaks_self, chrom_states, CHIP, motifs

// df$tgt %>% .[. %in% names(vec)]
// df$tgt %>% .[!. %in% names(vec)]



Formatting_csv_tables_channel = Formatting_csv_tables_channel
  .mix(Enrichment_results_for_formatting_table)




// note DT stands for Data_Type. It can be either 'genes_(BP|CC|MF|KEGG)', 
// 'CHIP', 'motif' or 'chrom_states'
// comp_order has this format: 
// Eri1Daf2_vs_Eri1|Eri1Eat2_vs_Eri1|Eri1Isp1_vs_Eri1|...


Enrichment_results_for_plotting
  // format: rds_file (key__enrich.rds)
  .flatten()
  // we need to flatten since the genes_sets
  .map{ [ it.name.replaceFirst(~/__enrich.rds/, ''), it ] }
  // format: key (ET__PA__FC__TV__COMP__DT), rds_file
  .tap{ Enrichment_results_for_plotting_barplots_1 }
  .combine(comparisons_grouped_for_heatmap)
  // format: key, rds_file, GRP, comp_order
  .dump(tag:'enrichment')
  .map{ [ it[0].split('__')[4], it ].flatten() }
  // format: COMP, key, rds_file, GRP, comp_order
  .filter{ it[0] in it[4].split('\\|') }
  // keeping only COMP that are in the group
  .map{ [ it[1].split('__'), it[2..4] ].flatten() }
  // format: ET, PA, FC, TV, COMP, DT, rds_file, GRP, comp_order
  .map{ [ it[0, 1, 3, 7, 5].join('__'), it[5, 8, 3, 6] ].flatten() }
  // format: key (ET__PA__TV__GRP__DT), DT, comp_order, TV, rds_file
  .groupTuple(by: [0, 1, 2, 3])
  // format: key, DT, comp_order, TV, rds_files
  .filter{ it[4].size() > 1 }
  // .map{ [ it[0], it[1], it[1].replaceAll('_(KEGG|GO_(BP|CC|MF))', ''), 
  .map{ [ it[0], it[1], it[1].replaceAll('_(KEGG|BP|CC|MF)', ''), 
          it[2], it[3], it[4] ]  }
  // format: key (ET__PA__FC__TV__COMP__DT), DT, DT_short, comp_order, TV, rds_files
  .map{ [ it[0], it[1], params.heatmaps_params[it[2]], 
            params.padj_breaks[it[2]], params.heatmaps_filter[it[2]], 
            it[3], it[4], it[5] ] }
  // format: key, DT, plots_params, padj_breaks, filters, comp_order, TV, rds_files
  .dump(tag:'heatmap')
  .set{ Enrichment_results_for_plotting_heatmaps }




Enrichment_results_for_plotting_barplots_1
  // format: key (ET__PA__FC__TV__COMP__DT), rds_file
  .map{ [ it[0], it[0].split('__')[5], it[1] ]  }
  // format: key (ET__PA__FC__TV__COMP__DT), DT, rds_file
  // .map{ [ it[0], it[1], it[1].replaceAll('_(KEGG|GO_(BP|CC|MF))', ''), it[2] ]  }
  .map{ [ it[0], it[1], it[1].replaceAll('_(KEGG|BP|CC|MF)', ''), it[2] ]  }
  // format: key (ET__PA__FC__TV__COMP__DT), DT, DT_short, rds_file
  .dump(tag: 'barplot')
  .map{ [ it[0], it[1], params.barplots_params[it[2]], params.padj_breaks[it[2]], it[3]]}
  // format: key (ET__PA__FC__TV__COMP__DT), DT, plots_params, padj_breaks, rds_file
  .set{ Enrichment_results_for_plotting_barplots_2 }



process Figures__making_enrichment_barplots {
  tag "${key}"

  label "figures"

  publishDir path: "${out_fig_indiv}/3_Enrichment/Barplots__${data_type}", 
             mode: "${pub_mode}"

  input:
    set key, data_type, plots_params, padj_breaks, 
      file(res_gene_set_enrichment_rds) \
      from Enrichment_results_for_plotting_barplots_2

  output:
    set val("Barplots__${data_type}"), val("3_Enrichment"), file("*.pdf") \
      optional true into Barplots_for_merging_pdfs

  shell:
    
    '''
    #!/usr/bin/env Rscript


    library(ggplot2)
    library(grid)
    library(gridExtra)
    library(RColorBrewer)
    library(magrittr)
    library(data.table)

    source('!{projectDir}/bin/get_new_name_by_unique_character.R')
    source('!{projectDir}/bin/functions_pvalue_plots.R')
    
    key         = '!{key}'
    data_type   = '!{data_type}'
    df1         = readRDS('!{res_gene_set_enrichment_rds}')
    padj_breaks = !{padj_breaks}
    plot_params = !{plots_params}
    

    # getting parameters
    if(grepl('func_anno', data_type)) data_type = 'func_anno'
    df             = df1
    padj_threshold = plot_params[1] %>% as.numeric
    signed_padj    = plot_params[2] %>% as.logical
    add_var        = plot_params[3] %>% as.character
    add_number     = plot_params[4] %>% as.logical
    max_characters = plot_params[5] %>% as.integer
    max_terms      = plot_params[6] %>% as.integer


    # quitting if there are no significant results to show
    if(all(df$padj > padj_threshold)) quit(save = 'no')

    # removing the genes_id column from func_anno enrichments 
    # (for easier debugging)
    df %<>% .[, names(.) != 'genes_id']

    # adding the loglog and binned padj columns
    df %<>% getting_padj_loglog_and_binned(signed_padj, padj_breaks)

    # adding the yaxis_terms column with shortened and unique names
    df$yaxis_terms = df$tgt %>% get_shorter_names(max_characters)

    # selecting lowest pvalues
    df = df[seq_len(min(nrow(df), max_terms)), ]
    df$yaxis_terms %<>% factor(., levels = rev(.))

    bed_data_types = c('CHIP', 'chrom_states', 'peaks_self', 'genes_self')
    is_bed_overlap = data_type %in% bed_data_types
    if(is_bed_overlap){
      xlab = paste0('Overlap (DA: ', df$tot_da[1], ', NDA: ', df$tot_nda[1], ')')
    } else {
      xlab = paste0('Overlap (DA: ', df$tot_da[1], ')')
    }

    p1 = ggplot(df, aes(x = yaxis_terms, y = ov_da)) + coord_flip() + 
      geom_bar(stat = 'identity') + ggtitle(key) + theme_bw() + 
      theme(axis.title.y = element_blank(), 
        axis.text = element_text(size = 11, color = 'black'), 
        plot.title = element_text(hjust = 0.9, size = 10), 
        legend.text = element_text(size = 7)
       ) + ylab(xlab)

    point_size = scales::rescale(c(nrow(df), seq(0, 30, len = 5)), c(6, 3))[1]
    p_binned = get_plot_binned(p1, signed_padj = signed_padj, 
                    add_var = add_var, add_number  = add_number, 
                    point_size  = point_size)

    pdf(paste0(key, '__barplot.pdf'), paper = 'a4r')
      print(p_binned)
    dev.off()

    '''
}

// # signed_padj = ifelse(data_type %in% c('CHIP', 'chrom_states'), T, F)

Merging_pdfs_channel = Merging_pdfs_channel.mix(Barplots_for_merging_pdfs
  .groupTuple(by: [0, 1]))


process Figures__making_enrichment_heatmap {
  tag "${key}"

  label "figures"

  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir path: "${out_fig_indiv}/3_Enrichment/Heatmaps__${data_type}", 
             mode: "${pub_mode}"

  input:
    set key, data_type, plots_params, padj_breaks, filters, comp_order, 
      tv, file('*') from Enrichment_results_for_plotting_heatmaps

  output:
    set val("Heatmaps__${data_type}"), val("3_Enrichment"), file("*.pdf") \
      optional true into Heatmaps_for_merging_pdfs

  shell:
    '''
    #!/usr/bin/env Rscript


    library(magrittr)
    library(ggplot2)
    library(RColorBrewer)
    library(data.table)

    source('!{projectDir}/bin/get_new_name_by_unique_character.R')
    source('!{projectDir}/bin/get_chrom_states_names_vec.R')
    source('!{projectDir}/bin/functions_pvalue_plots.R')
    source('!{projectDir}/bin/functions_grouped_plot.R')

    key         = '!{key}'
    comp_order  = '!{comp_order}'
    data_type   = '!{data_type}'
    padj_breaks = !{padj_breaks}
    plot_params = !{plots_params}
    filters     = !{filters}
    seed        = '!{params.heatmaps__seed}'


    # getting parameters
    if(grepl('func_anno', data_type)) data_type = 'func_anno'
    
    padj_threshold  = plot_params[1] %>% as.numeric
    signed_padj     = plot_params[2] %>% as.logical
    add_var         = plot_params[3] %>% as.character
    add_number      = plot_params[4] %>% as.logical
    max_characters  = plot_params[5] %>% as.integer
    up_down_pattern = plot_params[6] %>% as.character
    
    if(data_type %in% c('CHIP', 'motifs', 'func_anno')){
      n_total              = filters[1] %>% as.integer
      n_shared             = filters[2] %>% as.integer
      n_unique             = filters[3] %>% as.integer
      remove_similar       = filters[4] %>% as.logical
      remove_similar_n     = filters[5] %>% as.integer
      agglomeration_method = filters[6] %>% as.character
      select_enriched      = filters[7] %>% as.logical
    }

    # loading, merging and processing data
    rds_files = list.files(pattern = '*.rds')
    ldf = lapply(rds_files, readRDS)
    df = do.call(rbind, ldf)

    # filtering table
    if(data_type %in% c('genes_self', 'peaks_self')){
      key1 = paste(df[1, c('ET', 'PA', 'TV')], collapse = '__')
      df$tgt_key = sapply(strsplit(df$tgt, '__'), 
                          function(x) paste(x[c(1,2,4)], collapse = '__'))
      df %<>% .[.$tgt_key %in% key1, ]
      df$tgt_key <- NULL
    }

    ## quitting if there are no significant results to show
    if(all(df$padj > padj_threshold)) quit(save = 'no')

    # adding the comp_FC
    df$comp_FC = apply(df[, c('COMP', 'FC')], 1, paste, collapse = '_') %>% 
                 gsub('_vs_', '_', .)

    # ordering the x-axis comparisons
    comp_order_levels = get_comp_order_levels(comp_order, up_down_pattern)
    comp_order1 = comp_order_levels %>% .[. %in% unique(df$comp_FC)]

    # adding the yaxis terms column
    df$yaxis_terms = df$tgt

    # adding loglog and binned padj columns
    df %<>% getting_padj_loglog_and_binned(signed_padj = signed_padj, padj_breaks)

    purrr__map_chr <- function(x, c1) lapply(x, function(y) y[c1]) %>% 
                      as.character

    # if plotting self overlap keeping only targets in the group
    if(data_type %in% c('genes_self', 'peaks_self')) {
      terms_levels = rev(comp_order1)
      strs = strsplit(df$tgt, '__')
      tgt_ET   = purrr__map_chr(strs, 1)
      tgt_PA   = purrr__map_chr(strs, 2)
      tgt_FC   = purrr__map_chr(strs, 3)
      tgt_TV   = purrr__map_chr(strs, 4)
      tgt_COMP = purrr__map_chr(strs, 5)
      ET = unique(df$ET)
      PA = unique(df$PA)
      TV = unique(df$TV)
      df$tgt_comp_FC = paste0(tgt_COMP, '_', tgt_FC) %>% gsub('_vs_', '_', .)
      df = subset(df, tgt_comp_FC %in% comp_FC & tgt_ET == ET & tgt_PA == PA)
      df$yaxis_terms = df$tgt_comp_FC
    }

    # reformatting df to a matrix
    dt = data.table::copy(df) %>% setDT
    dt[, signed_padj_log_log := sign(L2OR) * padj_loglog]
    mat_dt = dcast(dt, yaxis_terms ~ comp_FC, 
                    value.var = 'signed_padj_log_log', fill = get_pval_loglog(1))
    mat = as.matrix(mat_dt[,-1]) %>% set_rownames(mat_dt$yaxis_terms)

    # selecting and ordering the y-axis terms
    if(data_type %in% c('CHIP', 'motifs', 'func_anno')){

      sel_rows = which(rowSums(abs(mat)) != 0)
      sel_cols = which(colSums(abs(mat)) != 0)
      if(length(sel_rows) < 2 | length(sel_cols) < 2) quit(save = 'no')
      mat = mat[sel_rows, sel_cols]
      comp_order1 = comp_order_levels %>% .[. %in% colnames(mat)]
      
      terms_levels = select_y_axis_terms_grouped_plot(mat, n_total = n_total, 
        n_shared = n_shared, n_unique = n_unique, 
        remove_similar = remove_similar, remove_similar_n = remove_similar_n, 
        seed = seed, agglomeration_method = agglomeration_method, 
        select_enriched = select_enriched)

    }

    if(data_type == 'chrom_states') {
      vec = unique(df$tgt)
      vec_order = vec %>% gsub('.', '', ., fixed = T) %>% gsub(' .*', '', .) %>% 
                  as.integer %>% order
      terms_levels = vec[vec_order] %>% rev
    }

    # making the matrix of selected terms and with the right order of conditions
    mat_final = mat[terms_levels, comp_order1]
    nrows = nrow(mat_final) ; ncols = ncol(mat_final)

    # shortening the long terms names while keeping them unique
    rownames(mat_final) %<>% get_shorter_names(max_characters)

    # creating a data.frame with row and column indexes for plotting
    df_final = add_matrix_indexes_to_df(mat_final, df, nrows, ncols)

    # making and saving plots
    p1 = getting_heatmap_base(df_final, nrows, ncols, title = key, 
      cur_mat = mat_final)
    point_size = scales::rescale(c(nrows, seq(0, 40, len = 5)), c(3, 0.8))[1]
    p_binned = get_plot_binned(p1, signed_padj = signed_padj, add_var, 
      add_number, point_size = point_size)

    pdf(paste0(key, '__heatmap.pdf'))
      print(p_binned)
    dev.off()


    '''
}

// note that: gtools::mixedsort(df$tgt) would be simpler for chromatin states 
// if gtools was in the container

// signed_padj = data_type %in% c('CHIP', 'chrom_states')

Merging_pdfs_channel = Merging_pdfs_channel.mix(Heatmaps_for_merging_pdfs
  .groupTuple(by: [0, 1]))


// 
// signed_padj = F;  add_loglog = F; add_L2OR = F
// signed_padj = T;  add_loglog = T; add_L2OR = T
// 
// 
// ## Some explainations on the selection of terms to display (with the example 
// of the algorithm for genes)
// # 26 terms at max will be plotted since it is the maximum to keep a readable 
// plot
// # the first 6 slots will be attributed to terms that are the most shared 
// amoung groups (if any term is shared)
// # then the top_N term for each group will be selected, aiming at 20 terms max
// # then the lowest pvalues overall will fill the rest of the terms (since 
// shared terms will leave empty gaps)
// # Note that there should be at maximum 10 comparisons in each grouping plot 
// (10 * 2 for up and down means at least 1 pvalue for each comparison)






Merging_pdfs_channel = Merging_pdfs_channel.dump(tag: 'merging_pdf')

process Figures__merging_pdfs {
  tag "${file_name}"

  label "pdftools"

  publishDir path: "${out_fig_merge}/${out_path}", mode: "${pub_mode}"

  input:
    set file_name, out_path, file("*") from Merging_pdfs_channel

  output:
    file("*.pdf") optional true

  script:
    """

    pdfunite `ls *pdf | sort` ${file_name}.pdf

    """
}





Formatting_csv_tables_channel = Formatting_csv_tables_channel
  .dump(tag: 'csv_tables')

process Tables__formatting_csv_tables {
  tag "${out_folder}__${data_type}"

  label "r_basic"

  publishDir path: "${out_tab_indiv}/${out_folder}/${data_type}", 
    mode: "${pub_mode}", enabled: params.tables__save_csv

  input:
    set data_type, out_folder, file(rds_file) from Formatting_csv_tables_channel

  output:
    set data_type, out_folder, file('*.csv') \
      into Formatted_csv_tables_for_merging optional true
    set val("Tables_Individual/${out_folder}/${data_type}"), file('*.csv') \
      into Formatted_csv_tables_for_saving_Excel_tables optional true

  shell:
    '''
    #!/usr/bin/env Rscript

    library(magrittr)

    source('!{projectDir}/bin/get_formatted_table.R')

    data_type = '!{data_type}'
    rds_file = '!{rds_file}'

    v_fdr_thresholds = 
    eval(parse(text = '!{params.tables__v_fdr_thresholds}'))


    # reading
    df = readRDS(rds_file)

    # filtering
    data_type1 = data_type
    if(grepl('func_anno', data_type)) data_type1 = 'func_anno'
    fdr_threshold = v_fdr_thresholds[data_type1]
    df = subset(df, padj <= fdr_threshold)

    # formating
    df %<>% get_formatted_table

    # saving
    if(nrow(df) > 0) {
      output_file_name = paste0(gsub('.rds', '', rds_file), '.csv')
      write.csv(df, output_file_name, row.names = F)
    }

    '''
}

// => this process allows to reformat tables without breaking the cache


Exporting_to_Excel_channel = Exporting_to_Excel_channel
  .mix(Formatted_csv_tables_for_saving_Excel_tables)


Formatted_csv_tables_for_merging
  .groupTuple(by: [0, 1])
  .dump(tag: 'merge_tables')
  .set{ Formatted_tables_grouped_for_merging_tables }

process Tables__merging_csv_tables {
  tag "${out_folder}__${data_type}"

  label "r_basic"

  publishDir path: "${out_tab_merge}/${out_folder}", mode: "${pub_mode}", 
             enabled: params.tables__save_csv

  input:
    set data_type, out_folder, file(csv_file) \
      from Formatted_tables_grouped_for_merging_tables

  output:
    set val("Tables_Merged/${out_folder}"), file("*.csv") \
      into Merged_csv_tables_for_Excel optional true

  shell:
    '''
    #!/usr/bin/env Rscript

    library(magrittr)
    library(dplyr)
    source('!{projectDir}/bin/get_formatted_table.R')

    data_type = '!{data_type}'


    # merging tables, 
    all_files = list.files(pattern = '*.csv')
    ldf = lapply(all_files, read.csv, stringsAsFactors = F, as.is = T)
    df = do.call(rbind, ldf)

    # formatting and saving merged table
    df %<>% get_formatted_table
    write.csv(df, paste0(data_type, '.csv'), row.names = F)

    '''
}

Exporting_to_Excel_channel = Exporting_to_Excel_channel
  .mix(Merged_csv_tables_for_Excel)



Exporting_to_Excel_channel = Exporting_to_Excel_channel.dump(tag: 'excel')

process Tables__saving_excel_tables {
  tag "${csv_file}"

  // label "openxlsx" => sh: : Permission denied ; Error: zipping up workbook 
  // failed. Please make sure Rtools is installed or a zip application is 
  // available to R.

  label "differential_abundance"

  publishDir path: "${res_dir}/${out_path}", mode: "${pub_mode}" 

  when: 
    params.tables__save_excel

  input:
    set out_path, file(csv_file) from Exporting_to_Excel_channel

  output:
    file("*.xlsx")

  shell:
    '''
    #!/usr/bin/env Rscript

    library(openxlsx)

    csv_file = '!{csv_file}'
    excel__add_conditional_formatting = !{params.excel__add_conditional_formatting}
    excel__max_width = !{params.excel__max_width}


    options(digits = 1)

    df = read.csv(csv_file, stringsAsFactors = T, as.is = T)
    output_file_name = paste0(gsub('.csv', '', csv_file), '.xlsx')

    nms = names(df)

    class(df$pval) = 'scientific'
    class(df$padj) = 'scientific'

    if('pt_da' %in% nms){
      class(df$pt_da) = 'percentage'
      class(df$pt_nda) = 'percentage'
      L2OR_Inf_up   = which(df$L2OR == 'Inf')
      L2OR_Inf_down = which(df$L2OR == '-Inf')
      L2OR_not_Inf  = which(abs(df$L2OR) != 'Inf')
      df$L2OR[L2OR_Inf_up]   = 1e99
      df$L2OR[L2OR_Inf_down] = -1e99
    }

    names_colors_1  = c( 'filter',  'target',    'fold',  'pvalue',   'da'   )
    # color type           red        green     purple     orange     gold     
    header_colors_1 = c('#963634', '#76933c', '#60497a', '#e26b0a', '#9d821f')
    body_colors_1   = c('#f2dcdb', '#ebf1de', '#e4dfec', '#fde9d9', '#f1edcb')

    names_colors_2  = c(  'nda',    'other',    'gene', 'coordinate')
    # color type       pale_blue    grey     darkblue   darkolive
    header_colors_2 = c('#31869b', '#808080', '#16365c',  '#494529' )
    body_colors_2   = c('#daeef3', '#f2f2f2', '#c5d9f1',  '#ddd9c4' )

    names_colors = c(names_colors_1, names_colors_2)
    header_colors = c(header_colors_1, header_colors_2)
    body_colors = c(body_colors_1, body_colors_2)

    names(header_colors) = names_colors
    names(body_colors)   = names_colors


    get_nms_type <- function(nms){
      nms_coordinates = c('chr','start', 'end',	'width', 'strand')
      nms_coordinates = c(nms_coordinates, paste0('gene_', nms_coordinates))

      if(nms %in% c('GE', 'ET', 'PA', 'FC', 'TV', 'COMP')) return('filter') else
      if(nms %in% c('gene_name', 'gene_id', 'entrez_id'))  return('gene') else
      if(nms %in% c('pval', 'padj'))                       return('pvalue') else
      if(nms %in% c('L2FC', 'L2OR'))                       return('fold') else
      if(nms %in% c('tgt'))                                return('target') else
      if(nms %in% c('pt_da', 'tot_da', 'ov_da'))           return('da') else
      if(nms %in% c('pt_nda', 'tot_nda', 'ov_nda'))        return('nda') else
      if(nms %in% nms_coordinates)                         return('coordinate') else
                                                           return('other')
    }
    nms_types = sapply(names(df), get_nms_type)
    nms_color_header = unname(header_colors[nms_types])
    nms_color_body   = unname(body_colors  [nms_types])

    sheet = 1
    cols = seq_len(ncol(df))
    rows = seq_len(nrow(df) + 1)

    # create the workbook
    wb = write.xlsx(df, output_file_name, borders = 'rows', keepNA = F)

    # add filter, set width and height
    addFilter(wb, sheet, 1, cols)
    setRowHeights(wb, sheet, 1, heights = 50)
    widths = apply(df, 2, function(x) {
      if(all(is.na(x))) return(5)
      width = max(nchar(x), na.rm = T) + 2.5
      width = ifelse(width > excel__max_width, excel__max_width, width)
      return(width)
    })
    setColWidths(wb, sheet, cols, widths = widths)

    for(col in cols) {
      col_nm = nms[col]
      halign = ifelse('GE' %in% nms & col_nm %in% c('tgt', 'genes_id'), 
                      'left', 'center')
      header_style = createStyle(fontColour = '#ffffff', 
        fgFill = nms_color_header[col], halign = halign, valign = 'center', 
        textDecoration = 'Bold', border = 'TopBottomLeftRight', wrapText = T)
      addStyle(wb, sheet, header_style, rows = 1, col)

      body_style = createStyle(halign = halign, valign = 'center', 
        fgFill = nms_color_body[col])
      addStyle(wb, sheet, body_style, rows = rows[-1], col)

      if(excel__add_conditional_formatting){
        if(col_nm == 'padj') conditionalFormatting(wb, sheet, cols = col, 
          rows = rows[-1], type = 'colourScale', 
          style = c('#e26b0a', '#fde9d9')) 
        vec = df[[col]]
        if(col_nm == 'L2FC') {
          conditionalFormatting(wb, sheet, cols = col, rows = rows[-1], 
            type = 'colourScale', 
            style = c(blue = '#6699ff', white = 'white', red = '#ff7c80'), 
            rule = c(min(vec), 0, max(vec)))
        }
        if(col_nm == 'L2OR') {
          if(length(L2OR_not_Inf) > 0) {
            vec1 = vec[L2OR_not_Inf]
            vec1 = vec1[!is.na(vec1)]
            conditionalFormatting(wb, sheet, cols = col, 
              rows = L2OR_not_Inf + 1, type = 'colourScale', 
              style = c(blue = '#6699ff', white = 'white', red = '#ff7c80'), 
              rule = c(min(vec1), 0, max(vec1)))
          }
          if(length(L2OR_Inf_up) > 0) addStyle(wb, sheet, 
            createStyle(fgFill = c(lightblue = '#ff7c80'), halign = 'center', 
            valign = 'center'), rows = L2OR_Inf_up + 1, col)
          if(length(L2OR_Inf_down) > 0) addStyle(wb, sheet, 
            createStyle(fgFill = c(lightred = '#6699ff'), halign = 'center', 
            valign = 'center'), rows = L2OR_Inf_down + 1, col)
        }
      }
    }

    # save the final workbook
    saveWorkbook(wb, output_file_name, overwrite = TRUE)

    '''
}

// documentation for openxlsx
// # https://www.rdocumentation.org/packages/openxlsx/versions/4.1.0.1
// # https://ycphs.github.io/openxlsx/articles/Introduction.html




////////////////////////////////////////////////////////////////////////////
//// THE END MY FRIEND


 workflow.onComplete
 {
  println ""
  println "Workflow completed on: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'Succeeded' : 'Failed' }"
  println "Workflow Duration: $workflow.duration"
  println ""
 }
 
 
 
