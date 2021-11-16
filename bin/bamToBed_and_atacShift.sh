#!/bin/bash

function getInsertSites {
    	# Files and directories to write to
    	INSERT_SITES_BED_FILE="${2}_1bp_shifted_reads.bed"
    	INSERT_SITES_BAM_PREFIX="${2}_1bp_shifted_reads"
    	INSERT_SITES_BAM_FILE="${2}_1bp_shifted_reads.bam"

    	bamToBedInsertLengthAsScore "${1}" | adjustBedTn5 | fivePrimeOfBed6Cols - > "${INSERT_SITES_BED_FILE}"
    	# Convert back to a bam, but now just the 1bp
    	bedToBam -i "${INSERT_SITES_BED_FILE}" -g "${3}" > temp.bam
    	samtools sort temp.bam -o "${INSERT_SITES_BAM_PREFIX}.bam"
    	rm temp.bam
    	samtools index "${INSERT_SITES_BAM_FILE}"

    	# And then a little clean up to save space
    	# gzip "${INSERT_SITES_BED_FILE}"
    }

    function bamToBedInsertLengthAsScore {
    	bamToBed -i "${1}"> temp.bed
    	samtools view  "${1}" | perl -lane 'print $F[0]."\t".$F[8];' > insert.temp

    	paste temp.bed insert.temp | \
    	perl -lane '$name=(split/\//,$F[3])[0]; if ($name =~ /$F[6]/){ print $F[0]."\t".$F[1]."\t".$F[2]."\t".$F[3]."\t".$F[7]."\t".$F[5]; }else{ print STDERR "NAMES DONT MATCH for $name $F[6]!!!!";}'
    	rm temp.bed insert.temp
    }

    function adjustBedTn5 {
    	awk -F $'\t' 'BEGIN {OFS = FS}{if ($6 == "+") { $2 = $2 + 4 } else if ($6 == "-") {$3 = $3 - 5} print $0}' "${1}"
    }


    function fivePrimeOfBed6Cols {
    	cat "${1}" | awk 'BEGIN{OFS="\t"}{if($6 == "-") $2=$3-1; print $1, $2, $2+1, $4, $5, $6}' | sortBed -i -
    }

  getInsertSites "${1}" "${2}" "${3}"

# $1 = bam file   $2 = path/prefix of new file   $3 = genome size

