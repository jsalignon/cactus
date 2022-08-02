

thousand_reads = 100
random_seed    = 38

indir          = params.indir
thousand_reads = params.thousand_reads
thousand_reads = params.thousand_reads
experiment     = params.experiment


n_K_reads = thousand_reads + 'K_reads'
total_reads = thousand_reads * 1000

paired_end = Channel
  .fromFilePairs("${indir}/*_{1,2}.fastq.gz")
  .map { it -> [ type: 'paired_end', id: it[0], path: it[1] ] }
  // .view()
  
single_end = Channel
  .fromPath("${indir}/*")
  .filter { !(it =~ /_(1|2).fastq.gz/) }
  .map { it -> [ it.name.take(it.name.lastIndexOf('.')), it] }
  .map { it -> [ it[0].take(it[0].lastIndexOf('.')), it[1]] }
  .map { it -> [ type: 'single_end', id: it[0], path: it[1] ] }
  // .view()

reads = paired_end.mix(single_end)
  // .filter{ it.id == "mrna_SRX11708668_SRR15406505_R2"}
  .filter{ it.id =~ /"atac_.*"/}
  .view()
  // .filter{read_type, experiment_id, fastq_file -> id experiment_id =~ /"${experiment}_.*"/}
  // .filter{ it[0] =~ /.*ATAC.*/ }

// .filter{ it[0] =~ /.*ATAC.*/ }
  
// process seqtk {
//   tag "${id} ${type}"
// 
//   container 'quay.io/biocontainers/seqtk:1.3--hed695b0_2'
// 
// 	publishDir indir + '_' + n_K_reads, mode: 'copy', overwrite:true
// 
// 	input:
// 	set val(type), val(id), file(fastq_file) from reads
// 	val total_reads 
// 
// 	output:
// 	file '*.fastq.gz'
// 
//   script:
// 
//     if( type == 'single_end' ) {
//   	"""
// 
//       seqtk sample -s ${random_seed} ${fastq_file} ${total_reads} | gzip > "sample_"${n_K_reads}_${id}".fastq.gz"
// 
//     """
//   }
//     else {
//     """
//       seqtk sample -s ${random_seed} ${id}_2.fastq.gz ${total_reads} | gzip > "sample_"${n_K_reads}_${id}"_1.fastq.gz"
//       seqtk sample -s ${random_seed} ${id}_1.fastq.gz ${total_reads} | gzip > "sample_"${n_K_reads}_${id}"_2.fastq.gz"
//     """
//   }
// 
// }



