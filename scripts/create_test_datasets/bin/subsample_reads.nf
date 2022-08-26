

random_seed    = 38
thousand_reads = params.thousand_reads
experiment     = params.experiment
specie         = params.specie

prepro_dir = "preprocessing/${specie}"
in_dir     = "${prepro_dir}/fastq"
out_dir    = "${specie}/fastq/${experiment}/"

n_K_reads = thousand_reads + 'K_reads'
total_reads = thousand_reads * 1000

paired_end = Channel
  .fromFilePairs("${in_dir}/*_{1,2}.fastq.gz")
  .map { it -> [ type: 'paired_end', id: it[0], path: it[1] ] }
  // .view()
  
single_end = Channel
  .fromPath("${in_dir}/*")
  .filter { !(it =~ /_(1|2).fastq.gz/) }
  .map { it -> [ it.name.take(it.name.lastIndexOf('.')), it] }
  .map { it -> [ it[0].take(it[0].lastIndexOf('.')), it[1]] }
  .map { it -> [ type: 'single_end', id: it[0], path: it[1] ] }
  // .view()

reads = paired_end.mix(single_end)
  .filter{ it.id =~ "${experiment}_.*"}
  // .view()
  

process seqtk {
  tag "${id} ${type} ${n_K_reads}"

  container 'quay.io/biocontainers/seqtk:1.3--hed695b0_2'

	publishDir "${out_dir}", mode: 'copy', overwrite:true

	input:
	set val(type), val(id), file(fastq_file) from reads
	val total_reads 

	output:
	file '*.fastq.gz'

  script:

    if( type == 'single_end' ) {
  	"""

      seqtk sample -s ${random_seed} ${fastq_file} ${total_reads} | gzip > "sample_"${n_K_reads}_${id}".fastq.gz"

    """
  }
    else {
    """
      seqtk sample -s ${random_seed} ${id}_2.fastq.gz ${total_reads} | gzip > "sample_"${n_K_reads}_${id}"_1.fastq.gz"
      seqtk sample -s ${random_seed} ${id}_1.fastq.gz ${total_reads} | gzip > "sample_"${n_K_reads}_${id}"_2.fastq.gz"
    """
  }

}



workflow.onComplete
{
 println ""
 println "Workflow completed on: $workflow.complete"
 println "Execution status: ${ workflow.success ? 'Succeeded' : 'Failed' }"
 println "Workflow Duration: $workflow.duration"
 println ""
}


