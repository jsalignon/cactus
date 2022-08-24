

// nextflow run jsalignon/cactus/download/download_references.nf --specie worm  -latest -r main
// nextflow run $cactus/download/download_references.nf --specie worm


process download_references {
  tag params.specie

  container = params.skewer_pigz

  publishDir path: "${params.reference_dir}", mode: "link"

  input:

  output:
    file("${params.specie}")  
    
  script:
  def figshare_path = "https://ndownloader.figshare.com/files"
  def local_file    = "{params.specie}_test_dataset.tar.gz"
  
  """
      local_file=${params.specie}_references.tar.gz
      wget ${figshare_path}/${params.reference_file}?access_token=${params.figshare_token} -O ${local_file}
      local_md5sum=\$(md5sum ${local_file} | awk '{print \$1}')
      
      if [[ "\$local_md5sum" == "${params.reference_md5sum}" ]] ; then
        pigz --decompress --keep --recursive --stdout --processes ${params.threads} ${local_file} | tar -xvf
      else
        echo "md5sum is wrong for file ${local_file}"
    	fi
       
  """

}

// tar --use-compress-program="pigz -p ${params.threads} -k -r" -xvf ${local_file} => the pigz option is not avaialble in busybox's tar
// pigz -d -k -r -c -p 3 $local_file | tar -xvf - // same command with abbreviations

