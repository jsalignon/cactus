

// nextflow run /home/jersal/workspace/cactus/download/download_references.nf --specie worm

process download_references {

  container = params.skewer_pigz

  publishDir path: "${params.reference_dir}", mode: "link"

  input:

  output:
    file("*")  
    
  script:
  def figshare_path = "https://ndownloader.figshare.com/files"
  
  """
      local_file=${params.reference_dir}/${params.specie}_references.tar.gz
      wget ${figshare_path}/${params.reference_file}?access_token=${params.figshare_token} -O \$local_file
      local_md5sum=\$(md5sum \$local_file | awk '{print \$1}')
      
      if [[ "\$local_md5sum" == "$params.reference_md5sum" ]] ; then
        tar -xvf \$local_file
      else
        echo "md5sum is wrong for file \$local_file"
        rm \$local_file
    	fi
       
  """

}

// wget -q -O- ${figshare_path}/${params.reference_file}?access_token=${params.figshare_token} | tar -xz




