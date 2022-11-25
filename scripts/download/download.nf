

process download_test_datasets {
  tag params.species

  label "skewer_pigz"

  publishDir path: "${launchDir}", mode: params.pub_mode

  when: params.test_datasets

  input:

  output:
    set file('data/'), file('parameters/'), file('design/')
    
  script:
  def figshare_path = "https://ndownloader.figshare.com/files"
  def local_file    = "${params.species}_test_dataset.tar.gz"
  
  """
  wget ${figshare_path}/${params.test_dataset_file} -O ${local_file}
  local_md5sum=\$(md5sum ${local_file} | awk '{print \$1}')
  
  if [[ "\$local_md5sum" == "${params.test_dataset_md5sum}" ]] ; then
    pigz --decompress --keep --recursive --stdout --processes ${params.threads} ${local_file} | tar -xvf -
  else
    echo "md5sum is wrong for file ${local_file}"
  fi
  """

}

// tar --use-compress-program="pigz -p ${params.threads} -k -r" -xvf ${local_file} => the pigz option is not avaialble in busybox's tar
// pigz -d -k -r -c -p 3 $local_file | tar -xvf - // same command with abbreviations

process download_references {
  tag params.species

  label "skewer_pigz"

  publishDir path: "${params.references_dir}", mode: params.pub_mode

  when: params.references

  input:

  output:
    file(params.species)
  
  script:
  def figshare_path = "https://ndownloader.figshare.com/files"
  def local_file    = "${params.species}_references.tar.gz"
  
  """
  wget ${figshare_path}/${params.references_file} -O ${local_file}
  local_md5sum=\$(md5sum ${local_file} | awk '{print \$1}')
  
  if [[ "\$local_md5sum" == "${params.references_md5sum}" ]] ; then
    pigz --decompress --keep --recursive --stdout --processes ${params.threads} ${local_file} | tar -xvf -
  else
    echo "md5sum is wrong for file ${local_file}"
  fi
  """

}


