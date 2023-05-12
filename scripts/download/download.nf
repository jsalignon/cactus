


Files_to_download = Channel
  .from(  
    [  
      ['test_dataset',  params.test_datasets, params.test_dataset_file, params.test_dataset_md5sum],
      ['references',    params.references,    params.references_file,   params.references_md5sum]
    ]
  )
  .filter{ it[1] == true}
  .map{ it[0, 2, 3] }
  // .view()
  .set{ Files_to_download_1 }

process download_files {
  tag "${params.species} ${params.figshare_version} ${file_type}"

  label "gnu_wget"

  input:
    set val(file_type), figshare_file, val(true_md5sum) from Files_to_download_1

  output:
    set val(file_type), "${params.species}_${file_type}_checked.tar.gz" into Files_to_extract
    
  script:
    def figshare_path = "https://ndownloader.figshare.com/files"
    def local_file    = "${params.species}_${file_type}.tar.gz"
    def local_file_1  = "${params.species}_${file_type}_checked.tar.gz"
    
    """
    wget ${figshare_path}/${figshare_file} -O ${local_file}
    local_md5sum=\$(md5sum ${local_file} | awk '{print \$1}')
    
    if [[ "\$local_md5sum" == "${true_md5sum}" ]] ; then
      mv ${local_file} ${local_file_1}
    else
      echo "md5sum is wrong for file ${local_file}"
    fi
    """

}


process extract_files {
  tag "${params.species} ${params.figshare_version} ${file_type}"

  label "skewer_pigz"

  // publishDir path: "${launchDir}", mode: params.pub_mode, enabled: file_type == "test_dataset"
  // publishDir path: "${params.references_dir}", mode: params.pub_mode, enabled: file_type == "references"
  
  publishDir path: "/", mode: params.pub_mode, saveAs: {
           if (file_type == 'test_dataset')     "${launchDir}/${it}"
      else if (file_type == 'references') "${params.references_dir}/${it}"
    }

  input:
    set val(file_type), downloaded_file from Files_to_extract

  output:
      set file('data/'), file('parameters/'), file('design/'), 
          file(params.species) optional true

  script:
    """
    pigz --decompress --keep --recursive --stdout \
      --processes ${params.threads} ${downloaded_file} | tar -xvf -
    """

}

// tar --use-compress-program="pigz -p ${params.threads} -k -r" -xvf ${local_file} => the pigz option is not avaialble in busybox's tar
// pigz -d -k -r -c -p 3 $local_file | tar -xvf - // same command with abbreviations

