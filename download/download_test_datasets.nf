


// nextflow run /home/jersal/workspace/cactus/download/download_test_datasets.nf --specie worm

process download_test_datasets {

  container = params.skewer_pigz

  publishDir path: "${launchDir}", mode: "link"

  input:

  output:
    file("*")  
    
  script:
  def figshare_path = "https://ndownloader.figshare.com/files"
  
  """
      local_file=${params.specie}_test_dataset.tar.gz
      wget ${figshare_path}/${params.test_dataset_file}?access_token=${params.figshare_token} -O \$local_file
      local_md5sum=\$(md5sum \$local_file | awk '{print \$1}')
      
      if [[ "\$local_md5sum" == "$params.test_dataset_md5sum" ]] ; then
        tar -xvf \$local_file
        rm \$local_file
      else
        echo "md5sum is wrong for file \$local_file"
        rm \$local_file
    	fi
       
  """

}

// wget -q -O- ${figshare_path}/${params.test_dataset_file}?access_token=${params.figshare_token} | tar -xz






// process download_test_datasets {
//   tag "${id}"
// 
//   container = params.skewer_pigz
// 
//   when: params.download_test_datasets
// 
//   publishDir path: "${launchDir}", mode: "${pub_mode}"
// 
//   input:
// 
//   output:
//     file("*") into Test_datasets_donwloaded
//     // val('true') into Test_datasets_donwloaded
// 
//   script:
//   """
// 
//     figshare_path="https://ndownloader.figshare.com/files"
//     wget -q -O- ${figshare_path}/${params.test_datasets_file}?access_token=${params.figshare_token} | tar -xz
// 
// 
//   """
//   // test_datasets_donwloaded = true
// }
// // wget ${figshare_path}/${data_file}?access_token=${figshare_token} -O ${data_dir}/${specie}_data.tar.gz"



