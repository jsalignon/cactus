
cd /home/jersal/workspace/cactus/figshare

compress_test_dataset () {
  specie=$1 ; ncores=$2
  tar --use-compress-program="pigz -p ${ncores} -k -r" -cvf ${specie}_test.tar.gz -C ../test_datasets/${specie} {conf,design,data}
}

compress_data () {
  specie=$1 ; ncores=$2
  tar --use-compress-program="pigz -p ${ncores} -k -r" -cvf ${specie}_data.tar.gz ../data/${specie}
}

compress_specie () {
  specie=$1 ; ncores=$2
  compress_test_dataset $specie $ncores
  compress_data         $specie $ncores
}


compress_specie worm 40
compress_specie fly 40
compress_specie mouse 40
compress_specie human 40

md5sum *.tar.gz > md5_sums.txt


lftp -u "user,password" ftps.figshare.com 
set ftp:ssl-force true
set ftp:ssl-protect-data true
set ftp:passive-mode true
set ftp:auto-passive-mode true
set ftp:ssl-force on
set ftp:ssl-protect-data on
set ftp:passive-mode on
set ftp:auto-passive-mode on
set ssl:verify-certificate no
# mkdir -p data/cactus
cd data/cactus
# put worm_data.tar.gz
# put worm_test.tar.gz
# put fly_data.tar.gz
# put fly_test.tar.gz
# put mouse_data.tar.gz
# put mouse_test.tar.gz
put human_data.tar.gz
put -c human_data.tar.gz
reput human_data.tar.gz

# put human_test.tar.gz
# put README.txt
# put manifest.txt
cat debug_log.txt 
## => this last commands allows to check if the put command worked. The upload often fails (i.e. " Incomplete upload, skipping..."). In this case, the command should be repeated. It should be written: "started" and then on the line below "finished" with a file_id.
## => commentting files that are already on the server and do not need re-uploading



# resume upload
# https://unix.stackexchange.com/questions/91602/is-there-any-ftp-command-to-resume-upload
# => reput command!








# # my private link:
# https://figshare.com/s/d943e9f8f7b014420b0a
# 
# # downloading the file
# TOKEN="*"
# FILE=36732768
# wget -c https://ndownloader.figshare.com/files/36067514?access_token=$TOKEN -O worm_data.tar.gz


