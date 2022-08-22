
cd /home/jersal/workspace/cactus/figshare

compress_specie <-
tar --use-compress-program="pigz -p 40 -k -r" -cvf worm_data.tar.gz ../data/worm
tar --use-compress-program="pigz -p 40 -k -r" -cvf worm_test.tar.gz -C ../test_datasets/worm {conf,design,data}
tar --use-compress-program="pigz -p 40 -k -r" -cvf worm_test.tar.gz -C ../test_datasets/worm {conf,design}




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
mkdir -p data/cactus
cd data/cactus
put worm_data.tar.gz
put worm_test.tar.gz
put README.txt
put manifest.txt
cat debug_log.txt 
# => this last commands allows to check if the put command worked. The upload often fail (i.e. " Incomplete upload, skipping...") Need to repeat it sometimes. It should be written: "started" and then on the line below "finished"














# # my private link:
# https://figshare.com/s/d943e9f8f7b014420b0a
# 
# # downloading the file
# TOKEN="*"
# FILE=36732768
# wget -c https://ndownloader.figshare.com/files/36067514?access_token=$TOKEN -O worm_data.tar.gz


