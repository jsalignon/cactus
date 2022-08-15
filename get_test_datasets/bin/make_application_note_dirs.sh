
dataset_name=application_note
AN_human=$dataset_name/human
AN_worm=$dataset_name/worm

mkdir -p $AN_human/conf $AN_worm/conf

#### worm
cp -r worm/design $AN_worm
sed -i 's/data.*K_reads_atac/..\/..\/preprocessing\/worm\/fastq\/atac/' application_note/worm/design/atac_fastq.tsv
sed -i 's/data.*K_reads_mrna/..\/..\/preprocessing\/worm\/fastq\/mrna/' application_note/worm/design/mrna_fastq.tsv

# run.config
cat > $AN_worm/conf/run.config <<EOL
params {
  specie            = 'worm'
  use_input_control = false
  save_bed_type     = 'last'
  chip_ontology     = 'all'
  chromatin_state   = 'iHMM.M1K16.worm_L3'
  threshold_type_for_splitting_subsets   = 'fdr' 
  threshold_values_for_splitting_subsets = [ 1.3, 3 ]
  
  memory_picard     = '20G'
  executor {
    queueSize = 50
    \$local { 
      memory  = '100 GB'
    }
  }
}
EOL



#### human
cp -r human/design $AN_human
sed -i 's/data.*K_reads_atac/..\/..\/preprocessing\/human\/fastq\/atac/' application_note/human/design/atac_fastq.tsv
sed -i 's/data.*K_reads_mrna/..\/..\/preprocessing\/human\/fastq\/mrna/' application_note/human/design/mrna_fastq.tsv

# run.config
cat > ${specie}/conf/run.config <<EOL
params {
  specie            = 'human'
  use_input_control = false
  save_bed_type     = 'all'
  chip_ontology     = 'cell_type.fibroblast'
  chromatin_state   = 'ENCFF941SVR'
  threshold_type_for_splitting_subsets   = 'fdr'  
  threshold_values_for_splitting_subsets = [ 1.3, 3 ]
}
EOL
