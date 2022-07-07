
params.ensembl_release = '102'

params.number_of_cores = 8   // how to initialize if missing?
number_of_cores = params.number_of_cores

Species_channel = Channel
	.from( 	
		[
			['worm',  'caenorhabditis_elegans',  'ce11', 'WBcel235', 'toplevel' ],
			['fly',   'drosophila_melanogaster', 'dm6',  'BDGP6.28', 'toplevel' ],
			['mouse', 'mus_musculus',            'mm10', 'GRCm38',   'primary_assembly' ],
		  ['human', 'homo_sapiens',            'hg38', 'GRCh38',   'primary_assembly' ]
		]
	)
	.dump(tag: 'start_channel')
	.multiMap { it ->
		    fasta: it[0, 1, 3, 4]
		blacklist: it[0, 2, 3]
		    orgdb: it[0, 1]
	}
	.set { Start_channel }





// channel_species = 
// 	Channel.from( [ 'worm', 'fly', 'mouse', 'human' ] )
// 	.flatten()
// 	.map{ [ it, file("${it}/genome/annotation/annotation.gff3"), file("${it}/genome/sequence/genome.fa") ].flatten() }
// 	.dump(tag:"species")
// 	.into{ channel_species_1 ; channel_species_2 ; channel_species_3 ; channel_species_4 }


process get_fasta_and_gff {
	tag "${specie}"

	container = params.samtools_bedtools_perl

	publishDir path: "${specie}/genome", mode: 'link', saveAs: {
    if (it.indexOf(".gff3") > 0) "annotation/${it}"
		else if (it.indexOf(".fa") > 0) "sequence/${it}"
    else if (it.indexOf(".txt") > 0) "${it}"
  }

	input:
		set specie, specie_long, genome, assembly from Start_channel.fasta

	output:
		set specie, file('annotation.gff3'), file('genome.fa') into (Genome_and_annotation, Genome_and_annotation_1, Genome_and_annotation_2, Genome_and_annotation_3)

	shell:
	'''
		
	specie='!{specie}'
	specie_long='!{specie_long}'
	genome='!{genome}'
	assembly='!{assembly}'
	release='!{params.ensembl_release}'
	
	URL=ftp://ftp.ensembl.org/pub/release-$release
	
	gff3_file=${specie_long^}.$genome.$release.gff3.gz
	fasta_file=${specie_long^}.$genome.dna_sm.$assembly.fa.gz
	
	wget -O annotation.gff3.gz $URL/gff3/$specie_long/$gff3_file 
	gunzip annotation.gff3.gz 
	
	wget -O genome.fa.gz $URL/fasta/$specie_long/dna/$fasta_file 
	gunzip genome.fa.gz
	
	echo "gff3 file : $gff3_file" > README.txt
	echo "fasta file: $fasta_file" >> README.txt
		
		
	'''
}




process get_blacklisted_regions {
	tag "${specie}"

	container = params.cvbio

	publishDir path: "${specie}/blacklisted_regions", mode: 'link'

	input:
		set specie, specie_code, ncbi_code from Start_channel.blacklist

	output:
		file("${specie_code}_blacklist_Ensembl.bed")

	shell:
	'''
			
		specie_code='!{specie_code}'
		ncbi_code='!{ncbi_code}'
		
		url_blacklist="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists"
		url_mapping="https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master"

		wget -O ${specie_code}_blacklist_NCBI.bed.gz $url_blacklist/${specie_code}-blacklist.v2.bed.gz
		gunzip  ${specie_code}_blacklist_NCBI.bed.gz
		
		ncbi_code_1=${ncbi_code%%.*}
				
		wget -O NCBI_to_Ensembl.txt $url_mapping/${ncbi_code_1}_UCSC2ensembl.txt
		
		cvbio UpdateContigNames -i ${specie_code}_blacklist_NCBI.bed -o ${specie_code}_blacklist_Ensembl.bed -m NCBI_to_Ensembl.txt --comment-chars '#' --columns 0 --skip-missing true
		
	'''
}


// this lines fails with a weird Nextflow bug of variable substitution:
// ncbi_code_1=${ncbi_code//\.*/}
// Script compilation error
// - file : /home/jersal/workspace/cactus/software/get_data/parse_genome_files.nf
// - cause: Unexpected character: '\'' @ line 96, column 4.


// AnnotationDbi::saveDb(org.Ce.eg.db, file = 'org_db_ce.sqlite')
// saveRDS(cur_seq_info, '~/workspace/cactus/data/worm/genome/sequence/cur_seqinfo.rds')
// 
// cur_seq_info = rtracklayer::SeqinfoForUCSCGenome('ce11')
// cur_seq_info@seqnames %<>% gsub('chr', '', .)


process filtering_annotation_file {
  tag "${specie}"

  container = params.samtools_bedtools_perl

  publishDir path: "${specie}/genome/annotation", mode: 'link'

  input:
    set specie, file(gff3), file(fasta) from Genome_and_annotation

  output:
		set specie, file('anno_genes_only.gff3'), file('gff3_without_pseudogenes_and_ncRNA.gff3') into channel_gff3_filtered

  shell:
  '''

		grep -i protein_coding !{gff3} | awk  -F"\t"  '{if($3 == "gene") print $0}' - > anno_genes_only.gff3

		grep -viE "pseudogen|ncRNA" !{gff3} > gff3_without_pseudogenes_and_ncRNA.gff3


  '''

}

// Notes

// grep -i protein_coding ${gff} > protein_coding_mRNA_and_genes.gff3

// anno_without_pseudogenes.gff3 is used for creating a txdb database without non coding transcripts, for the annotation of ATAC seq peaks

// protein_coding_genes.gff3 is used to create the gene metadata table for all analysis


process generating_bowtie2_indexes {
  tag "${specie}"

  container = params.bowtie2_samtools

  publishDir path: "${"${specie}/genome/sequence"}/bowtie2_indexes", mode: 'link'

  input:
    set specie, file(gff3), file(fasta) from Genome_and_annotation_1

  output:
		file("*.bt2")

  shell:
  '''

		bowtie2-build --threads !{number_of_cores} !{fasta} genome


  '''

}


process indexing_genomes {
  tag "${specie}"

  container = params.bowtie2_samtools

  publishDir path: "${"${specie}/genome/sequence"}", mode: 'link'

  input:
    set specie, file(gff3), file(fasta) from Genome_and_annotation_2

  output:
		set specie, file('genome.fa.fai') into Fasta_fai

  shell:
  '''

		samtools faidx !{fasta}


  '''

}

Genome_and_annotation_3
	.join(Fasta_fai)
	.set{ Fasta_fai_gff3 }



process getting_chromosome_length_and_transcriptome {
  tag "${specie}"

	container = params.gffread

  publishDir path: "${"${specie}/genome/sequence"}", mode: 'link'

  input:
    set specie, file(gff3), file(fasta), file(fasta_indexes) from Fasta_fai_gff3

  output:
    file('chromosome_size.txt')
		set specie, file('transcriptome.fa') into Transcriptome_for_building_kallisto_indexes

  script:
  """
      cut -f1-2 ${fasta_indexes} > chromosome_size.txt
      gffread -C ${gff3} -g ${fasta} -w transcriptome.fa
  """

}





process generating_kallisto_transcriptome_indexes {
  tag "${specie}"

  container = params.kallisto

  publishDir path: "${"${specie}/genome/sequence"}", mode: 'link'

  input:
		set specie, file(transcriptome) from Transcriptome_for_building_kallisto_indexes

  output:
    file('transcriptome_kallisto_index')

  script:
  """
      kallisto index --make-unique -i transcriptome_kallisto_index ${transcriptome}
  """

}




process getting_orgdb {
  tag "${specie}"

  container = params.annotationhub

  publishDir path: "${specie}/orgdb", mode: 'link'

  input:
		set specie, specie_long from Start_channel.orgdb

  output:
    file('orgdb.sqlite')

  shell:
  '''

		#!/usr/bin/env Rscript

		library(AnnotationHub)

		specie = '!{specie}'
		specie_long = '!{specie_long}'
		annotationhub_cache = '!{params.annotationhub_cache}'

		specie_initial =  sapply(strsplit(specie_long, '_')[[1]], substr, 1, 1) %>% {paste0(toupper(.[1]), .[2])}

		ah = AnnotationHub(cache = annotationhub_cache)

		orgdb = query(ah, paste0('ncbi/standard/3.14/org.', specie_initial, '.eg.db.sqlite'))[[1]]

		AnnotationDbi::saveDb(orgdb, 'orgdb.sqlite')

  '''

}

// orgdb = AnnotationDbi::loadDb('orgdb.sqlite')

// AnnotationDbi::saveDb(org.Ce.eg.db, file = 'org_db_ce.sqlite')
// saveRDS(cur_seq_info, '~/workspace/cactus/data/worm/genome/sequence/cur_seqinfo.rds')
// 


// ah = AnnotationHub(cache = annotationhub_cache)
// Error: Corrupt Cache: index file
//   See AnnotationHub's TroubleshootingTheCache vignette section on corrupt cache
//   cache: /home/jersal/workspace/cactus/util/annotationhub_cache
//   filename: annotationhub.index.rds

// => this error was solved by creating the cache in the cache folder direclty with the container and the command:
// ah = AnnotationHub(cache = '.')

// https://bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHub-HOWTO.html

// singularity shell bioconductor-annotationhub:3.2.0--r41hdfd78af_0
// R
// library(AnnotationHub)
// ah = AnnotationHub()
// specie_initial = 'Ce'
// orgdb_sqlite_file = paste0('org.', specie_initial, '.eg.db.sqlite')
// mcols(query(ah, 'OrgDb', '2021-10-08'))[1,] %>% as.data.frame
// mcols(query(ah, 'ncbi/standard/3.14/org.Ce.eg.db.sqlite'))[1,] %>% as.data.frame
// orgdb = query(ah, 'ncbi/standard/3.14/org.Ce.eg.db.sqlite')[[1]]

// mcols(query(ah, 'ncbi/standard/3.14/org.Ce.eg.db.sqlite'))[1,] %>% as.data.frame %>% t
//                    AH95964
// title              "org.Ce.eg.db.sqlite"
// dataprovider       "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/"
// species            "Caenorhabditis elegans"
// taxonomyid         6239
// genome             "NCBI genomes"
// description        "NCBI gene ID based annotations about Caenorhabditis elegans"
// coordinate_1_based 1
// maintainer         "Bioconductor Maintainer <maintainer@bioconductor.org>"
// rdatadateadded     "2021-10-08"
// preparerclass      "OrgDbFromPkgsImportPreparer"
// tags               character,3
// rdataclass         "OrgDb"
// rdatapath          "ncbi/standard/3.14/org.Ce.eg.db.sqlite"
// sourceurl          "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/, ftp://ftp.ensembl.org/pub/current_fasta"
// sourcetype         "NCBI/ensembl"

// ah[ah$dataprovider == 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/' & ah$rdataclass == 'OrgDb',]
// mcols(ah[ah$dataprovider == 'Ensembl', ]) %>% .[nrow(.),] %>% as.data.frame %>% t
// 
// mcols(ah[ah$dataprovider == 'Ensembl' & ah$sourcetype == 'FASTA' & , ])
// 
// mcols(ah[ah$dataprovider == 'Ensembl' & ah$rdatadateadded == '2021-10-20' & , ])
// ""

// mcols(ah[ah$dataprovider == 'Ensembl', ])$sourcetype %>% table
//   FASTA     GTF ensembl
//   15163   10156    2962


// Playing around with this package!

// > mcols(ah[ah$dataprovider == 'Ensembl', ]) %>% .[nrow(.),] %>% as.data.frame %>% t
//                    AH100303
// title              "Zosterops_lateralis_melanops.ASM128173v1.ncrna.2bit"
// dataprovider       "Ensembl"
// species            "zosterops lateralis_melanops"
// taxonomyid         1220523
// genome             "ASM128173v1"
// description        "TwoBit ncRNA sequence for zosterops lateralis_melanops"
// coordinate_1_based 1
// maintainer         "Bioconductor Maintainer <maintainer@bioconductor.org>"
// rdatadateadded     "2021-10-20"
// preparerclass      "EnsemblTwoBitPreparer"
// tags               character,5
// rdataclass         "TwoBitFile"
// rdatapath          "ensembl/release-105/fasta/zosterops_lateralis_melanops/ncrna/Zosterops_lateralis_melanops.ASM128173" [truncated]
// sourceurl          "ftp://ftp.ensembl.org/pub/release-105/fasta/zosterops_lateralis_melanops/ncrna/Zosterops_lateralis_" [truncated]
// sourcetype         "FASTA"

// > 
// mcols(ah[ah$dataprovider == 'Ensembl', ]) %>% .[nrow(.),] %>% as.data.frame %>% t
// 
// mcols(query(ah, 'ensembl/release-105'))[1, ] %>% as.data.frame %>% t
// mcols(query(ah, 'ensembl/release-105', 'caenorhabditis elegans'))
//  ] %>% as.data.frame %>% t
// 
// df = mcols(query(ah[ah$species == 'caenorhabditis elegans' & ah$description == 'TwoBit DNA sequence for caenorhabditis elegans'], 'ensembl/release-105', 'dna_sm')) 
// 
// mcols(query(ah[ah$species == 'caenorhabditis elegans'], 'ensembl/release-105', 'dna_sm')) 
// mcols(query(query(ah[ah$species == 'caenorhabditis elegans'], 'ensembl/release-105'), 'dna_sm')) %>% as.data.frame %>% t
// mcols(query(query(query(ah, 'caenorhabditis elegans'), 'ensembl/release-105'), 'dna_sm')) %>% as.data.frame %>% t
// query(query(query(ah, 'caenorhabditis elegans'), 'ensembl/release-105'), 'dna_sm')[[1]]





// process getting_R_annotation_files {
//   tag "${specie}"
// 
//   container = params.bioconductor
// 
//   publishDir path: "${specie}/genome/annotation", mode: 'link'
// 
//   input:
// 		set specie, file(gff3_genes_only), file(gff3_without_pseudogenes_and_ncRNA) from channel_gff3_filtered
// 
//   output:
//     file('*')
// 
//   shell:
//   '''
// 
// 		#!/usr/bin/env Rscript
// 
// 		library(magrittr)
// 		library(AnnotationDbi)
// 
// 		gff3_genes_only = '!{gff3_genes_only}'
// 		gff3_without_pseudogenes_and_ncRNA = '!{gff3_without_pseudogenes_and_ncRNA}'
// 		specie = '!{specie}'
// 		genomedir = '!{genomedir}'
// 
// 		specie_initial = switch(genome, 
// 			WBcel235 = 'ce', 
// 			BDGP6.28 = 'dm', 
// 			GRCm38 = 'mm', 
// 			GRCh38 = 'hs'
// 		)
// 
// 		orgdb = AnnotationDbi::loadDb(paste0(genomedir, '/../org_db/org_db_', specie_initial, '.sqlite'))
// 
// 		# export annotation dataframe
// 		anno_df = rtracklayer::readGFF(gff3_genes_only) 
// 		entrez_id = mapIds(orgdb, keys = anno_df$gene_id,  column = 'ENTREZID', keytype = 'ENSEMBL', multiVals = 'first')
// 		anno_df %<>% dplyr::rename(gene_name = Name)
// 		anno_df %<>% dplyr::mutate(width = end - start)
// 		anno_df$entrez_id = entrez_id
// 		anno_df$chr = as.character(anno_df$seqid)
// 		anno_df %<>% dplyr::select(chr, start, end, width, strand, gene_name, gene_id, entrez_id)
// 		saveRDS(anno_df, 'df_genes_metadata.rds', version = 2)
// 
// 		# export txdb
// 		txdb = GenomicFeatures::makeTxDbFromGFF(gff3_without_pseudogenes_and_ncRNA)
// 		AnnotationDbi::saveDb(txdb, file="txdb_without_pseudogenes_and_ncRNA.sqlite")
// 
// 		# export gene vs transcripts 
// 		df_genes_transcripts = select(txdb, keys = keys(txdb), columns = c('GENEID', 'TXNAME'), keytype = 'GENEID')
// 		saveRDS(df_genes_transcripts, 'df_genes_transcripts.rds', version = 2)
// 
// 		# extract promoters
// 		promoters = GenomicFeatures::promoters(GenomicFeatures::genes(txdb), upstream = 1500, downstream = 500)
// 
// 		# export promoters as dataframe
// 		promoters_df = as.data.frame(promoters, stringsAsFactors = F)
// 		promoters_df$start[promoters_df$start < 0] = 0
// 		promoters_df$end[promoters_df$end < 0] = 0
// 		promoters_df %<>% dplyr::rename(gene_name = gene_id)
// 		promoters_df %<>% dplyr::inner_join(anno_df %>% dplyr::select(gene_id, gene_name, entrez_id), by = 'gene_name')
// 		promoters_df %<>% dplyr::arrange(seqnames, start)
// 		promoters_df %<>% dplyr::rename(chr = seqnames)
// 		saveRDS(promoters_df, file = 'promoters_df.rds', version = 2)
// 
// 		# export promoters as bed file
// 		promoters_df1 = promoters_df
// 		promoters_df1$score = 0
// 		promoters_df1 %<>% dplyr::select(chr, start, end, gene_name, score, strand, gene_id)
// 		export_df_to_bed(promoters_df1, 'promoters.bed')
// 
// 
//   '''
// 
// }

// source('!{cactusdir}/src/R_analysis/general_functions.R')

// cur_seq_info = rtracklayer::SeqinfoForUCSCGenome('ce11')
// cur_seq_info@seqnames %<>% gsub('chr', '', .)



// rtracklayer::export(promoters, 'promoters.bed')

// unfortunately, the "promoters function" extract promoters with coordinate below zero (I found this bug for 3 mitochondrial genes in C elegans). Thus we need to put them at zero.

// promoters_df %<>% dplyr::select(chr:strand, gene_name, gene_id, entrez_id)
// promoters_df %<>% dplyr::mutate(seqnames = as.character(seqnames), strand = as.character(strand))


// Rle(promoters_df$seqnames)
// factor-Rle of length 20191 with 7 runs
//   Lengths:  2908  3520  2683  3285  4994  2789    12
//   Values :     I    II   III    IV     V     X MtDNA



// df_genes_transcripts = select(txdb, keys = keys(txdb), columns = c('GENEID', 'TXNAME'), keytype = "GENEID")
// dim(df_genes_transcripts)
// [1] 33552     2

// a = select(txdb, keys = keys(txdb)[1:5], columns = c('GENEID', 'TXID', 'TXNAME', 'TXSTART', 'TXEND', 'TXSTRAND'), keytype = "GENEID")
//   GENEID  TXID               TXNAME TXSTRAND  TXSTART    TXEND
// 1 2L52.1  5291 transcript:2L52.1a.1        +     1848     4717
// 2 2L52.1  5292 transcript:2L52.1b.1        +     3410     4663
// 3 4R79.2 22496 transcript:4R79.2b.1        - 17480396 17482329
// 4 4R79.2 22497 transcript:4R79.2a.1        - 17480396 17483332

// columns(txdb)
// [1] "CDSCHROM"   "CDSEND"     "CDSID"      "CDSNAME"    "CDSPHASE"
// [6] "CDSSTART"   "CDSSTRAND"  "EXONCHROM"  "EXONEND"    "EXONID"
// [11] "EXONNAME"   "EXONRANK"   "EXONSTART"  "EXONSTRAND" "GENEID"
// [16] "TXCHROM"    "TXEND"      "TXID"       "TXNAME"     "TXSTART"
// [21] "TXSTRAND"   "TXTYPE"

// vec = anno_df$gene_name %>% setNames(., anno_df$gene_id)
// saveRDS(vec, 'gene_id_to_name_vector.rds', version = 2)
// 
// vec = anno_df$entrez_id %>% setNames(., anno_df$gene_id)
// saveRDS(vec, 'gene_ensembl_id_to_entrez_id_vector.rds', version = 2)


// saveDfFromGRangesAsBed(all_promoters, file = 'promoters.bed')
// this line of code was replaced by:
// rtracklayer::export(promoters, 'promoters.bed')

 
 
// a = GenomicFeatures::promoters(GenomicFeatures::genes(txdb), upstream = 1500, downstream = 500)
// rtracklayer::export(a, 'test.bed')
// rtracklayer::export(a, 'test.bed', format = 'bed')
// 
// txdb = 
// a = annotatePeak(gr, c(-200, 200), level = 'gene')



// txdb = GenomicFeatures::makeTxDbFromGFF("/home/jersal/lluis/atac/bin/igenomes/references/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gff3")
// 
// txdb = GenomicFeatures::makeTxDbFromGFF("genes.gff3")
// 
// anno_df1 = rtracklayer::readGFF(anno_file) 
// gr = GenomicRanges::makeGRangesFromDataFrame(anno_df1)



// => in fact I don't need a txdb object
// GenomicFeatures::promoters(a, upstream = 500, downstream = 500)
// GenomicFeatures::promoters(a, upstream = 1, downstream = 1)


// library(biomaRt)

// specie_initial = paste(sapply(strsplit(tolower(specie), '_')[[1]], function(x) substr(x, 1,1)), collapse = '')


// dataset = switch(specie, 
// 	Caenorhabditis_elegans = 'celegans_gene_ensembl', 
// 	Drosophila_melanogaster = 'dmelanogaster_gene_ensembl', 
// 	Mus_musculus = 'mmusculus_gene_ensembl', 
// 	Homo_sapiens = 'hsapiens_gene_ensembl'
// )
// 
// mart = useMart("ensembl", dataset = dataset)
// 
// ensembl_gene_id = gff3_df$gene_id
// df_biomart = getBM( filters = 'ensembl_gene_id', attributes = c('ensembl_gene_id', 'entrezgene_id'), values = ensembl_gene_id, mart = mart )
// 
// gff3_df$entrez_gene_id = df_biomart$entrezgene_id[match(gff3_df$gene_id, df_biomart$ensembl_gene_id)]



// mart = useMart("ensembl")
// a = listDatasets(mart)
// a[grepl('elegans|hsapiens|musculus|melanogaster', a$dataset), ]
// 									 dataset                              description    version
// 48       hsapiens_gene_ensembl                 Human genes (GRCh38.p13) GRCh38.p13
// 53  dmelanogaster_gene_ensembl Drosophila melanogaster genes (BDGP6.28)   BDGP6.28
// 104      celegans_gene_ensembl  Caenorhabditis elegans genes (WBcel235)   WBcel235
// 130     mmusculus_gene_ensembl                  Mouse genes (GRCm38.p6)  GRCm38.p6

// ## for C elegans
// dim(gff3_df) # 31293
// gff3_df1 = gff3_df[, c('gene_id', 'entrez_gene_id')] %>% .[!duplicated(.), ]
// dim(gff3_df1) # 20447
// length(which(is.na(gff3_df1$entrez_gene_id))) # => 563 id were lost out of 20K. Not great but still reasonnable.



// on completion
 workflow.onComplete
 {
  println ""
  println "Workflow completed on: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'Succeeded' : 'Failed' }"
  println "Workflow Duration: $workflow.duration"
  println ""
 }
 
