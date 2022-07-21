

params.number_of_cores = 8   // how to initialize if missing?
number_of_cores = params.number_of_cores

params.ensembl_release  = '102'
ncbi_orgdb_version      = '3.14'
homer_genomes_version   = '6.4'
homer_organisms_version = '6.3'
homer_promoters_version = '5.5'

params.homer_odd_score_threshold = 0.5

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
			  	homer: it[0, 2]
				   pwms: it[0, 1]
      blacklist: it[0, 2, 3]
	  			orgdb: it[0, 1]
	}
	.set { Start_channel }





process get_homer_data {
	tag "${specie}"

	container = params.samtools_bedtools_perl

	publishDir path: "${specie}/homer_data", mode: 'link'
	input:
		set specie, specie_code from Start_channel.homer

	output:
		file('*')

	shell:
	'''
		
		specie='!{specie}'
		specie_code='!{specie_code}'
		homer_genomes_version='!{homer_genomes_version}'
		homer_organisms_version='!{homer_organisms_version}'
		homer_promoters_version='!{homer_promoters_version}'
		
		URL=http://homer.ucsd.edu/homer/data
		
		wget -nc $URL/genomes/${specie_code}.v${homer_genomes_version}.zip
		wget -nc $URL/organisms/${specie}.v${homer_organisms_version}.zip
		wget -nc $URL/promoters/${specie}.v${homer_promoters_version}.zip
		
		for z in *.zip; do unzip "$z"; done
		mv data/* .
		rm -r data *.zip
		 
	'''
}



process get_pwms {

	container = params.r_basic

	publishDir path: "util/motifs_PWMS", mode: 'link'

	output:
		file('dt_cisbp_encode.rds') into cisbp_motifs_all_species
		file('*')
	
	shell:
	'''
		#!/usr/bin/env Rscript
		
		homer_odd_score_threshold = !{params.homer_odd_score_threshold}
		
		url="http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/Bulk_downloads/EntireDataset"
		f1 = 'PWMs.zip'               ; download.file(paste0(url, '/', f1), f1)
		f1 = 'TF_Information.txt.zip' ; download.file(paste0(url, '/', f1), f1)
		
		unzip('TF_Information.txt.zip')
		unzip('PWMs.zip')
		
		system("sed -i 's/#/|/g' TF_Information.txt")
		
		library(data.table)
		library(magrittr)
		library(purrr)
	
		dt_all = fread('TF_Information.txt')
		
		encode_species = c('Caenorhabditis_elegans', 'Drosophila_melanogaster', 'Mus_musculus', 'Homo_sapiens')
		dt_cisbp = dt_all[TF_Species %in% encode_species]
		dt_cisbp = dt_cisbp[Motif_ID != '.']
		dt_cisbp = dt_cisbp[!duplicated(TF_ID)]

		# removing emtpy motifs
		nrow(dt_cisbp) # 2936
		nrows = sapply(dt_cisbp$Motif_ID, function(x) nrow(read.table(paste0('pwms/', x, '.txt'), sep = '\t')))
		dt_cisbp = dt_cisbp[nrows != 1]
		nrow(dt_cisbp) # 2779

		# loading motifs and determining odd scores, thresholds and consensus sequences
		dt_cisbp$motif = lapply(dt_cisbp$Motif_ID, function(x) read.table(paste0('pwms/', x,'.txt'), sep = '\t', header = T, stringsAsFactors = F)[, -1])
		dt_cisbp[, log2_odd_score := sapply(motif, function(motif1) sum(apply(motif1, 1, function(cur_row) log2(max(cur_row)/0.25))))]
		dt_cisbp[, homer_threshold := log2_odd_score * homer_odd_score_threshold]
		dt_cisbp[, consensus_sequence := sapply(motif, function(motif1) paste0(colnames(motif1)[apply(motif1, 1, function(x) which.max(x))], collapse = ''))]

		saveRDS(dt_cisbp, 'dt_cisbp_encode.rds')
	
		
	'''
	
}

Start_channel.pwms
	.combine(cisbp_motifs_all_species)
	.set{ cisbp_motifs_all_species_1 }


process split_pwms {

	container = params.r_basic

	publishDir path: "${specie}/motifs_PWMs", mode: 'link'

	input:
		set specie, specie_long, file(dt_cisbp_encode_rds) from cisbp_motifs_all_species_1

	output:
		file('homer_motifs.txt')
	
	shell:
	'''
		#!/usr/bin/env Rscript
		
		library(data.table)
		
		dt_cisbp_encode = readRDS('!{dt_cisbp_encode_rds}')
		specie_long = '!{specie_long}'
		
		specie_long_1 = paste0(toupper(substr(specie_long, 1, 1)), substr(specie_long, 2, nchar(specie_long)))
		dt_cisbp = dt_cisbp_encode[TF_Species == specie_long_1]
		
		
		sink('homer_motifs.txt')
		
			apply(dt_cisbp, 1, function(x) {
					cat(paste0('>', x$consensus_sequence, '\t', x$TF_Name, '\t', x$homer_threshold, '\n'))
					cat(paste0(apply(x$motif, 1, paste0, collapse = '\t'), '\n'))
				})

		sink()
		
	'''
	
}




// bioconductor-universalmotif:1.8.3--r40h399db7b_0

// dt_encode = dt_all[TF_Species %in% encode_species]
// dt_encode = dt_encode[Motif_ID != '.']
// dt_encode1 = dt_encode[!duplicated(TF_ID)]
// > dt_encode$TF_Species %>% table
// .
//  Caenorhabditis_elegans Drosophila_melanogaster            Homo_sapiens
//                     402                    1301                    5422
//            Mus_musculus
//                    1503
// >
// > dt_encode1$TF_Species %>% table
// .
//  Caenorhabditis_elegans Drosophila_melanogaster            Homo_sapiens
//                     373                     425                    1200
//            Mus_musculus
//                     938

// However, there are still multiple entries per TF
// dt_encode[TF_Species == 'Caenorhabditis_elegans']$TF_Name %>% table %>% sort %>% rev %>% head(10)
  // hlh-1   skn-1  snpc-4   pha-4   efl-1   tra-1 nhr-182   mab-3  lin-39   eor-1
  //     7       4       3       3       3       2       2       2       2       2
// dt_encode[TF_Species == 'Caenorhabditis_elegans' & TF_Name == 'skn-1']
// dt_encode[TF_Species == 'Caenorhabditis_elegans' & TF_Name == 'hlh-1']
// dt_encode[TF_Species == 'Caenorhabditis_elegans' & TF_Name == 'snpc-4']
// dt_encode[TF_Species == 'Caenorhabditis_elegans' & TF_Name == 'pha-4']

// I think I understand what happened: the entries that are duplicated are all from the similarity regression method. And they all have the same score. So instead of picking one arbitrarily, CISBP kept all of these entries with the same score. I could keep only the CHIP.
// Conclusions: these are probably exactly the same motifs (or very very similar). So we can just keep the first of them!


// # Identifying and removing empty motifs files
// file_is_correct = map_lgl(paste0('pwms/', dt_encode$Motif_ID, '.txt'), function(f1) {ifelse(nrow(fread(f1)) == 1, F, T)})
// file_is_correct %>% table # => all TRUE
// TRUE
// 8628

// https://github.com/TedCCLeung/PhMotif/blob/9ace03a196e02375b3024af3506b9f58e98f1a1d/data-raw/processCISBPfiles.R


// nice explanations about the cisbp files here:
// https://github.com/Danko-Lab/rtfbs_db/blob/a29f73b24e8984aad8bcd1227b06446ebef44b08/rtfbsdb/man/CisBP.getTFinformation.Rd
// Three TF information files in CisBP dataset.\cr\cr
// 1: TF_Information.txt : (direct motifs) or (no direct but inferred motifs with 90\%)\cr
// 2: TF_Information_all_motifs.txt: (direct motifs) and (inferred motifs above the threshold)\cr
// 3: F_Information_all_motifs_plus.txt: All motifs\cr

// 'TF_Information_all_motifs.txt' is a superset of 'TF_Information.txt'.  It also includes any motif that can be inferred for a given TF, given the TF family-specific threshold.  For example, if a TF has a directly determined motif, and two TFs with motifs with 90% and 80% identical DBDs to it, TF_Information.txt will only include the direct motif, but 
// TF_Information_all_motifs.txt will include all three motifs.  Likewise, if a TF does not have a direct motif, but has two TFs with 90% and 80% identical DBDs to it, TF_Information.txt will only include the motif from the 90% indentical TF, while TF_Information_all_motifs.txt would include both.\cr\cr



// f1 = 'TF_Information_all_motifs.txt.zip' ; download.file(paste0(url, '/', f1), f1)
// unzip('TF_Information_all_motifs.txt.zip')
// system("sed -i 's/#/|/g' TF_Information_all_motifs.txt")
// dt_all = fread('TF_Information_all_motifs.txt')







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





process get_chip_seq {
	tag "${specie}"

	container = params.encodeexplorer

	// publishDir path: "${specie}/blacklisted_regions", mode: 'link'

	input:
		// set specie, specie_code, ncbi_code from Start_channel.blacklist

	output:
		file("${specie_code}_blacklist_Ensembl.bed")

	shell:
	'''

		#!/usr/bin/env Rscript
		
		annotationhub_cache = '!{params.annotationhub_cache}'

		library(ENCODExplorer)
		library(AnnotationHub)
		library(magrittr)
	
		ah = AnnotationHub(cache = annotationhub_cache)
		query(ah, "ENCODExplorerData")


		https://www.encodeproject.org/search/?type=File&file_format=bed&status=released&output_type=optimal+IDR+thresholded+peaks&assay_title=TF+ChIP-seq
		
		url_encode = 'https://www.encodeproject.org/' 
		url_search = paste0(url_encode, 'search/?' 
		url_append = '&frame=object&format=json&limit=all'

		my_query = 'type=File&file_format=bed&status=released&output_type=optimal+IDR+thresholded+peaks&assay_title=TF+ChIP-seq&assembly=ce11'
		my_url = paste0(url_encode, my_query, url_append)
		RCurl::url.exists(my_url)
		res1 = jsonlite::fromJSON(my_url)
		df = res1[['@graph']]
		df$href[1:2]
		df[, c('href', 'md5sum')]
		
		fileName = strsplit(x = file_url, split = "@@download/", fixed = TRUE)[[1]][2]
		fileName <- paste0(dir,"/", fileName, sep="")
		
		# Identify the target URL.
		href <- as.character(file_url)
		
		# Calculate the md5 of the file if it already exists.
		md5sum_file <- tools::md5sum(fileName)
		md5sum_encode <- as.character(file_md5)



	'''
}








library(data.table)

url_encode = 'https://www.encodeproject.org/' 
url_search = paste0(url_encode, 'search/?')
url_append = '&frame=object&format=json&limit=all'

my_query = 'type=File&file_format=bed&status=released&output_type=optimal+IDR+thresholded+peaks&assay_title=TF+ChIP-seq&assembly=ce11'
my_url = paste0(url_search, my_query, url_append)
RCurl::url.exists(my_url)
res1 = jsonlite::fromJSON(my_url)
df = res1[['@graph']]

dt = data.table(df[, c('href', 'md5sum')])
dt = dt[1:5]
dt[, file_name := gsub('.*download/', '', href)]

sapply(1:nrow(dt), function(c1) download.file(url = paste0(url_encode, dt$href[c1]), quiet = T, destfile = dt$file_name[c1], method = 'curl', extra = '-L' ))

dt[, md5sum_dl := tools::md5sum(file_name)]
if(any(dt$md5sum != dt$md5sum_dl)) stop('not all md5 sums are equal')





### FILE
url_append_test = '&frame=embedded&format=json&limit=all'
my_query = 'type=File&accession=ENCFF209GZO'
my_url = paste0(url_search, my_query, url_append_test)
RCurl::url.exists(my_url)
res1 = jsonlite::fromJSON(my_url)
df = res1[['@graph']]
df1 = df[1,]


my_query = 'type=File&file_format=bed&status=released&output_type=optimal+IDR+thresholded+peaks&assay_title=TF+ChIP-seq&assembly=ce11'
my_url = paste0(url_search, my_query, url_append_test)
RCurl::url.exists(my_url)
res1 = jsonlite::fromJSON(my_url)
df = res1[['@graph']]
df1 = df[1,]

## => the @id present in all category indicates the path:
# https://www.encodeproject.org/{cur_id}/?format=json
# i.e. @id=	"/organisms/celegans/"
# path: https://www.encodeproject.org/genes/176069/?format=json

# https://www.encodeproject.org/genes/176069/?format=json
# https://www.encodeproject.org/files/ENCDO503XTK/?format=json
# https://www.encodeproject.org/files/ENCSR342TEL/?format=json
# https://www.encodeproject.org/files/ENCFF209GZO/?format=json
df1[, c('accession', 'target', 'simple_biosample_summary')]
dt1 = data.table(accession = df1$accession, target_name = df1$target$label, target_id = df1$target$genes, classification = df1$biosample_ontology$classification, assembly = df1$assembly, biosample_summary = df1$simple_biosample_summary)

dt1 = data.table(file = df$accession, target_name = df$target$label, target_id = gsub(df$target$genes, pattern = '/.*/(.*)/', replacement = '\\1'), classification = df$biosample_ontology$classification, assembly = df$assembly, biosample_summary = df$simple_biosample_summary, donor = gsub(df$donor, pattern = '/.*/(.*)/', replacement = '\\1'))

# combining with the donor information to have more information on the sample
my_query = 'type=Donor&organism.scientific_name=Caenorhabditis+elegans'
my_url = paste0(url_search, my_query, url_append_test)
RCurl::url.exists(my_url)
res2 = jsonlite::fromJSON(my_url)
df2 = res2[['@graph']]

dt2 = data.table(donor = df2$accession, strain_name = df2$strain_name, genotype = df2$genotype, description = purrr::map_chr(df2$genetic_modifications, ~ifelse(nrow(.x) > 0, ifelse('description' %in% colnames(.x), .x$description, 'NA'), 'NA')))

dt3 = dt1[dt2, , on = 'donor']

dt3[1,] %>% t




my_query = 'type=Gene&organism.scientific_name=Caenorhabditis+elegans'
my_url = paste0(url_search, my_query, url_append_test)
RCurl::url.exists(my_url)
res1 = jsonlite::fromJSON(my_url)
df = res1[['@graph']]
df1 = df[1,]

scientific_name	"Caenorhabditis elegans"

get_encode_df <- function(my_query){
	my_url = paste0(url_search, my_query, url_append)
	if(!RCurl::url.exists(my_url)) stop('url doesn\'t exist')
	res = jsonlite::fromJSON(my_url)
	df = res[['@graph']]
	return(df)
}

// Here are example of json links associated to a given bed file 

// of a worm experiment:
// https://www.encodeproject.org/files/ENCFF667MVT/  => bigBed narrowPeak
// https://www.encodeproject.org/experiments/ENCSR956TMJ/?format=json
// https://www.encodeproject.org/biosamples/ENCBS172XOM/?format=json
// https://www.encodeproject.org/worm-donors/ENCDO164EGI/?format=json
// https://www.encodeproject.org/genetic-modifications/ENCGM416HQF/?format=json
// https://www.encodeproject.org/targets/npax-4-celegans/?format=json
// https://www.encodeproject.org/genes/182466/?format=json

// https://www.encodeproject.org/help/data-organization/ 
// => we don't need the Biosample level: replicates (biosamples) are merged in the donor levels => note 2: we do need the biosample level, as it is there that the ontologies are stored (i.e. cell line, tissue, cell type)
// cell line: biosample_ontology.term_name
// cell type: biosample_ontology.cell_slims
//    tissue: biosample_ontology.organ_slims

// of a human experiment:
// https://www.encodeproject.org/files/ENCFF739GZC/?format=json
// https://www.encodeproject.org/experiments/ENCSR468DVP/?format=json
// https://www.encodeproject.org/human-donors/ENCDO000AAD/?format=json
// https://www.encodeproject.org/biosample-types/cell_line_EFO_0002067/?format=json

// https://www.encodeproject.org/biosample-types/tissue_UBERON_0000160/?format=json


//// Test in worm

df_chip_files = get_encode_df('type=File&file_format=bed&status=released&output_type=optimal+IDR+thresholded+peaks&assay_title=TF+ChIP-seq&assembly=ce11')

df_experiments = get_encode_df('type=Experiment&assembly=ce11')
df_biosamples = get_encode_df('type=Biosample&organism.scientific_name=Caenorhabditis+elegans')
df_donor = get_encode_df('type=Donor&organism.scientific_name=Caenorhabditis+elegans')
df_genetic_modifications = get_encode_df('type=GeneticModification&modified_site_by_target_id.organism=/organisms/celegans/')
df_targets = get_encode_df('type=Target&organism.scientific_name=Caenorhabditis+elegans')
df_genes = get_encode_df('type=Gene&organism.scientific_name=Caenorhabditis+elegans')

get_id <- function(cur_id) exp_dataset <- gsub(cur_id, pattern="/.*/(.*)/", replacement = "\\1")
s
dt_experiment = data.table(df_experiments[, c('accession', 'target', 'life_stage_age', 'biosample_summary', 'description', 'simple_biosample_summary')])

dt_experiment = data.table(
	accession = 
	)

dt_donor = data.table(df_donor[, c('accession', 'genotype', 'strain_name')])

https://www.encodeproject.org/search/?type=File&biosample_ontology.cell_slims=leukocyte&biosample_ontology.term_name=K562&biosample_ontology.organ_slims=blood


//// Test in human

pnrow <- function(x) print(nrow(x))

df_chip_files = get_encode_df('type=File&file_format=bed&status=released&output_type=optimal+IDR+thresholded+peaks&assay_title=TF+ChIP-seq&assembly=GRCh38') %T>% pnrow # 1554
df_experiments = get_encode_df('type=Experiment&assembly=GRCh38') %T>% pnrow # 15362
df_biosamples = get_encode_df('type=Biosample&organism.scientific_name=Homo+sapiens') %T>% pnrow # 15828
// df_donor = get_encode_df('type=Donor&organism.scientific_name=Caenorhabditis+elegans')
// df_genetic_modifications = get_encode_df('type=GeneticModification&modified_site_by_target_id.organism=/organisms/celegans/')
// df_targets = get_encode_df('type=Target&organism.scientific_name=Caenorhabditis+elegans')
// df_genes = get_encode_df('type=Gene&organism.scientific_name=Caenorhabditis+elegans')



/// let's first get all the data and annotation that we need in one process and do the detailed annotation in annother process




library(data.table)
library(magrittr)
library(purrr)

url_encode = 'https://www.encodeproject.org/' 
url_search = paste0(url_encode, 'search/?')
url_append = '&frame=object&format=json&limit=all'

get_encode_df <- function(my_query){
	my_url = paste0(url_search, my_query, url_append)
	if(!RCurl::url.exists(my_url)) stop('url doesn\'t exist')
	res = jsonlite::fromJSON(my_url)
	df = res[['@graph']]
	return(df)
}

pnrow <- function(x) print(nrow(x))


df_chip_files = get_encode_df('type=File&file_format=bed&output_type=optimal+IDR+thresholded+peaks&assembly=GRCh38&assembly=ce11&assembly=dm6&assembly=mm10&assay_title=TF+ChIP-seq&status=released') %T>% pnrow # 2716
df_experiments = get_encode_df('type=Experiment&assay_title=TF+ChIP-seq&status=released') %T>% pnrow # 4425
df_biosample_types = get_encode_df('type=BiosampleType') %T>% pnrow # 936
df_targets = get_encode_df('type=Target&investigated_as=transcription+factor') %T>% pnrow # 6492
df_targets$genes %<>% unlist
df_targets = df_targets[!map_lgl(df_targets$genes, is.null), ] %T>% pnrow # 6490
df_targets = df_targets[!map_lgl(df_targets$genes, is.na), ] %T>% pnrow # 6490
df_genes = get_encode_df('type=Gene') %T>% pnrow # 179933

map_int(df_targets$genes, length) %>% table
//    1    2   23
// 6464   25    1




dt0 = data.table(df_chip_files[, c('accession', 'dataset', 'href')])
dt = dt0[, 1:2]

dt_experiment = data.table(df_experiments[, c('@id', 'biosample_summary', 'biosample_ontology', 'life_stage_age', 'target')])

dt_targets = data.table(df_targets[, c('@id', 'genes')])
dt_targets[, genes := unlist(genes)]
dt_genes = data.table(df_genes[, c('@id', 'geneid', 'symbol')])

dt_genes[dt_targets, , on = c('@id' = 'genes')]

dt_genes[dt_targets, , on = c('genes' = '@id')]

dt_targets[dt_genes, ]
unlist(df_genes$genes) %>% str

// characterization:
df_chip_files$assembly %>% table
  // ce11    dm6 GRCh38   mm10
  //  474    532   1554    156


 df_biosamples[1,] %>% t

// df_biosamples = data.table(df_experiments[, c('accession', 'genotype', 'life_stage_age', 'biosample_summary', 'description', 'simple_biosample_summary')])


df1[, c('title', 'dataset')]

df1$title # [1] "ENCFF209GZO"
df1$dataset # [1] "ENCFF209GZO"



### FILE
url_append_test = '&frame=embedded&format=json&limit=all'
my_query = 'type=File&accession=ENCFF209GZO'
my_url = paste0(url_search, my_query, url_append_test)
RCurl::url.exists(my_url)
res1 = jsonlite::fromJSON(my_url)
df = res1[['@graph']]
df1 = df[1,]


my_query = 'type=File&file_format=bed&status=released&output_type=optimal+IDR+thresholded+peaks&assay_title=TF+ChIP-seq&assembly=ce11'
my_url = paste0(url_search, my_query, url_append_test)
RCurl::url.exists(my_url)
res1 = jsonlite::fromJSON(my_url)
df = res1[['@graph']]
df1 = df[1,]


# https://www.encodeproject.org/files/ENCDO503XTK/?format=json
# https://www.encodeproject.org/files/ENCSR342TEL/?format=json
# https://www.encodeproject.org/files/ENCFF209GZO/?format=json
df1[, c('accession', 'target', 'simple_biosample_summary')]
dt1 = data.table(accession = df1$accession, target_name = df1$target$label, target_id = df1$target$genes, classification = df1$biosample_ontology$classification, assembly = df1$assembly, biosample_summary = df1$simple_biosample_summary)

df1[, c('accession', 'target', 'simple_biosample_summary')]
dt1 = data.table(accession = df$accession, target_name = df$target$label, target_id = df$target$genes, classification = df$biosample_ontology$classification, assembly = df$assembly, biosample_summary = df$simple_biosample_summary)

# => need to combine with the donor information to have the full details on the strain and the notes and so on
https://www.encodeproject.org/worm-donors/ENCDO503XTK/?format=json

# 


https://www.encodeproject.org/search/?type=File&accession=ENCFF209GZO
https://www.encodeproject.org/search/?type=Experiment&accession=ENCSR342TEL
https://www.encodeproject.org/search/?type=Biosample&accession=ENCBS802NRN
https://www.encodeproject.org/search/?type=Donor&accession=ENCDO503XTK
https://www.encodeproject.org/search/?type=Analysis&accession=ENCAN948JEG



### FILE
my_query = 'type=File&accession=ENCFF209GZO'
my_url = paste0(url_search, my_query, url_append)
RCurl::url.exists(my_url)
res1 = jsonlite::fromJSON(my_url)
df = res1[['@graph']]
df1 = df[1,]

df1[, c('accession', 'target', 'simple_biosample_summary')]


## BIOSAMPLE (NOT NEEDED)
my_query = 'type=Biosample&accession=ENCBS802NRN'
my_url = paste0(url_search, my_query, url_append)
RCurl::url.exists(my_url)
res2 = jsonlite::fromJSON(my_url)
df = res2[['@graph']]

## DONOR
my_query = 'type=Donor&accession=ENCDO503XTK'
my_url = paste0(url_search, my_query, url_append)
RCurl::url.exists(my_url)
res2 = jsonlite::fromJSON(my_url)
df = res2[['@graph']]
df1 = df[1,]
df1[, c('accession', 'genotype', 'strain_name', 'notes')]

df1$genetic_modifications
[1] "/genetic-modifications/ENCGM568HEQ/"


## GeneticModification
my_query = 'type=GeneticModification&accession=ENCGM568HEQ'
my_url = paste0(url_search, my_query, url_append)
RCurl::url.exists(my_url)
res2 = jsonlite::fromJSON(my_url)
df = res2[['@graph']]
df1 = df[1,]




my_query = 'dataset=/experiments/ENCSR342TEL/'
my_url = paste0(url_search, my_query, url_append)
RCurl::url.exists(my_url)

res2 = jsonlite::fromJSON(my_url)
df = res2[['@graph']]
df1 = df[1,]

df1[, c('accession', 'target', 'description', 'biosample_summary', 'simple_biosample_summary')]

Biosample {ENCBS590PQG|/biosamples/ENCBS590PQG/}



my_query = 'biosample_ontology.term_id=UBERON_0000468'
my_url = paste0(url_search, my_query, url_append)
RCurl::url.exists(my_url)


my_query = 'type=Biosample&biosample_ontology.term_id=UBERON_0000468'

my_query = 'type=Biosample&term_id=UBERON_0000468'
my_url = paste0(url_search, my_query, url_append)
RCurl::url.exists(my_url)



res2 = jsonlite::fromJSON(my_url)
df = res2[['@graph']]
df



my_query = 'type=Experiment&accession=ENCSR342TEL'
my_url = paste0(url_search, my_query, url_append)
RCurl::url.exists(my_url)

res2 = jsonlite::fromJSON(my_url)
df = res2[['@graph']]

https://www.encodeproject.org/search/?type=Biosample&biosample_ontology=UBERON_0000468



my_query = 'type=Biosample&term_id=UBERON_0000468/'
my_url = paste0(url_search, my_query, url_append)
RCurl::url.exists(my_url)


https://www.encodeproject.org/experiments/ENCSR342TEL/?format=json

type=Replicate&experiment.accession=ENCSR000AKS&format=json&frame=embedded




















https://www.encodeproject.org/search/?type=Experiment&control_type%21=%2A&status=released&perturbed=false&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=TF+ChIP-seq&files.file_type=bed+idr_ranked_peak


# Step 2 : Handling unavailable files via ENCODE rest-api

    filter <- "/&frame=object&format=json&limit=all"
    
    #Simple rest-api query to
    
    if(length(unavail) > 0) {
      
      for(i in 1:length(unavail)){
        #Step 2.1 : Handling files
        if(RCurl::url.exists(paste0(url_file,unavail[[i]],filter))){ 
          res <- jsonlite::fromJSON(paste0(url_file,unavail[[i]],filter))
          if (res[["notification"]] == "Success") {
            results <- res[["@graph"]]
						
url_encode = 'https://www.encodeproject.org/search/?' 
url_append = '&format=json&limit=all'
url_search = paste0(url_encode, 'searchTerm=', url_append)
url_file = paste0(url_encode, 'type=file&title=', url_append)
url_ds = paste0(url_encode, 'type=file&dataset=/', url_append)

RCurl::url.exists(paste0(url_search, unavail[[i]], '&format=json&limit=all'))

https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=TF+ChIP-seq&files.file_type=bed+idr_ranked_peak&assembly=GRCh38

my_url = paste0('https://www.encodeproject.org/search/?type=Experiment&control_type%21=%2A&status=released&perturbed=false&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=TF+ChIP-seq&files.file_type=bed+idr_ranked_peak', url_append)
RCurl::url.exists(my_url)
res = jsonlite::fromJSON(my_url)
exp_dataset <- res[["@graph"]][["@id"]]
// exp_dataset <- gsub(exp_dataset, pattern="/(.*)/.*/", replacement = "\\1")
exp_dataset <- gsub(exp_dataset, pattern="/.f*/(.*)/", replacement = "\\1")

res[["@graph"]]$files[[1]]


my_url = 'https://www.encodeproject.org/search/?type=Experiment&searchTerm=CTCF/&frame=object&format=json&limit=all'
RCurl::url.exists(my_url)
res = jsonlite::fromJSON(my_url)

my_query = '?type=Experiment&control_type%21=%2A&status=released&perturbed=false&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=TF+ChIP-seq&files.file_type=bed+idr_ranked_peak&type=File'

https://www.encodeproject.org/search/?type=Experiment&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens

https://www.encodeproject.org/search/?type=Experiment&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.file_type=bed+idr_ranked_peak

df_chip = encode_df %>% .[.$assay == 'TF ChIP-seq', ] %T>% pnrow # 141085
df_chip1 = df_chip %>% .[.$file_type == 'bed narrowPeak', ] %T>% pnrow # 30776
df_chip2 = df_chip1 %>% .[.$output_type == 'optimal IDR thresholded peaks'


url_append = '&format=json&limit=all'
my_query = 'assay_term_name=ChIP-seq&files.file_type=bed+idr_ranked_peak&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens'
my_url = paste0(url_encode, my_query, url_append)
RCurl::url.exists(my_url)


url_append = '&format=json&limit=all'
my_query = 'assay_term_name=ChIP-seq&files.file_type=bed+narrowPeak&files.output_type=optimal+IDR+thresholded+peaks&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens'
my_url = paste0(url_encode, my_query, url_append)
RCurl::url.exists(my_url)
res1 = jsonlite::fromJSON(my_url)
df = res1[['@graph']]



url_append = '&format=json&limit=all'
my_query = 'type=File&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens'
my_url = paste0(url_encode, my_query, url_append)
RCurl::url.exists(my_url)

res1 = jsonlite::fromJSON(my_url)
df = res1[['@graph']]


// https://www.encodeproject.org/search/?type=File&file_format_type=narrowPeak&file_format=bed&output_type=optimal+IDR+thresholded+peaks&assay_title=TF+ChIP-seq&assembly=GRCh38&status=released



url_append = '&format=json&limit=all'
my_query = 'type=File&file_format_type=narrowPeak&file_format=bed&output_type=optimal+IDR+thresholded+peaks&assay_title=TF+ChIP-seq&assembly=GRCh38&status=released'
my_url = paste0(url_encode, my_query, url_append)
RCurl::url.exists(my_url)

res1 = jsonlite::fromJSON(my_url)
df = res1[['@graph']]





https://www.encodeproject.org/search/?type=File&file_format=bed&status=released&output_type=optimal+IDR+thresholded+peaks&assay_title=TF+ChIP-seq&assembly=ce11

https://www.encodeproject.org/files/ENCFF734ALE/?format=json

// > df$status %>% table
// released
//     3665
// > df$control_type %>% table
// .
// control
//       6

 str(res,1)
 res[["@graph"]] %>% str(1)


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



process filtering_annotation_file {
  tag "${specie}"

  container = params.samtools_bedtools_perl

  publishDir path: "${specie}/genome/annotation", mode: 'link'

  input:
    set specie, file(gff3), file(fasta) from Genome_and_annotation

  output:
		set specie, file('anno_genes_only.gff3'), file('gff3_without_pseudogenes_and_ncRNA.gff3') into Channel_gff3_filtered

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
    set specie, file('chromosome_size.txt') into Chromosome_size_for_seqinfo
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
    set specie, specie_long, file('orgdb.sqlite') into Orgdb_for_annotation

  shell:
  '''

		#!/usr/bin/env Rscript

		library(AnnotationHub)

		specie = '!{specie}'
		specie_long = '!{specie_long}'
		ncbi_orgdb_version = '!{ncbi_orgdb_version}'
		annotationhub_cache = '!{params.annotationhub_cache}'

		specie_initial =  sapply(strsplit(specie_long, '_')[[1]], substr, 1, 1) %>% {paste0(toupper(.[1]), .[2])}

		ah = AnnotationHub(cache = annotationhub_cache)

		orgdb = query(ah, paste0('ncbi/standard/', ncbi_orgdb_version, '/org.', specie_initial, '.eg.db.sqlite'))[[1]]

		AnnotationDbi::saveDb(orgdb, 'orgdb.sqlite')

  '''

}


Orgdb_for_annotation
	.join(Chromosome_size_for_seqinfo)
	.join(Channel_gff3_filtered)
	.set{Channel_gff3_filtered_1}

process getting_R_annotation_files {
  tag "${specie}"

  container = params.bioconductor

  publishDir path: "${specie}/genome/annotation", mode: 'link'

  input:
		set specie, specie_long, orgdb, file(chr_size), file(gff3_genes_only), file(gff3_without_pseudogenes_and_ncRNA) from Channel_gff3_filtered_1

  output:
    file('*')

  shell:
  '''

		#!/usr/bin/env Rscript

		library(magrittr)
		library(AnnotationDbi)

		source('!{params.cactus_dir}/software/bin/export_df_to_bed.R')
		gff3_genes_only = '!{gff3_genes_only}'
		gff3_without_pseudogenes_and_ncRNA = '!{gff3_without_pseudogenes_and_ncRNA}'
		specie = '!{specie}'
		orgdb = AnnotationDbi::loadDb('!{orgdb}')
		df_chr_size = read.table('!{chr_size}', sep = '\t')
		specie_long = '!{specie_long}'

		# export annotation dataframe
		anno_df = rtracklayer::readGFF(gff3_genes_only) 
		entrez_id = mapIds(orgdb, keys = anno_df$gene_id,  column = 'ENTREZID', keytype = 'ENSEMBL', multiVals = 'first')
		anno_df %<>% dplyr::rename(gene_name = Name)
		anno_df %<>% dplyr::mutate(width = end - start)
		anno_df$entrez_id = entrez_id
		anno_df$chr = as.character(anno_df$seqid)
		anno_df %<>% dplyr::select(chr, start, end, width, strand, gene_name, gene_id, entrez_id)
		saveRDS(anno_df, 'df_genes_metadata.rds', version = 2)

		# create the seqinfo object
		nr = nrow(df_chr_size)
		seqinfo = GenomeInfoDb::Seqinfo(seqnames = df_chr_size[,1], seqlengths = df_chr_size[,2], isCircular = rep(F, nr) , genome = rep(specie_long, nr))
		saveRDS(seqinfo, 'seqinfo.rds')

		# export txdb
		txdb = GenomicFeatures::makeTxDbFromGFF(gff3_without_pseudogenes_and_ncRNA, chrominfo = seqinfo)
		AnnotationDbi::saveDb(txdb, file="txdb_without_pseudogenes_and_ncRNA.sqlite")

		# export gene vs transcripts 
		df_genes_transcripts = select(txdb, keys = keys(txdb), columns = c('GENEID', 'TXNAME'), keytype = 'GENEID')
		saveRDS(df_genes_transcripts, 'df_genes_transcripts.rds', version = 2)

		# extract promoters
		promoters = GenomicFeatures::promoters(GenomicFeatures::genes(txdb), upstream = 1500, downstream = 500)

		# export promoters as dataframe
		promoters_df = as.data.frame(promoters, stringsAsFactors = F)
		promoters_df$start[promoters_df$start < 0] = 0
		promoters_df$end[promoters_df$end < 0] = 0
		promoters_df %<>% dplyr::rename(gene_name = gene_id)
		colnames(promoters_df)[colnames(promoters_df) == 'gene_name'] = 'gene_id'
		promoters_df %<>% dplyr::inner_join(anno_df %>% dplyr::select(gene_id, gene_name, entrez_id), by = 'gene_id')
		promoters_df %<>% dplyr::arrange(seqnames, start)
		promoters_df %<>% dplyr::rename(chr = seqnames)
		saveRDS(promoters_df, file = 'promoters_df.rds', version = 2)

		# export promoters as bed file
		promoters_df1 = promoters_df
		promoters_df1$score = 0
		promoters_df1 %<>% dplyr::select(chr, start, end, gene_name, score, strand, gene_id)
		export_df_to_bed(promoters_df1, 'promoters.bed')


  '''

}

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
 
