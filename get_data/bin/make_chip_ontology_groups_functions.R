

get_ontology_chip_names <- function(dt, ontology_id, ontology_name){
  
  ontology_values = dt[[ontology_id]]
  
  # making a dictionary of all bed chip_names for the ontology
  all_ontologies = strsplit(unique(ontology_values), ', ') %>% unlist %>% unique
  if(length(all_ontologies) == 0) return(list())
  all_ontologies %<>% setNames(., paste0(ontology_name, '.', .) %>% gsub(' ', '_', .) %>% gsub('-', '', .))
  l_ontology_chip_names = map(all_ontologies, ~dt[grepl(.x, ontology_values)]$local_file)
  
  # keeping only ontologies with at least 10 chip_names 
  ontologies_to_keep = which(map_int(l_ontology_chip_names, length) > 10) %>% names
  l_ontology_chip_names %<>% .[ontologies_to_keep]
  
  return(l_ontology_chip_names)
}
  