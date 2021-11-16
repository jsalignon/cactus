
# this function allows to read arrays into R

read_from_nextflow <- function(x) substr(x, 2, nchar(x) - 1) %>% gsub(' ', '', .) %>% strsplit(., ',') %>% .[[1]]
