
get_encode_df <- function(my_query){
  my_url = paste0(url_search, my_query, url_append)
  if(!RCurl::url.exists(my_url)) stop('url doesnt exist')
  res = jsonlite::fromJSON(my_url)
  df = res[['@graph']]
  # type = gsub(my_query, pattern = 'type=(.*?)&.*', replacement = '\\1') %>% tolower
  type = gsub(df[1,'@id'], pattern="/(.*)/.*/", replacement = "\\1") %>% gsub('-', '_', .)
  colnames(df)[1] = type
  return(df)
}

pnrow <- function(x) print(nrow(x))

collapse_slims <- function(x) map_chr(x, ~paste(.x, collapse = ', '))

capture_group <- function(x, pattern) gsub(pattern, '\\1', x)
