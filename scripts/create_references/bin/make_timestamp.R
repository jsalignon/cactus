

wl <- function(x) writeLines(x)
wl1 <- function(x) writeLines(paste0(x, "\n"))
wl2 <- function(x) writeLines(paste0("\n", x, "\n"))

sink("timestamp.txt")

  paste0("This reference was built on ", Sys.Date(), ".") %>% wl2

  "The version of the workflow and package manager tools used are:" %>% wl

  system("nextflow -version", intern = T) %>% wl
  system("singularity --version", intern = T) %>% wl1

sink()
