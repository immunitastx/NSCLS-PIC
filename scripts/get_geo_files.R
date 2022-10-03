#### Download data from GEOQuery ####
library(GEOquery)
library(tidyverse)
gse_data <- "GSE160903"

geo_dat <- GEOquery::getGEO(gse_data)
View(geo_dat$`GSE160903-GPL24247_series_matrix.txt.gz`@phenoData@data$geo_accession)
gsm <- lapply(geo_dat, function(x){
  x@phenoData@data$geo_accession
}) %>%
  unlist() %>%
  unname() %>% unique()

#GEOquery::getGEOSuppFiles()

for (cur_gsm in gsm) {
  GEOquery::getGEOSuppFiles(GEO = cur_gsm,baseDir = "output/umi.tab/",makeDirectory = FALSE)
}
Sys.chmod("output/umi.tab/*",mode = "0777",)
system("ls output/umi.tab/")
system("gzip -d * ls -1 | awk -F'_' '{print $0,$2}' | xargs -n 2 mv output/umi.tab/ output/umi.tab/")
