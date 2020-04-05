BiocManager::install("biomaRt")
library(biomaRt)
mart = useEnsembl("ENSEMBL_MART_SNP",dataset = "hsapiens_snp",GRCh=37)
attributePages(mart)
listAttributes(mart)
listFilters(mart)
getBM(attributes = c("refsnp_id","polyphen_prediction", "polyphen_score", 
                     "sift_prediction", "sift_score","chrom_start","chrom_end"),
      filters = "snp_filter",
      values = c("rs138601174","rs782085786"),
      mart = mart)

#listDatasets(mart, verbose = FALSE)
#listMarts()
#listEnsembl(mart = NULL,host = "www.ensembl.org", version = NULL,
            #GRCh = 37, mirror = NULL, verbose = FALSE)
