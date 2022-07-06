library(magrittr)
library(pROC)
library(dplyr)
library(biomaRt)

#####

data <- read.csv('../data/epi_top2000.csv')
gene_names <- colnames(data[3:2002])



########
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "useast.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol",
    "entrezgene_id",
    "ensembl_gene_id",
    "gene_biotype"),
  filter = "hgnc_symbol",
  values = gene_names,
  uniqueRows=TRUE)

selected_genes <- annotLookup %>% dplyr::filter(gene_biotype=='protein_coding')


write.csv(selected_genes$hgnc_symbol, '../data/selected_genes.csv', col.names=F, row.names=F)
