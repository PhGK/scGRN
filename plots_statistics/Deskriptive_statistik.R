setwd('/mnt/scratch2/mlprot/Projekte/singlecell')
library(ggplot2)
library(stringr)
library(magrittr)
library(tidyr)
library(dplyr)
library(DescTools)
library(gplots)
library(ComplexHeatmap)
library(pbmcapply)
library(circlize)
library(data.table)
library(Hmisc)
library(stringi)
library(mvMORPH)
library(stats)
library(xtable)
library(plyr)
library(Rtsne)
library(igraph)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

expr_data <- fread('../data/epi_top2000.csv') %>% dplyr::select(patient_id, cell_id, cell_type_epi)
LRP_data <- fread('../results/high_values_concat.csv') %>% dplyr::select(patient_id, cell_id, cell_type_epi)
celltypes <- expr_data$cell_type_epi %>%
  as.factor()
summary(celltypes)

patients <- LRP_data$patient_id %>% unique()

sapply(patients, function(idx) LRP_data %>% filter(patient_id == idx)%>% unique %>% .$cell_type_epi %>% as.factor() %>% summary())

sapply(patients, function(idx) expr_data %>% filter(patient_id == idx)%>% unique %>% .$cell_type_epi %>% as.factor() %>% summary())

# check cell numbers of LRP values

LRP_data$cell_id %>% unique() %>% length()

