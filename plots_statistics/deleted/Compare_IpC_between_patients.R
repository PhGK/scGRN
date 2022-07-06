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
###
# - compute modules
# - build patient-wise heatmap corrected for dropout 


quant_edges <- 0.8
test_data_tumor <- fread(paste('../usedata/high_values_concat.csv', sep = ""))

filtered_test_data <- test_data_tumor %>% filter(target_gene<source_gene) %>% dplyr::select(-c(V1, cell_id)) %>%
  #mutate('patient_id' = stringr::str_split(sample_name, '_', simplify = T)[,1])
  dplyr::mutate('patient_id' = substr(sample_name, 1,4)) %>%
  dplyr::filter(patient_id!= 'p029', patient_id != 'p028') %>% # check below for cell count (no tumor cells?)
  filter(source_gene!='meanpercell', target_gene != 'meanpercell') %>%
  filter(cell_type_epi != 'Neuroendocrine')

quantile_now_all_patients <- filtered_test_data %>% group_by(sample_name) %>% do(data.frame(t(quantile(.$LRP, c(quant_edges)))))
colnames(quantile_now_all_patients)[2] <- 'quant'
quantile_filtered_test_data_all_patients <- filtered_test_data %>% left_join(quantile_now_all_patients) %>% filter(LRP>quant)

check_patients <- quantile_filtered_test_data_all_patients %>% dplyr::select(patient_id, cell_type_epi) %>%
  group_by(patient_id, cell_type_epi) %>%
  unique() 

cell_counts <-quantile_filtered_test_data_all_patients %>% dplyr::select(sample_name,cell_type_epi) %>%
  unique() %>%
  group_by(cell_type_epi) %>%
  dplyr::summarize('counts' = n())

######################
######################

patient_ids <- quantile_filtered_test_data_all_patients %>% dplyr::select(patient_id) %>% unique() 

get_relative_counts <- function(patient, cell_type) {
  
  patient_set <- quantile_filtered_test_data_all_patients %>% filter(patient_id==patient)
  
  ncells_perid <- patient_set$sample_name %>%
    unique %>%
    length() 
    #filter(cell_type_epi == cell_type) %>%
    #dplyr::select(patient_id, sample_name) #%>%
    #unique() %>%
    #dplyr::summarize('all_cellcount' = n())
  
  
  count_data <- patient_set %>% 
    filter(cell_type_epi == cell_type) %>%
    dplyr::select(patient_id, sample_name, source_gene, target_gene) %>%
    group_by(source_gene, target_gene) %>%
    dplyr::summarize('spec_cell_count' = n())
    
    count_data$all_cellcount <- ncells_perid 
    countdata <- count_data %>% dplyr::mutate('ratio' = spec_cell_count/all_cellcount, 'patient_id' = patient)
    countdata  
    
}

get_relative_counts(patient_ids$patient_id[1], 'Tumor')

get_all_patients <- function(cell_type) {
    all_patients <- rbindlist(lapply(patient_ids$patient_id, function(patient) get_relative_counts(patient, cell_type)))
    all_patients$cell_type <- cell_type
    all_patients
}




patient_mutations <- data.frame(patient_id = c('p018', 'p019','p023','p024','p027','p030','p031','p032','p033','p034'),
                                'KRAS' = c(1,0,1,0,0,1,1,1,1,0),
                                'TP53' = c(0,0,1,0,1,0,0,0,0,0),
                                'PIK3CA' = c(0,0,0,0,0,0,1,0,0,0),
                                'EML4' = c(0,0,0,1,0,0,0,0,0,0))


all_data <- rbindlist(lapply(cell_counts$cell_type_epi, get_all_patients))
all_data <- left_join(all_data, patient_mutations)

current_patient <- all_data %>% dplyr::filter(cell_type == 'Tumor', patient_id == patient_ids$patient_id[11])  %>%
  dplyr::arrange(desc(ratio))
