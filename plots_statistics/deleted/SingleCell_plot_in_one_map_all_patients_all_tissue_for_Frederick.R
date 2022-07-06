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


choose_color <- function(tissue_type, mode='all') {
  if (tissue_type == 'Tumor') {color <- 'red'}
  if (tissue_type == 'AT1') {color <- 'steelblue'}
  if (tissue_type == 'AT2') {color <- 'steelblue'}
  if (tissue_type == 'Ciliated') {color <- 'steelblue'}
  if (tissue_type == 'Neuroendocrine') {color <- 'steelblue'}
  if (tissue_type == 'Club') {color <- 'steelblue'}
  color
}
choose_color <- Vectorize(choose_color)

quant_edges <- 0.9
test_data_tumor <- fread(paste('/mnt/scratch2/mlprot/Projekte/singlecell/data/high_values_concat.csv', sep = ""))

filtered_test_data <- test_data_tumor %>% filter(target_gene<source_gene) %>% dplyr::select(-c(V1, cell_id)) %>%
  #mutate('patient_id' = stringr::str_split(sample_name, '_', simplify = T)[,1])
  dplyr::mutate('patient_id' = substr(sample_name, 1,4)) %>%
  dplyr::filter(patient_id!= 'p029') %>%
  filter(source_gene!='MALAT1', target_gene != 'MALAT1')




testa<- filtered_test_data %>% group_by(cell_type_epi, patient_id) %>%
  dplyr::summarize('meanv' = mean(LRP))

############
quantile_now_all_patients <- filtered_test_data %>% group_by(sample_name) %>% do(data.frame(t(quantile(.$LRP, c(quant_edges)))))
colnames(quantile_now_all_patients)[2] <- 'quant'
quantile_filtered_test_data_all_patients <- filtered_test_data %>% left_join(quantile_now_all_patients) %>% filter(LRP>quant)

################
#quantify interations
cell_counts <- quantile_filtered_test_data_all_patients %>%dplyr::select(cell_type_epi, sample_name) %>% 
  unique() %>% 
  dplyr::group_by(cell_type_epi) %>% 
  dplyr::summarize('cell_count' = n())

quantification <- quantile_filtered_test_data_all_patients %>% 
  group_by(source_gene, target_gene, cell_type_epi) %>%
  dplyr::summarize('interactioncount' = n()) %>%
  left_join(cell_counts) %>%
  dplyr::mutate('relativecount' = interactioncount/cell_count)



cell_counts <- quantile_filtered_test_data_all_patients$sample_name %>% unique %>% length()
quantification_all <- quantile_filtered_test_data_all_patients %>% 
  group_by(source_gene, target_gene) %>%
  dplyr::summarize('interactioncount' = n()) %>%
  dplyr::mutate('relative_count' = interactioncount /cell_counts)

####
#quantify JUN-FOS network
junfos <- c('HSP90AA1','HSP90AB1', 'EGR1', 'RHOB', 'HSPH1', 'HSPB1', 'JUN', 'JUNB', 'FOS', 'FOSB', 'FUS', 'ATF3', 'IER2', 'IER3', 'INTS6', 'ID2', 'ID1', 'KLF4', 'KLF6')

quantification_junfos <- quantification_all %>% dplyr::filter((source_gene %in% junfos)|  (target_gene %in% junfos))

####
####
#quantify ID2-ID4 network
ID2ID4 <- c('ID2', 'ID4', 'IER3', 'HEXIM1', 'INTS6', 'BRD2', 'GADD45G', 'KLF4')

quantificationID2ID4 <- quantification %>% dplyr::filter((source_gene %in% ID2ID4) &  (target_gene %in% ID2ID4))
####
#quantify ID2-ID4 network
RSPH1 <-c('RSPH1', 'SPA17', 'ZMYND10', 'SOX2', 'HSP90AA1', 'GSTA1', 'SSBP4', 'RUVBL1', 'RUVBL2', 'IGFBP7', 'SLC44A4', 'TSPAN1', 'MAPK10', 'MLF1', 'MUC16', 'IGFBP7','MAPK10', 'AKAP9', 'IGFBP5')

quantificationRSPH1 <- quantification %>% dplyr::filter((source_gene %in% RSPH1) &  (target_gene %in% RSPH1))

####
cxcl <- c('CXCL1', 'CXCL2', 'CXCL3', 'NFKBIA')

quantification_cxcl <- quantification %>% dplyr::filter((source_gene %in% cxcl)&  (target_gene %in% cxcl))

#####
#quantification tumor interactions in different patients
cell_counts_tumor <- quantile_filtered_test_data_all_patients %>%dplyr::select(cell_type_epi, sample_name, patient_id) %>% 
  unique() %>% 
  dplyr::filter(cell_type_epi == 'Tumor') %>% 
  dplyr::group_by(patient_id) %>%
  dplyr::summarize('cell_count_tumor' = n())

quantification_tumor <- quantile_filtered_test_data_all_patients %>% 
  group_by(source_gene, target_gene, cell_type_epi, patient_id) %>%
  dplyr::summarize('interactioncount' = n()) %>%
  left_join(cell_counts_tumor) %>%
  dplyr::mutate('relativecount' = interactioncount/cell_count_tumor) %>%
  dplyr::filter(cell_count_tumor>50, interactioncount >5,relativecount>0.2, cell_type_epi=='Tumor')

######

ntissue <- quantile_filtered_test_data_all_patients %>% dplyr::select(sample_name, cell_type_epi)%>%
  unique() %>%
  .$cell_type_epi %>%as.factor %>% summary()
ntissue

return_sample <- function(ratio, otherdataset) {
  nsamples <- dim(otherdataset)[1]
  if (ratio <=1) {
    new_data = otherdataset[sample(nsamples,as.integer(nsamples*ratio), replace=F),]
  }
  else {
    if (nsamples<1){
      new_data <- otherdataset
    }
    else {
    new_data = otherdataset[sample(nsamples,as.integer(nsamples*ratio), replace=T),]
    }
  }
  new_data
}

return_sample(0.1, quantile_filtered_test_data_all_patients %>% dplyr::group_by(cell_type_epi)) %>%
  .$cell_type_epi %>% as.factor %>% summary()
    

sample_n <- function(dataset, N) {
  TUMORdata <- dataset%>% filter(cell_type_epi == 'Tumor')
  AT1data <- dataset%>% filter(cell_type_epi == 'AT1')
  AT2data <- dataset %>% filter(cell_type_epi == 'AT2')
  Ciliateddata <- dataset%>% filter(cell_type_epi == 'Ciliated')
  Clubdata <- dataset %>% filter(cell_type_epi == 'Club')
  Neuroendocrinedata <- dataset%>% filter(cell_type_epi == 'Neuroendocrine')
  
  
  nTUMOR <- dim(TUMORdata)[1]
  nAT1 <- dim(AT1data)[1]
  nAT2 <- dim(AT2data)[1]
  nCiliated <- dim(Ciliateddata)[1]
  nClub <- dim(Clubdata)[1]
  nNeuroendocrine <- dim(Neuroendocrinedata)[1]
  
  ratioTUMOR <- 5*N/nTUMOR
  ratioAT1 <- N/nAT1
  ratioAT2 <- N/nAT2
  ratioCiliated <- N/nCiliated
  ratioClub <- N/nClub
  ratioNeuroendocrine <- N/nNeuroendocrine
  
  new_TUMORdata <- TUMORdata %>% 
    group_by(source_gene, target_gene) %>%
    return_sample(ratioTUMOR,.)
  
  new_AT1data <- AT1data %>% 
    group_by(source_gene, target_gene) %>%
    return_sample(ratioAT1,.)
  
  new_AT2data <- AT2data %>% 
    group_by(source_gene, target_gene) %>%
    return_sample(ratioAT2,.)
  
  new_Ciliateddata <- Ciliateddata %>% 
    group_by(source_gene, target_gene) %>%
    return_sample(ratioCiliated,.)
  
  new_Clubdata <- Clubdata %>% 
    group_by(source_gene, target_gene) %>%
    return_sample(ratioClub,.)
  
  new_Neuroendocrinedata <- Neuroendocrinedata %>% 
    group_by(source_gene, target_gene) %>%
    return_sample(ratioNeuroendocrine,.)
  
  rbind(new_TUMORdata, new_AT1data, new_AT2data, new_Ciliateddata, new_Clubdata, new_Neuroendocrinedata)
  } 
  
resampled_data <- sample_n(quantile_filtered_test_data_all_patients, 5000)

resampled_data %>% .$cell_type_epi %>% as.factor() %>% summary()

####
graph_all_patients_resampled <- resampled_data %>%
  dplyr::select(from =source_gene, to = target_gene, sample_name = sample_name, tissue_type = cell_type_epi) %>%
  graph_from_data_frame(directed=F)

E(graph_all_patients_resampled)$color <- choose_color(E(graph_all_patients_resampled)$tissue_type, 'all')
#E(graph_all_patients_resampled)$color[!which_multiple(graph_all_patients_resampled)] <- NA

all_edges <- ends(graph_all_patients_resampled, E(graph_all_patients_resampled))  %>% data.frame()
colnames(all_edges) <- c('source_gene', 'target_gene')
all_edges <-all_edges %>%
  dplyr::mutate('dummy' = 1) %>%
  group_by (source_gene, target_gene) %>%
  dplyr::mutate('countvs' = n())

E(graph_all_patients_resampled)$weight <- 0.01
E(graph_all_patients_resampled)$color[all_edges$countvs<7] <- NA
set.seed(3)
l <- layout_with_fr(graph_all_patients_resampled)
png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_allpatients_equal_proportion.png', width=4000, height=4000, res=100)
par(bg='black')
plot(graph_all_patients_resampled, vertex.size=0.1, edge.width=0.1, 
     edge.curved = curve_multiple(graph_all_patients_resampled, start = 0.03),
     vertex.label=ifelse( (degree(graph_all_patients_resampled) > 50) | (V(graph_all_patients_resampled)$name %in% c('CXCL1','CXCL2','CXCL3', 'NFKBIA','HSP90AA1','HSP90AB1', 'TSPAN1', 'SRI', 'IGFB7', 'MANF', 'CALR', 'HSPA5')), V(graph_all_patients_resampled)$name, NA),
     vertex.label.color = 'white',
     vertex.color=NA,
     vertex.label.cex= 1.0,
     layout=l,
     xlim=c(-0.5,0.1), ylim=c(-0.4,0.2))
dev.off()

nodes <- delete_edges(graph_all_patients_resampled, E(graph_all_patients_resampled))
V(nodes) %>% length()
###############################################################################
###############################################################################


patient_ids <- filtered_test_data$patient_id %>% unique()
patient_ids
id <- patient_ids[3]
png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_single_cells_concat_patients_alltisue.png', width=2000, height=2000, res=100)
#par(bg = 'black')
par(mfrow = c(3,3))
for (id in patient_ids[1:9]) {
print(id)
quantile_filtered_test_data_onepatient <- quantile_filtered_test_data_all_patients %>% filter(patient_id == id)

###
'
graph_one_patient_resampled <- quantile_filtered_test_data_onepatient %>% 
  sample_n(20000) %>%
  dplyr::select(from =source_gene, to = target_gene, sample_name = sample_name, tissue_type = tissue_type) #%>%
  graph_from_data_frame(directed=F)
'
quantile_filtered_test_data_onepatient$cell_type_epi %>% as.factor() %>% summary()

frame_one_patient_resampled <- quantile_filtered_test_data_onepatient %>% 
  sample_n(1000) %>%
  dplyr::select(from =source_gene, to = target_gene, sample_name = sample_name, tissue_type = cell_type_epi)

graph_one_patient_resampled <- nodes

for (i in seq(dim(frame_one_patient_resampled)[1])) {
  try({
  graph_one_patient_resampled <- graph_one_patient_resampled %>% add_edges(c(frame_one_patient_resampled$from[i], 
                                                                           frame_one_patient_resampled$to[i]),
                                                                           sample_name = frame_one_patient_resampled$sample_name[i], 
                                                                           tissue_type = frame_one_patient_resampled$tissue_type[i])
  })
}

### 
E(graph_one_patient_resampled)$color <- choose_color(E(graph_one_patient_resampled)$tissue_type, 'both')

one_patient_edges <- ends(graph_one_patient_resampled, E(graph_one_patient_resampled))  %>% data.frame()
colnames(one_patient_edges) <- c('source_gene', 'target_gene')
one_patient_edges <-one_patient_edges %>%
  dplyr::mutate('dummy' = 1) %>%
  group_by (source_gene, target_gene) %>%
  dplyr::mutate('countvs' = n())

E(graph_one_patient_resampled)$color[one_patient_edges$countvs<2] <- NA

#png(paste('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_single_cells_one_patient_colored_',id, '.png', sep=''), width=2000, height=2000, res=100)
#par(bg='black')

plot(graph_one_patient_resampled, vertex.size=0.1, edge.width=0.1, 
     edge.curved = curve_multiple(graph_one_patient_resampled, start = 0.02),
     vertex.label=ifelse( (degree(graph_one_patient_resampled) > 230), V(graph_one_patient_resampled)$name, NA),
     vertex.label.color = 'white',
     vertex.label.cex= 1.5,
     layout=l,
     xlim=c(-0.5,0.5), ylim=c(-0.6,0.4))

#dev.off()

}
dev.off()

