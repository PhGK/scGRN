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
  if (tissue_type == 'AT1') {color <- 'green'}
  if (tissue_type == 'AT2') {color <- 'steelblue1'}
  if (tissue_type == 'Ciliated') {color <- 'yellow'}
  if (tissue_type == 'Club') {color <- 'pink'}
  color
}
choose_color <- Vectorize(choose_color)

quant_edges <- 0.8
test_data_tumor <- fread(paste('/mnt/scratch2/mlprot/Projekte/singlecell/data/high_values_concat.csv', sep = ""))

filtered_test_data <- test_data_tumor %>% filter(target_gene<source_gene) %>% dplyr::select(-c(V1, cell_id)) %>%
  #mutate('patient_id' = stringr::str_split(sample_name, '_', simplify = T)[,1])
  dplyr::mutate('patient_id' = substr(sample_name, 1,4)) %>%
  dplyr::filter(patient_id!= 'p029', patient_id != 'p028') #%>% # check below for cell count (no tumor cells?)
  filter(source_gene!='MALAT1', target_gene != 'MALAT1')

testa<- filtered_test_data %>% group_by(cell_type_epi, patient_id) %>%
  dplyr::summarize('meanv' = mean(LRP))

############
quantile_now_all_patients <- filtered_test_data %>% group_by(sample_name) %>% do(data.frame(t(quantile(.$LRP, c(quant_edges)))))
colnames(quantile_now_all_patients)[2] <- 'quant'
quantile_filtered_test_data_all_patients <- filtered_test_data %>% left_join(quantile_now_all_patients) %>% filter(LRP>quant)

graph_all_patients <- quantile_filtered_test_data_all_patients %>%
  dplyr::select(from =source_gene, to = target_gene, sample_name = sample_name, tissue_type = cell_type_epi) %>%
  graph_from_data_frame(directed=F)

####
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
  
resampled_data <- sample_n(quantile_filtered_test_data_all_patients, 3000)

#resampled_data %>% .$cell_type_epi %>% as.factor() %>% summary()

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
E(graph_all_patients_resampled)$color[all_edges$countvs<5] <- NA
set.seed(0)
l <- layout_with_fr(graph_all_patients_resampled, niter=5000)
png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_single_all_patients_all_tissue_resampled.png', width=4000, height=4000, res=100)
par(bg='black')
plot(graph_all_patients_resampled, vertex.size=0.1, edge.width=0.1, 
     edge.curved = curve_multiple(graph_all_patients_resampled, start = 0.03),
     vertex.label=ifelse( (degree(graph_all_patients_resampled) > 200) | (V(graph_all_patients_resampled)$name %in% c('CXCL1','CXCL2','CXCL3', 'NFKBIA', 'TSPAN1', 'SRI', 'IGFB7', 'MANF', 'CALR', 'HSPA5')), V(graph_all_patients_resampled)$name, NA),
     vertex.label.color = 'white',
     vertex.color=NA,
     vertex.label.cex= 1.0,
     layout=l,
     xlim=c(-0.7,0.3), ylim=c(-0.5,0.3))
dev.off()

nodes <- delete_edges(graph_all_patients_resampled, E(graph_all_patients_resampled))
V(nodes) %>% length()
###############################################################################
###############################################################################

patient_ids <- filtered_test_data$patient_id %>% unique()
patient_ids
id <- patient_ids[3]
png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_single_cells_concat_patients_alltisue.png', width=6000, height=6000)#, res=100)
par(bg = 'black')
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
  sample_n(5000) %>%
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

E(graph_one_patient_resampled)$color[one_patient_edges$countvs<5] <- NA

#png(paste('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_single_cells_one_patient_colored_',id, '.png', sep=''), width=2000, height=2000, res=100)
#par(bg='black')

plot(graph_one_patient_resampled, vertex.size=0.1, edge.width=0.1, 
     edge.curved = curve_multiple(graph_one_patient_resampled, start = 0.04),
     vertex.label=ifelse( (degree(graph_one_patient_resampled) > 100) | (V(graph_one_patient_resampled)$name %in% c('CXCL3', 'NFKBIA', 'TSPAN1', 'SRI', 'IGFB7', 'MANF', 'CALR', 'HSPA5')), V(graph_one_patient_resampled)$name, NA),
     vertex.label.color = 'white',
     vertex.label.cex= 1.0,
     layout=l,
     xlim=c(-0.8,0.4), ylim=c(-0.6,0.6))

#dev.off()

}
dev.off()
##########################################################
##########################################################


# analyze interaction numbers for all and individual patients
ribogroup <- c('H3F3A', 'FAU', 'RPL26', 'RPL22', 'RPL10', 'RPL5', 'RPS27', 'NACA', 'RACK1', 'EEF1A1', 'PTMA')
fosjungroup <- c('HSP90AA1', 'EGR1', 'RHOB', 'HSPH1', 'HSPB1', 'JUN', 'JUNB', 'FOS', 'FOSB', 'FUS', 'ATF3', 'IER2', 'IER3', 'INTS6', 'ID2', 'ID1', 'KLF4', 'KLF6')
rsph1group <- c('RSPH1', 'NQO1', 'SLC44A4', 'TSPAN1', 'MAPK10', 'MLF1', 'MUC16', 'IGFBP7', 'IGFBP2', 'IGF2BP2', 'IGFBP5')

ribodata <- resampled_data %>% dplyr::mutate('istumor' = ifelse(cell_type_epi == 'Tumor', 'Tumor', 'Normal')) %>%
  dplyr::filter(source_gene %in% ribogroup, target_gene %in% ribogroup)

fosjundata <- resampled_data %>% dplyr::mutate('istumor' = ifelse(cell_type_epi == 'Tumor', 'Tumor', 'Normal')) %>%
  dplyr::filter(source_gene %in% fosjungroup, target_gene %in% fosjungroup)

rsph1data <- resampled_data %>% dplyr::mutate('istumor' = ifelse(cell_type_epi == 'Tumor', 'Tumor', 'Normal')) %>%
  dplyr::filter(source_gene %in% rsph1group, target_gene %in% rsph1group)

ggplot(ribodata, aes(istumor, fill = patient_id, group = patient_id), color = 'black') + geom_bar(stat = 'count', position = 'dodge')
ggplot(fosjundata, aes(istumor, fill = patient_id, group = patient_id), color = 'black') + geom_bar(stat = 'count', position = 'dodge')
ggplot(rsph1data, aes(istumor, fill = patient_id, group = patient_id), color = 'black') + geom_bar(stat = 'count', position = 'dodge')

ggplot(ribodata, aes(cell_type_epi, fill = patient_id, group = patient_id), color = 'black') + geom_bar(stat = 'count', position = 'dodge')
ggplot(fosjundata, aes(cell_type_epi, fill = patient_id, group = patient_id), color = 'black') + geom_bar(stat = 'count', position = 'dodge')
ggplot(rsph1data, aes(cell_type_epi, fill = patient_id, group = patient_id), color = 'black') + geom_bar(stat = 'count', position = 'dodge')

png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/RSPH1values.png', height=800, width=800)
ggplot(rsph1data, aes(cell_type_epi, fill = patient_id, group = patient_id), color = 'black') + 
  geom_bar(stat = 'count', position = 'dodge') +
  scale_x_discrete(drop=FALSE)
dev.off()

png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/ribodatavalues.png', height=800, width=800)
ggplot(ribodata, aes(cell_type_epi, fill = patient_id, group = patient_id), color = 'black') + 
  geom_bar(stat = 'count', position = 'dodge') +
  scale_x_discrete(drop=FALSE)
dev.off()

png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/fosjundatavalues.png', height=800, width=800)
ggplot(fosjundata, aes(cell_type_epi, fill = patient_id, group = patient_id), color = 'black') + 
  geom_bar(stat = 'count', position = 'dodge') +
  scale_x_discrete(drop=FALSE)
dev.off()

ggplot(ribodata, aes(patient_id, fill = istumor), color = 'black') + geom_bar(stat = 'count', position = 'dodge')
ggplot(fosjundata, aes(patient_id, fill = istumor), color = 'black') + geom_bar(stat = 'count', position = 'dodge')
ggplot(rsph1data, aes(patient_id, fill = istumor), color = 'black') + geom_bar(stat = 'count', position = 'dodge')

######################################################
#cell counts, validity
unique_cells <- test_data_tumor %>% 
  group_by(cell_id, cell_type_epi) %>% 
  dplyr::summarize('meanLRP' = mean(LRP))

test_data_tumor %>% select(-c(LRP, source_gene, target_gene, V1, 'Unnamed: 0')) %>% unique() %>% dim()

(test_data_tumor %>% 
    dplyr::select(-c(LRP, source_gene, target_gene)) %>% 
    unique %>% dim())[1] /500 

count_data <- fread('/mnt/scratch2/mlprot/Projekte/singlecell/data/epi_top2000.csv')[,-c(5:2000)]
head(count_data)
count_data$cell_id %>% unique() %>% length
count_data$Row.names %>% unique() %>% length
###################################
#cell_count
unique_cells <- filtered_test_data %>% 
  group_by(sample_name) %>% 
  dplyr::summarize('meanLRP' = mean(LRP), patient_id, cell_type_epi)

ggplot(unique_cells, aes(patient_id, fill = cell_type_epi, group = cell_type_epi), color = 'black') + 
  geom_bar(stat = 'count', position = 'dodge') +
  scale_x_discrete(drop=FALSE)
###################################  
ggplot(rsph1data, aes(cell_type_epi, fill = patient_id, group = patient_id), color = 'black') + 
  geom_bar(stat = 'count', position = 'dodge') +
  scale_x_discrete(drop=FALSE)
######################
######################
#overlap between normal cells and tumor cells
library(vecsets)
resampled_data$cell_type_epi %>% as.factor() %>% summary()
cell_types <- resampled_data$cell_type_epi %>% unique
interaction_data <- resampled_data %>% unite('interaction', c(source_gene, target_gene))

compare_overlap <- function(cell_type) {
  tumor_data <- interaction_data %>% filter(cell_type_epi == 'Tumor') %>% .$interaction
  print(tumor_data %>% dim)
  other_data <- interaction_data %>% filter(cell_type_epi == cell_type) %>% .$interaction
  c(cell_type, dim(intersect(tumor_data, other_data))[1])
  #dim(tumor_data[tumor_data %in% intersect(tumor_data, other_data)])[1]
  c(cell_type, length(tumor_data[tumor_data %in%intersect(tumor_data, other_data)]))
  
}

lapply(cell_types, compare_overlap)
