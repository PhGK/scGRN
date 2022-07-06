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


choose_color <- function(tissue_type, mode='both') {
  if (mode == 'both') {colors <- c('red', 'steelblue1')}
  if (mode == 'Tumor') {colors <- c('red', 'black')}
  if (mode == 'Normal') {colors <- c('black', 'steelblue')}
  ifelse(tissue_type=='Tumor', colors[1], colors[2])
}



quant_edges <- 0.9
test_data_tumor <- fread(paste('/mnt/scratch2/mlprot/Projekte/singlecell/data/high_values_concat.csv', sep = ""))

filtered_test_data <- test_data_tumor %>% filter(target_gene<source_gene) %>% dplyr::select(-c(V1, cell_id)) %>%
  #mutate('patient_id' = stringr::str_split(sample_name, '_', simplify = T)[,1])
  dplyr::mutate('patient_id' = substr(sample_name, 1,4)) %>%
  dplyr::mutate(tissue_type = ifelse(cell_type_epi == 'Tumor', 'Tumor', 'Normal')) %>%
  dplyr::filter(patient_id!= 'p029')


testa<- filtered_test_data %>% group_by(tissue_type, patient_id) %>%
  dplyr::summarize('meanv' = mean(LRP))

count_edges <- filtered_test_data %>% group_by(source_gene, target_gene) %>%
  dplyr::summarize('counts' = n()) %>% filter( (source_gene =='CXCL1') |  (target_gene == 'CXCL1'))

############
quantile_now_all_patients <- filtered_test_data %>% group_by(sample_name) %>% do(data.frame(t(quantile(.$LRP, c(quant_edges)))))
colnames(quantile_now_all_patients)[2] <- 'quant'
quantile_filtered_test_data_all_patients <- filtered_test_data %>% left_join(quantile_now_all_patients) %>% filter(LRP>quant)

graph_all_patients <- quantile_filtered_test_data_all_patients %>%
  dplyr::select(from =source_gene, to = target_gene, sample_name = sample_name, tissue_type = tissue_type) %>%
  graph_from_data_frame(directed=F)

E(graph_all_patients)$color <- choose_color(E(graph_all_patients)$tissue_type, 'both')
E(graph_all_patients)$color[!which_multiple(graph_all_patients)] <- NA

l <- layout_with_fr(graph_all_patients)
png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_single_all_patients.png', width=4000, height=4000, res=100)
par(bg='black')
plot(graph_all_patients, vertex.size=0.1, edge.width=0.1, 
     edge.curved = curve_multiple(graph_all_patients, start = 0.1),
     vertex.label=ifelse( (degree(graph_all_patients) > 600) | (V(graph_all_patients)$name %in% c('CXCL3', 'NFKBIA', 'TSPAN1', 'SRI', 'IGFB7', 'MANF', 'CALR', 'HSPA5')), V(graph_all_patients)$name, NA),
     vertex.label.color = 'white',
     vertex.label.cex= 1.2,
     layout=l)
dev.off()

####
ntissue <- quantile_filtered_test_data_all_patients %>% dplyr::select(sample_name, tissue_type)%>%
  unique() %>%
  .$tissue_type %>%as.factor %>% summary()
ntissue

return_sample <- function(ratio, otherdataset) {
  nsamples <- dim(otherdataset)[1]
  if (ratio <=1) {
    new_data = otherdataset[sample(nsamples,as.integer(nsamples*ratio), replace=F),]
  }
  else {
    new_data = otherdataset[sample(nsamples,as.integer(nsamples*ratio), replace=T),]
  }
  new_data
}

    

sample_n <- function(dataset, N) {
  TUMORdata <- dataset%>% filter(tissue_type == 'Tumor')
  NORMALdata <- dataset %>% filter(tissue_type == 'Normal')
  nTUMOR <- dim(TUMORdata)[1]
  nNORMAL <- dim(NORMALdata)[1]
  ratioTUMOR <- N/nTUMOR
  ratioNORMAL <- N/nNORMAL
  
  new_TUMORdata <- TUMORdata %>% 
    group_by(source_gene, target_gene) %>%
    return_sample(ratioTUMOR,.)
  
  new_NORMALdata <- NORMALdata %>% 
    group_by(source_gene, target_gene) %>%
    return_sample(ratioNORMAL,.)
  
  rbind(new_TUMORdata, new_NORMALdata)
  } 

set.seed(0) 
resampled_data <- sample_n(quantile_filtered_test_data_all_patients, 10000)
####
graph_all_patients_resampled <- resampled_data %>%
  dplyr::select(from =source_gene, to = target_gene, sample_name = sample_name, tissue_type = tissue_type) %>%
  graph_from_data_frame(directed=F)

E(graph_all_patients_resampled)$color <- choose_color(E(graph_all_patients_resampled)$tissue_type, 'both')
E(graph_all_patients_resampled)$color[!which_multiple(graph_all_patients_resampled)] <- NA

all_edges <- ends(graph_all_patients_resampled, E(graph_all_patients_resampled))  %>% data.frame()
colnames(all_edges) <- c('source_gene', 'target_gene')
all_edges <-all_edges %>%
  dplyr::mutate('dummy' = 1) %>%
  group_by (source_gene, target_gene) %>%
  dplyr::mutate('countvs' = n())

E(graph_all_patients_resampled)$color[all_edges$countvs<3] <- NA
set.seed(2)
l <- layout_with_fr(graph_all_patients_resampled)
png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_single_all_patients_resampled.png', width=4000, height=4000, res=100)
par(bg='black')
plot(graph_all_patients_resampled, vertex.size=0.1, edge.width=0.2, 
     edge.curved = curve_multiple(graph_all_patients_resampled, start = 0.03),
     vertex.label=ifelse( (degree(graph_all_patients_resampled) > 150) | (V(graph_all_patients_resampled)$name %in% c('CXCL1','CXCL2','CXCL3', 'NFKBIA', 'TSPAN1', 'SRI', 'IGFB7', 'MANF', 'CALR', 'HSPA5')), V(graph_all_patients_resampled)$name, NA),
     vertex.label.color = 'white',
     vertex.color='white',
     vertex.shape = 'none',
     vertex.label.cex= 1.0,
     layout=l,
     xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
dev.off()

nodes <- delete_edges(graph_all_patients_resampled, E(graph_all_patients_resampled))
V(nodes) %>% length()
###############################################################################
#histogram of highest interactions
histo_data <- resampled_data %>% unite('interaction', c(source_gene, target_gene)) %>%
  dplyr::select(interaction, cell_type_epi) %>% 
  dplyr::mutate('istumor' = ifelse(cell_type_epi == 'Tumor','Tumor', 'Normal')) %>% 
  group_by(interaction, istumor) %>%
  dplyr::mutate(ncount = n()) #%>%
  filter(ncount>1)

ggplot(histo_data, aes(reorder(ncount, -ncount),fill=istumor, group = istumor)) + 
  geom_histogram(stat='count', position = 'stack', alpha = 0.6) 

###############################################################################

patient_ids <- filtered_test_data$patient_id %>% unique()
patient_ids
last_patient <- quantile_filtered_test_data_all_patients %>% filter(patient_id == patient_ids[9]) %>%
  dplyr::filter(cell_type_epi == 'Tumor') %>%
  dplyr::select(sample_name) %>%
  unique()
  
png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_single_cells_concat_patients.png', width=3000, height=3000, res=100)
par(bg = 'black')
par(mfrow = c(3,4))
for (id in patient_ids) {

quantile_filtered_test_data_onepatient <- quantile_filtered_test_data_all_patients %>% filter(patient_id == id)

###
'
graph_one_patient_resampled <- quantile_filtered_test_data_onepatient %>% 
  sample_n(20000) %>%
  dplyr::select(from =source_gene, to = target_gene, sample_name = sample_name, tissue_type = tissue_type) #%>%
  graph_from_data_frame(directed=F)
'
frame_one_patient_resampled <- quantile_filtered_test_data_onepatient %>% 
  sample_n(5000) %>%
  dplyr::select(from =source_gene, to = target_gene, sample_name = sample_name, tissue_type = tissue_type)

graph_one_patient_resampled <- nodes

for (i in seq(dim(frame_one_patient_resampled)[1])) {
  print(i)
  
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
     edge.curved = curve_multiple(graph_one_patient_resampled, start = 0.03),
     vertex.label=ifelse( (degree(graph_one_patient_resampled) > 200) | (V(graph_one_patient_resampled)$name %in% c('CXCL3', 'NFKBIA', 'TSPAN1', 'SRI', 'IGFB7', 'MANF', 'CALR', 'HSPA5')), V(graph_one_patient_resampled)$name, NA),
     vertex.label.color = 'white',
     vertex.label.cex= 1,
     vertex.shape = 'none',
     layout=l,
     xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))

#dev.off()

}
dev.off()

##############################################
##############################################
# analyze interaction numbers for all and individual patients
resampled_data$cell_type_epi <- as.factor(resampled_data$cell_type_epi)
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

ggplot(ribodata, aes(cell_type_epi, fill = patient_id, group = patient_id), color = 'black') + geom_bar(stat = 'count', position = 'dodge')+
  scale_x_discrete(drop=FALSE)
ggplot(fosjundata, aes(cell_type_epi, fill = patient_id, group = patient_id), color = 'black') + geom_bar(stat = 'count', position = 'dodge')+
  scale_x_discrete(drop=FALSE)
png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/RSPH1values.png', height=800, width=800)
ggplot(rsph1data, aes(cell_type_epi, fill = patient_id, group = patient_id), color = 'black') + 
  geom_bar(stat = 'count', position = 'dodge') +
  scale_x_discrete(drop=FALSE)
dev.off()

###

###

ggplot(ribodata, aes(patient_id, fill = istumor), color = 'black') + geom_bar(stat = 'count', position = 'dodge')
ggplot(fosjundata, aes(patient_id, fill = istumor), color = 'black') + geom_bar(stat = 'count', position = 'dodge')
ggplot(rsph1data, aes(patient_id, fill = istumor), color = 'black') + geom_bar(stat = 'count', position = 'dodge')

