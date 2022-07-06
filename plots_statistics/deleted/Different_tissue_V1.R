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
  if (tissue_type == 'Neuroendocrine') {color <- 'orange'}
  if (tissue_type == 'Club') {color <- 'pink'}
  color
}
choose_color <- Vectorize(choose_color)

quant_edges <- 0.8
#test_data_tumor <- fread(paste('~/Projekte/singlecell/data/high_values_concat.csv', sep = ""))
test_data_tumor <- fread(paste('/mnt/scratch2/mlprot/Projekte/singlecell/data/high_values_concat.csv', sep = ""))

filtered_test_data <- test_data_tumor %>% filter(target_gene<source_gene) %>% dplyr::select(-c(V1, cell_id)) %>%
  #mutate('patient_id' = stringr::str_split(sample_name, '_', simplify = T)[,1])
  dplyr::mutate('patient_id' = substr(sample_name, 1,4)) %>%
  dplyr::filter(patient_id!= 'p029', patient_id != 'p028') %>% # check below for cell count (no tumor cells?)
  filter(source_gene!='meanpercell', target_gene != 'meanpercell') %>%
  
  filter(cell_type_epi != 'Neuroendocrine')

quantile_now_all_patients <- filtered_test_data %>% group_by(sample_name) %>% do(data.frame(t(quantile(.$LRP, c(quant_edges)))))
colnames(quantile_now_all_patients)[2] <- 'quant'
quantile_filtered_test_data_all_patients <- filtered_test_data %>% left_join(quantile_now_all_patients) %>% filter(LRP>quant)

###############################################
#balance tissue types without regard to patient
###############################################
ntissue <- quantile_filtered_test_data_all_patients %>% dplyr::select(sample_name, cell_type_epi)%>%
  unique() %>%
  .$cell_type_epi %>%as.factor %>% summary()
ntissue

cells <- quantile_filtered_test_data_all_patients %>% dplyr::select(sample_name, cell_type_epi)%>%
  unique()


return_sample <- function(tissue_type, cells, ncells_per_tissue, overweight_tumor) {
  cells_one_tissue <- cells %>% filter(cell_type_epi == tissue_type)
  old_ncells <- dim(cells_one_tissue)[1]
  ratio = ncells_per_tissue/old_ncells
  if ((tissue_type == 'Tumor')& overweight_tumor) ncells_per_tissue <- 5*ncells_per_tissue

  if (ratio <=1) {
    new_data = cells_one_tissue[sample(ncells_per_tissue,as.integer(ncells_per_tissue), replace=F),]
  }
  else {
    if (dim(cells_one_tissue)[1]<1){
      new_data <- cells_one_tissue
    }
    else {
      new_data = cells_one_tissue[sample(old_ncells,as.integer(ncells_per_tissue), replace=T),]
    }
  }
  new_data
}

resample_cells <- function(ncells_per_tissue, cells, overweight_tumor=T) {
  ntissue <- cells$cell_type_epi %>%as.factor %>% summary()
  tissues <- cells$cell_type_epi %>% unique()
  print(tissues)
  print(ntissue)
  resampled_cells <- lapply(tissues, function(tissue_type) {return_sample(tissue_type, cells, ncells_per_tissue, overweight_tumor=overweight_tumor)}) %>%
    rbindlist()
  resampled_cells
}

#return_sample('Neuroendocrine', cells, ncells_per_tissue)

ncells_per_tissue <-20
resampled_cells <- resample_cells(ncells_per_tissue, cells)

resampled_data <- resampled_cells %>% left_join(quantile_filtered_test_data_all_patients)

########
#correct for dropouts

expr_data <- fread('/mnt/scratch2/mlprot/Projekte/singlecell/data/epi_top2000.csv')[,4:2003]
cell_type_epi <- fread('/mnt/scratch2/mlprot/Projekte/singlecell/data/epi_top2000.csv')$cell_type_epi
#zerorate_source <- data.frame('source_gene' = colnames(expr_data), 'zerorate_source' = colMeans(expr_data==0))
#zerorate_target <- data.frame('target_gene' = colnames(expr_data), 'zerorate_target' = colMeans(expr_data==0))

cell_frame <- data.frame('cell_type_epi' = cell_type_epi, expr_data)

zerorate_source <-  aggregate(cell_frame[,-(1)]==0, list(cell_frame$cell_type_epi), mean) %>%
  pivot_longer(!Group.1, names_to ='source_gene', values_to = 'zerorate_source')

colnames(zerorate_source)[1] <- 'cell_type_epi'
zerorate_target <-zerorate_source
colnames(zerorate_target)[2:3] <- c('target_gene', 'zerorate_target') 

resampled_data <- left_join(resampled_data, zerorate_source) %>%
  left_join(zerorate_target) %>%
  mutate(random = runif(dim(.)[1]), zerorate = (zerorate_target+zerorate_source)) %>%
  mutate(selected_rows = zerorate>(random-0.01))
resampled_data <- resampled_data %>% filter(selected_rows)

########
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

E(graph_all_patients_resampled)$weight <- 0.00001
E(graph_all_patients_resampled)$color[all_edges$countvs<4] <- NA
set.seed(1)
l <- layout_with_fr(graph_all_patients_resampled, niter=500)
#png('~/Projekte/singlecell/figures/singlecell_single_all_patients_all_tissue_resampled.png', width=6000, height=6000, res=100)
png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_single_all_patients_all_tissue_resampled.png', width=2000, height=2000, res=50)
par(bg='black')
plot(graph_all_patients_resampled, vertex.size=0.1, edge.width=0.2, 
     edge.curved = curve_multiple(graph_all_patients_resampled, start = 0.04),
     vertex.label=ifelse( (degree(graph_all_patients_resampled) > 200) | (V(graph_all_patients_resampled)$name %in% c('RSPH1','CXCL3','HSP90AA1','HSP90AB1', 'CXCL1', 'CXCL2','ZFP36','NFKBIA', 'TSPAN1', 'SRI', 'IGFB7', 'MANF', 'CALR', 'HSPA5')), V(graph_all_patients_resampled)$name, NA),
     vertex.label.color = 'white',
     vertex.color=NA,
     vertex.label.cex= 2,
     layout=l,
     xlim=c(-0.4,0.0), ylim=c(-0.4,-0.1))
dev.off()

nodes <- delete_edges(graph_all_patients_resampled, E(graph_all_patients_resampled))
V(nodes) %>% length()
######################
######################
#quantify interactions

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
quantification_junfosal_all <- quantification_all %>% dplyr::filter((source_gene %in% junfos)&  (target_gene %in% junfos))
quantification_junfos <- quantification %>% dplyr::filter((source_gene %in% junfos)&  (target_gene %in% junfos))

####
####
#quantify ID2-ID4 network
ID2ID4 <- c('ID2', 'ID4', 'IER3', 'HEXIM1', 'INTS6', 'BRD2', 'GADD45G', 'KLF4')

quantificationID2ID4_all <- quantification_all %>% dplyr::filter((source_gene %in% ID2ID4) &  (target_gene %in% ID2ID4))
quantificationID2ID4 <- quantification%>% dplyr::filter((source_gene %in% ID2ID4) &  (target_gene %in% ID2ID4))
all_interaction_count_ID2ID3<- quantile_filtered_test_data_all_patients %>% 
  group_by(cell_type_epi) %>%
  dplyr::filter((source_gene %in% ID2ID4) &  (target_gene %in% ID2ID4)) %>%
  dplyr::summarize('interactioncount' = n()) %>%
  left_join(cell_counts) %>%
  dplyr::mutate('relativecount' = interactioncount/cell_count)


####
#quantify ID2-ID4 network
RSPH1 <-c('RSPH1', 'SPA17', 'ZMYND10', 'SOX2', 'HSP90AA1', 'GSTA1', 'SSBP4', 'RUVBL1', 'RUVBL2', 'IGFBP7', 'SLC44A4', 'TSPAN1', 'MAPK10', 'MLF1', 'MUC16', 'IGFBP7','MAPK10', 'AKAP9', 'IGFBP5')

quantificationRSPH1 <- quantification_all %>% dplyr::filter((source_gene %in% RSPH1) &  (target_gene %in% RSPH1))

all_interaction_count_RSPH1<- quantile_filtered_test_data_all_patients %>% 
  group_by(cell_type_epi) %>%
  dplyr::filter((source_gene %in% RSPH1) &  (target_gene %in% RSPH1)) %>%
  dplyr::summarize('interactioncount' = n()) %>%
  left_join(cell_counts) %>%
  dplyr::mutate('relativecount' = interactioncount/cell_count)

####
cxcl <- c('CXCL1', 'CXCL2', 'CXCL3', 'NFKBIA')

quantification_cxcl <- quantification_all %>% dplyr::filter((source_gene %in% cxcl)&  (target_gene %in% cxcl))

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


###############################################################################
##############################################################################
# plot Heatmap with network expression
patient_mutations <- data.frame(patient_id = c('p018', 'p019','p023','p024','p027','p030','p031','p032','p033','p034'),
                                'KRAS' = c(1,0,1,0,0,1,1,1,1,0),
                                'TP53' = c(0,0,1,0,1,0,0,0,0,0),
                                'PIK3CA' = c(0,0,0,0,0,0,1,0,0,0),
                                'EML4' = c(0,0,0,1,0,0,0,0,0,0))

p_mat <- patient_mutations %>% dplyr::select(-patient_id) %>% as.matrix()
rownames(p_mat) <- patient_mutations$patient_id
mutant_heatmap <- p_mat %>% Heatmap()

quantification_heatmap <- quantile_filtered_test_data_all_patients %>% 
  filter(cell_type_epi == 'Tumor') %>%
  dplyr::select(source_gene, target_gene, patient_id) %>%
  dplyr::mutate('junfos_int' = ifelse((source_gene %in% junfos) & (target_gene %in% junfos),1,0)) %>%
  dplyr::mutate('RSPH1_int' = ifelse((source_gene %in% RSPH1) & (target_gene %in% RSPH1),1,0)) %>%
  dplyr::mutate('ID2ID4_int' = ifelse((source_gene %in% ID2ID4) & (target_gene %in% ID2ID4),1,0)) %>%
  dplyr::mutate('cxcl_int' = ifelse((source_gene %in% cxcl) & (target_gene %in% cxcl),1,0)) %>%
  dplyr::select(-c(source_gene, target_gene))

ncells_perid <- quantile_filtered_test_data_all_patients %>% 
  dplyr::select(patient_id, sample_name) %>%
  unique() %>%
  group_by(patient_id) %>%
  dplyr::summarize('cellcount' = n())
  

summarized_quant_hm <- quantification_heatmap %>% group_by(patient_id) %>%
                    dplyr::summarize('meanjunfos_int' = sum(junfos_int),
                      'meanRSPH1_int' = sum(RSPH1_int),
                      'meanID2ID4_int' = sum(ID2ID4_int),
                      'meancxcl_int' = sum(cxcl_int)) %>%
  left_join(ncells_perid) %>%
  dplyr::mutate(meanjunfos_int = meanjunfos_int/cellcount, 
                meanRSPH1_int = meanRSPH1_int/cellcount, 
                meanID2ID4_int = meanID2ID4_int/cellcount, 
                meancxcl_int = meancxcl_int/cellcount)

library(ComplexHeatmap)
heatmap_matrix <- summarized_quant_hm %>%dplyr::select(-c(patient_id, cellcount)) %>% as.matrix() 
rownames(heatmap_matrix) <- summarized_quant_hm$patient_id
heatmap_matrix <- apply(heatmap_matrix, 2, function(x) (x-min(x))/max(x-min(x)))

heatmap_matrix <- log(1+heatmap_matrix)

png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/patient_heatmap.png', width=1000, height=1000)#, res=100)
Heatmap(heatmap_matrix) + mutant_heatmap
dev.off()
###############################################################################
###############################################################################
#look for different interactions between KRAS and ~KRAS
cells_per_patient <- quantile_filtered_test_data_all_patients %>% 
  dplyr::select(patient_id, sample_name) %>%
  group_by(patient_id) %>%
  dplyr::summarize('patient_cell_count' = n())
  
quantile_data_with_mutations <- quantile_filtered_test_data_all_patients %>% 
  left_join(patient_mutations) %>%
  dplyr::mutate('KRASm' = ifelse(KRAS==1, 'y','n')) %>%
 group_by(KRASm, source_gene, target_gene) %>% 
  dplyr::select(source_gene,target_gene, patient_id, KRASm) %>%
  dplyr::mutate('counts' = n()) %>%
  left_join(cells_per_patient) %>%
  dplyr::mutate('norm_count' = counts / patient_cell_count) %>%
  unique()
###
#compare to expression data

expr_data <- fread('/mnt/scratch2/mlprot/Projekte/singlecell/data/epi_top2000.csv')[,4:2003] %>%
  dplyr::select(colnames(.)[colnames(.) %in% cxcl])

cell_type_epi <- fread('/mnt/scratch2/mlprot/Projekte/singlecell/data/epi_top2000.csv') %>%
  dplyr::select(cell_type_epi, patient_id)
#zerorate_source <- data.frame('source_gene' = colnames(expr_data), 'zerorate_source' = colMeans(expr_data==0))
#zerorate_target <- data.frame('target_gene' = colnames(expr_data), 'zerorate_target' = colMeans(expr_data==0))

cell_frame <- cbind(cell_type_epi, expr_data) %>% left_join(patient_mutations) %>% filter(cell_type_epi=='Tumor')
x <- cell_frame %>%filter(KRAS==1) %>%.$CXCL3
y <- cell_frame %>%filter(KRAS==0) %>%.$CXCL3

wilcox.test(x,y)

ggplot(cell_frame %>% filter(CXCL3!=0), aes(y=NFKBIA, fill = KRAS)) + geom_boxplot()

###############################################################################
###############################################################################
patient_ids <- filtered_test_data$patient_id %>% unique()
png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_individual_patitents_resampled.png', width=2000, height=2000)#, res=100)
par(bg = 'black')
par(mfrow = c(5,2))
for (id in patient_ids[1:10]) {
print(id)
quantile_filtered_test_data_onepatient <- quantile_filtered_test_data_all_patients %>% filter(patient_id == id)

####
#resample patient
####
cells <- quantile_filtered_test_data_onepatient %>% dplyr::select(sample_name, cell_type_epi)%>%
  unique()

ncells_per_tissue <- 20
resampled_cells <- resample_cells(ncells_per_tissue, cells)

resampled_data_onepatient <- resampled_cells %>% left_join(quantile_filtered_test_data_onepatient)

########
#correct for dropouts

#expr_data <- fread('/mnt/scratch2/mlprot/Projekte/singlecell/data/epi_top2000.csv')[,3:2002]
#cell_type_epi <- fread('/mnt/scratch2/mlprot/Projekte/singlecell/data/epi_top2000.csv')$cell_type_epi

#cell_frame <- data.frame('cell_type_epi' = cell_type_epi, expr_data)

#zerorate_source <-  aggregate(cell_frame[,-(1)]==0, list(cell_frame$cell_type_epi), mean) %>%
#  pivot_longer(!Group.1, names_to ='source_gene', values_to = 'zerorate_source')
#
#colnames(zerorate_source)[1] <- 'cell_type_epi'
#zerorate_target <-zerorate_source
#colnames(zerorate_target)[2:3] <- c('target_gene', 'zerorate_target') 
#
#resampled_data_onepatient <- left_join(resampled_data_onepatient, zerorate_source) %>%
#  left_join(zerorate_target) %>%
#  mutate(random = runif(dim(.)[1]), zerorate = zerorate_target * zerorate_source) %>%
#  mutate(selected_rows = zerorate>(random-0.01))
#resampled_data_onepatient <- resampled_data_onepatient %>% filter(selected_rows)

########

####

frame_one_patient_resampled <- resampled_data_onepatient %>%
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

E(graph_one_patient_resampled)$color[one_patient_edges$countvs<3] <- NA

#png(paste('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_single_cells_one_patient_colored_',id, '.png', sep=''), width=2000, height=2000, res=100)
#par(bg='black')

plot(graph_one_patient_resampled, vertex.size=0.1, edge.width=0.1, 
     edge.curved = curve_multiple(graph_one_patient_resampled, start = 0.01),
     vertex.label=ifelse( (degree(graph_one_patient_resampled) > 100) | (V(graph_one_patient_resampled)$name %in% c('CXCL3', 'NFKBIA', 'TSPAN1', 'SRI', 'IGFB7', 'MANF', 'CALR', 'HSPA5')), V(graph_one_patient_resampled)$name, NA),
     vertex.label.color = 'white',
     vertex.label.cex= 1.0,
     layout=l,
     xlim=c(-0.4,0.0), ylim=c(-0.4,-0.1))

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

ncells_per_tissue <- 10000
resampled_cells <- resample_cells(ncells_per_tissue, cells, overweight_tumor=F)

resampled_data <- resampled_cells %>% left_join(quantile_filtered_test_data_all_patients)

resampled_data$cell_type_epi %>% as.factor() %>% summary()
cell_types <- resampled_data$cell_type_epi %>% unique
interaction_data <- resampled_data %>% unite('interaction', c(source_gene, target_gene))

compare_overlap <- function(cell_type) {
  tumor_data <- interaction_data %>% filter(cell_type_epi == 'Tumor') %>% .$interaction
  print(tumor_data %>% dim)
  other_data <- interaction_data %>% filter(cell_type_epi == cell_type) %>% .$interaction
  c(cell_type, as.numeric(length(other_data[other_data %in%intersect(tumor_data, other_data)]))/length(other_data))
  
}

lapply(cell_types, compare_overlap)
