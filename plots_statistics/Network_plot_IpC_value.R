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
library(plotrix)
library(car)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

choose_color <- function(tissue_type, mode='all') {
  if (tissue_type == 'Tumor') {color <- 'red'}
  if (tissue_type == 'AT1') {color <- 'green'}
  if (tissue_type == 'AT2') {color <- 'steelblue1'}
  if (tissue_type == 'Ciliated') {color <- 'yellow'}
  if (tissue_type == 'Club') {color <- 'pink'}
  color
}
choose_color <- Vectorize(choose_color)

choose_color2 <- function(tissue_type, mode='all') {
  if (tissue_type == 'Tumor') {color <- 'red'} #'#FF8C00'
  if (tissue_type == 'AT1') {color <- 'green'}
  if (tissue_type == 'AT2') {color <- 'steelblue1'}
  if (tissue_type == 'Ciliated') {color <- 'yellow'}
  if (tissue_type == 'Club') {color <- 'pink'}
  color
}
choose_color2 <- Vectorize(choose_color)

quant_edges <- 0.8
test_data_tumor <- fread(paste('../results/high_values_concat.csv', sep = ""))


testa<- test_data_tumor %>% group_by(cell_type_epi, patient_id) %>%
  dplyr::summarize('meanv' = mean(LRP))

tumor_specific_edges <- read.csv('./figures/tumor_edges.csv') %>%
  dplyr::mutate('tumor_interaction' = paste0(source_gene, target_gene)) %>%.$tumor_interaction
############
quantile_data <- test_data_tumor %>% group_by(sample_name) %>% do(data.frame(t(quantile(.$LRP, c(quant_edges)))))
colnames(quantile_data)[2] <- 'quant'
quantile_data <- test_data_tumor %>% left_join(quantile_data) %>% filter(LRP>quant)

####
ntissue <- quantile_data %>% dplyr::select(sample_name, cell_type_epi)%>%
  unique() %>%
  .$cell_type_epi %>%as.factor %>% summary()
ntissue
####
#compute IpC
cell_counts <- quantile_data %>% 
  dplyr::select(cell_type_epi, cell_id) %>%
  unique() %>%
  dplyr::group_by(cell_type_epi) %>%
  dplyr::summarize('cell_count' = n()) %>% 
  dplyr::select(cell_type_epi, cell_count)

quantile_data_IpC_all <- quantile_data %>%
  dplyr::group_by(source_gene, target_gene, cell_type_epi) %>%
  dplyr::summarize(InteractionCount = n()) %>%
    left_join(cell_counts) %>%
  dplyr::mutate(IpC = InteractionCount/cell_count)

highest_quantile_data <- quantile_data_IpC_all %>% group_by(cell_type_epi) %>% 
  dplyr::mutate('rank' = rank(1/IpC)) %>%
  dplyr::filter(rank<300)# %>% #300
  #dplyr::mutate(IpC = ifelse(cell_type_epi ==  'Tumor', 5*IpC, IpC))
summary(highest_quantile_data$cell_type_epi %>% as.factor())

highest_quantile_data$spec <- ifelse(paste0(highest_quantile_data$source_gene, highest_quantile_data$target_gene) %in% tumor_specific_edges,1,0)
####
graph_all_patients <- highest_quantile_data %>%
  dplyr::select(from =source_gene, to = target_gene,  tissue_type = cell_type_epi, weight = IpC, spec=spec) %>%
  graph_from_data_frame(directed=F)

E(graph_all_patients)$color <- choose_color(E(graph_all_patients)$tissue_type, 'all')
#E(graph_all_patients)[E(graph_all_patients)$weight < 1e-2]$color <-NA
E(graph_all_patients)$width <- E(graph_all_patients)$weight*10
#interesting_genes <- c('CXCL1','CXCL2','CXCL3', 'NFKBIA', 'TSPAN1', 'SRI', 'IGFB7', 'MANF', 'CALR', 'HSPA5',
#                       'CKS1B', 'CKS2', 'STMN1', 'PBX1', 'EP300', 'RBX1', 'PCM1', 'HNRNPA2B1', 'PACSIN2', 'ST13', 'SORBS2', 'XRCC6', 'CDC42EP1',
#                       'ID2', 'GADD45G', 'INTS6', 'ID4', 'HEXIM1', 'DDIT3', 'BRD2')

#interesting_genes <- c('CKS1B', 'CKS2', 'STMN1', 'RAD21')
#interesting_genes <- c('KLF4', 'AREG', 'PTGS2', 'ERRFI1', 'CITED2', 'GDF15', 'EPS8', 'TACSTD2')
#interesting_genes <- c('ETS2', 'NDRG1', 'BNIP3L', 'BNIP3', 'VEGFA', 'DDIT4', 'IGFBP3')
interesting_genes <- c('NNMT', 'MSN', 'AKRIC3')
#interesting_genes <- c('XRCC6', 'EP300', 'PACSIN2', 'CDC32EP1', 'ST13', 'RBX1')
#interesting_genes <- c('CLU', 'WIFI', 'PACSIN2', 'PTK7')
#interesting_genes <- c('JUN', 'FOS', 'HSPA8', 'HSPD1', 'HSP90AA1', 'DNAJB1')
#interesting_genes <- c('JUN', 'FOS', 'HSPA8', 'HSPD1', 'HSP90AA1', 'DNAJB1')
#interesting_genes <- c('NFKBIA', 'CXCL1', 'CXCL2')
#interesting_genes <- c('SYNE1', 'SYNE2', 'PCM1', 'SETD2', 'MACF1', 'AKAP9')
#interesting_genes <- c('PRDX1', 'NQO1', 'TALDO1', 'TXN', 'GDF15', 'MDM2', 'COX6CC', 'AKR1C3','CDKN1A')
#interesting_genes <- c('WIF1', 'CLU', 'PTK7')

E(graph_all_patients)$width

E(graph_all_patients)[E(graph_all_patients)$color == 'red']$width <- 3*E(graph_all_patients)[E(graph_all_patients)$color == 'red']$width

E(graph_all_patients)[E(graph_all_patients)$spec==1]$color <- 'white'
E(graph_all_patients)[E(graph_all_patients)$spec==1]$width <- 5* E(graph_all_patients)[E(graph_all_patients)$spec==1]$width
set.seed(0)
E(graph_all_patients)$color
E(graph_all_patients)$weight <- E(graph_all_patients)$weight*0.1
l <- layout_with_fr(graph_all_patients, niter=500)



# plot average graph # this also serves for the correct positining of cluster labels 
png('./figures/graph_all.png', width=4000, height=4000) #(50,50 for pdf)
par(bg='black')
#par(bg='white')

plot(graph_all_patients, vertex.size=0.1,
     edge.curved = curve_multiple(graph_all_patients, start = 0.1),
     edge.width = E(graph_all_patients)$width,
     vertex.label=ifelse((degree(graph_all_patients) > 1000) | (V(graph_all_patients)$name %in% interesting_genes), V(graph_all_patients)$name, NA),
     vertex.label.color = 'white',
     vertex.color=NA,
     vertex.label.cex= 5.0,
     layout=l)
rect(xleft=-0.93, ybottom=0.42, xright=-0.7, ytop=0.65, border='white', lwd=8)
text(x= -0.95, y = 0.43, 'a', cex = 6, col='white')
#draw.circle(-0.60, -0.04, 0.05, border='white', lwd= 8)
text(x= -0.3, y = 0.43, 'T3', cex = 6, col='white')
text(x= 0.15, y = 0.0, 'T3', cex = 6, col='white')


rect(xleft=-0.93, ybottom=-0.32, xright=-0.65, ytop=-0.08, border='white', lwd=8)
text(x= -0.95, y = -0.31, 'b', cex = 6, col='white')
rect(xleft=-0.6, ybottom=-0.85, xright=-0.25, ytop=-0.45, border='white', lwd=8)
text(x= -0.62, y = -0.84, 'c', cex = 6, col='white')

rect(xleft=-0.3, ybottom=-1.05, xright=0.05, ytop=-0.8, border='white', lwd=8)
text(x= -0.32, y = -1.04, 'd', cex = 6, col='white')
ellipse(c(-0.25,-0.75), shape= matrix(c(1,-0.5,-0.5,1), nrow=2),radius=0.3, col='white', lwd=6)
ellipse(c(-0.8,-0.2), shape= matrix(c(1,-0.1,-0.1,1), nrow=2),radius=0.15, col='white', lwd=6)
ellipse(c(-0.8,0.5), shape= matrix(c(1,0,0,1), nrow=2),radius=0.1, col='white', lwd=6)
ellipse(c(0.55,-0.55), shape= matrix(c(0.5,-0.3,-0.3,1.4), nrow=2),radius=0.25, col='white', lwd=3)

dev.off()

nodes <- delete_edges(graph_all_patients, E(graph_all_patients))
#edges <- paste0(highest_quantile_data$source_gene, highest_quantile_data$target_gene)
nodenames <- c(highest_quantile_data$source_gene, highest_quantile_data$target_gene) %>% unique()
V(nodes) %>% length()
###############################################################################
###############################################################################

patient_ids <- test_data_tumor$patient_id %>% unique()
patient_ids
cell_types <- test_data_tumor$cell_type_epi %>% unique()
font_size = 4


# loop for combined image of GRNs for different patients and tissues
# choose manuscriptof supplement patients below
png('./figures/graph_stratified_by_patient.png', width=2500, height=2500)#, res=100)
par(mfrow = c(6,6), mar = c(0,0,0,0))
par(bg = 'black')

plot(x=0:10, y=0:10)

for (cell_type in cell_types) {
  plot(x=0:10, y=0:10)
  text(x = 5,y = 2,cell_type, cex = 7, col = 'white')
}

manuscript_ids <- c(3,2,4,7,1)
supplement_ids <- c(5,6,8,9,10)
for (id in patient_ids[manuscript_ids]) {
print(id)
quantile_data_onepatient <- quantile_data %>% filter(patient_id == id)

quantile_data_onepatient %>% dplyr::select(cell_id, cell_type_epi) %>%
  unique() %>%
  .$cell_type_epi %>% as.factor() %>% summary()


plot(x=0:10, y=0:10)
text(x = 8,y = 5,id, cex = 7, col = 'white')

for (tumor_normal in cell_types) {
  quantile_data_one_patient_one_tissue <- quantile_data_onepatient %>%
    dplyr::filter(cell_type_epi == tumor_normal) %>%
    dplyr::filter(source_gene %in% nodenames, target_gene %in% nodenames)
  
  
  cell_counts <- quantile_data_one_patient_one_tissue %>% 
    dplyr::select(cell_type_epi, cell_id) %>%
    unique() %>%
    dplyr::group_by(cell_type_epi) %>%
    dplyr::summarize('cell_count' = n()) %>% 
    dplyr::select(cell_type_epi, cell_count)
  
  quantile_data_IpC_oneone <- quantile_data_one_patient_one_tissue %>%
    dplyr::filter(cell_type_epi == tumor_normal) %>%
    dplyr::group_by(source_gene, target_gene, cell_type_epi) %>%
    dplyr::summarize(InteractionCount = n()) %>%
    left_join(cell_counts) %>%
    dplyr::mutate(IpC = InteractionCount/cell_count)
  
  high_quantile_data_subtype <- quantile_data_IpC_oneone %>% 
    dplyr::mutate('rank' = rank(1/IpC)) %>%
    #dplyr::filter(rank<=100) %>%
    dplyr::select(from =source_gene, to = target_gene,  tissue_type = cell_type_epi,  weight = IpC) %>%
    dplyr::mutate(spec= ifelse(paste0(from, to) %in% tumor_specific_edges,1,0))
  
  
  graph_one_patient <- nodes
  
  for (i in seq(dim(high_quantile_data_subtype)[1])) {
    try({
      graph_one_patient <- graph_one_patient %>% add_edges(c(high_quantile_data_subtype$from[i], 
                                                             high_quantile_data_subtype$to[i]),
                                                             tissue_type = high_quantile_data_subtype$tissue_type[i],
                                                             weight = high_quantile_data_subtype$weight[i],
                                                           spec = high_quantile_data_subtype$spec[i])
    })
  }
  
  ### 
  if (dim(high_quantile_data_subtype)[1]>0) {
    print(dim(high_quantile_data_subtype)[1]>0)
  E(graph_one_patient)$color <- choose_color2(E(graph_one_patient)$tissue_type, 'both')
  E(graph_one_patient)[E(graph_one_patient)$weight < 1e-6]$color <-NA
  E(graph_one_patient)$width <- E(graph_one_patient)$weight*7
  E(graph_one_patient)[E(graph_one_patient)$spec==1]$color <- 'white'
  E(graph_one_patient)[E(graph_one_patient)$spec==1]$width <- E(graph_one_patient)[E(graph_one_patient)$spec==1]$width*2.0
  
  }
  
  if (dim(quantile_data_IpC_oneone)[1]>0 & quantile_data_IpC_oneone$cell_count[1]<5)   E(graph_one_patient)$width <- E(graph_one_patient)$weight*0

  plot(graph_one_patient, vertex.size=0.1, edge.width=E(graph_one_patient)$width, 
       edge.curved = curve_multiple(graph_one_patient, start = 0.0),
       vertex.label=ifelse( (degree(graph_one_patient) > 100) & (V(graph_one_patient)$name %in% c('CXCL3', 'NFKBIA', 'TSPAN1', 'SRI', 'IGFB7', 'MANF', 'CALR', 'HSPA5')), V(graph_one_patient)$name, NA),
       vertex.label.color = 'white',
       vertex.label.cex= 1.0,
       layout=l)
  
  if ((tumor_normal=='Tumor') &(id %in% c('p023','p027', 'p030', 'p019', 'p033', 'p034')) ) {
    #rect(xleft=-0.65, ybottom=0.65, xright=-0.25, ytop=1.0, border='white', lwd=2)
    text(x= -0.65, y = 0.9, 'T2', cex = font_size, col='white')
    draw.circle(-0.60, -0.04, 0.05, border='white', lwd= 2)
  }
    
  if (tumor_normal == 'Tumor' & (id %in% c('p032'))) {
    #rect(xleft=0.65, ybottom=0.1, xright=1.0, ytop=0.8, border='white', lwd=2)
    text(x= 0.75, y = 0.7, 'T1', cex = font_size, col='white')
    text(x= 0.4, y = 0.0, 'T6', cex = font_size, col='white')
    
  }
  
  if (tumor_normal == 'Tumor' & (id %in% c('p023','p031', 'p024', 'p019', 'p034'))) {
    #rect(xleft=0.65, ybottom=0.1, xright=1.0, ytop=0.8, border='white', lwd=2)
    text(x= -0.4, y = 0.55, 'T3', cex = font_size, col='white')
  }
  
  if (tumor_normal == 'Tumor' & (id %in% c('p031'))) {
    #rect(xleft=0.65, ybottom=0.1, xright=1.0, ytop=0.8, border='white', lwd=2)
    text(x= 0.2, y = 0.65, 'T4', cex = font_size, col='white')
  }
  
  if (tumor_normal == 'Tumor' & (id %in% c('p018'))) {
    #rect(xleft=0.65, ybottom=0.1, xright=1.0, ytop=0.8, border='white', lwd=2)
    text(x= 0.4, y = 0.0, 'T5', cex = font_size, col='white')
  }

  
  if (dim(quantile_data_IpC_oneone)[1]>0 & quantile_data_IpC_oneone$cell_count[1]>=5) { #((dim(quantile_data_IpC_oneone)[1]>0 & quantile_data_IpC_oneone$cell_count[1]>=5) & (tumor_normal != 'Tumor'))  {
  ellipse(c(-0.25,-0.75), shape= matrix(c(1,-0.5,-0.5,1), nrow=2),radius=0.3, col='white', lwd=2, center.pch=F)
  text(x= -0.7, y = -0.8, 'C3', cex = font_size, col='white')
  ellipse(c(-0.8,-0.2), shape= matrix(c(1,-0.1,-0.1,1), nrow=2),radius=0.15, col='white', lwd=2, center.pch=F)
  text(x= -0.9, y = 0.05, 'C4', cex = font_size, col='white')
  
  ellipse(c(-0.8,0.5), shape= matrix(c(1,0,0,1), nrow=2),radius=0.1, col='white', lwd=2, center.pch=F)
  text(x= -0.9, y = 0.75, 'C6', cex = font_size, col='white')
  }
  
  if (dim(quantile_data_IpC_oneone)[1]>0 & quantile_data_IpC_oneone$cell_count[1]>=5) {# ((dim(quantile_data_IpC_oneone)[1]>0 & quantile_data_IpC_oneone$cell_count[1]>=5) & (tumor_normal == 'Ciliated')) {
    ellipse(c(0.55,-0.55), shape= matrix(c(0.5,-0.3,-0.3,1.4), nrow=2),radius=0.25, col='white', lwd=2, center.pch=F)
    text(x= 0.85, y = -0.55, 'C8', cex = font_size, col='white') 
  }
  #if (tumor_normal == 'Tumor') rect(xleft=0.15, ybottom=-0.9, xright=0.5, ytop=-0.7, border='white', lwd=3)
  }
}
dev.off()
##########################################################
##########################################################

######################
#overlap between normal cells and tumor cells
#library(vecsets)
quantile_data$cell_type_epi %>% as.factor() %>% summary()
  cell_types <- quantile_data$cell_type_epi %>% unique
#interaction_data <- quantile_data %>% unite('interaction', c(source_gene, target_gene))
interaction_data <- quantile_data %>% unite('interaction', c(source_gene, target_gene)) %>%
  dplyr::group_by(interaction, cell_type_epi) %>%
  dplyr::summarize(counts = n())


compare_overlap <- function(cell_type) {
  tumor_data <- interaction_data %>% filter(cell_type_epi == 'Tumor') %>%
    dplyr::mutate(tumorcounts = counts)
  other_data <- interaction_data %>% filter(cell_type_epi == cell_type)  %>%
    dplyr::mutate(othercounts = counts)
  combined_data <- inner_join(tumor_data, other_data, by= 'interaction')
  combined_data$tumorcounts[is.na(combined_data$tumorcounts)] <- 0
  combined_data$othercounts[is.na(combined_data$othercounts)] <- 0

  
  print(cell_type)
  cor(combined_data$tumorcounts, combined_data$othercounts, method = 'spearman' )
  #combined_data
  
}

lapply(cell_types, compare_overlap)
  
