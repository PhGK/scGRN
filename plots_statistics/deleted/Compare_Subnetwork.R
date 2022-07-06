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
#
choose_color <- function(tissue_type, mode='all') {
  if (tissue_type == 'Tumor') {color <- 'red'}
  if (tissue_type == 'AT1') {color <- 'green'}
  if (tissue_type == 'AT2') {color <- 'steelblue1'}
  if (tissue_type == 'Ciliated') {color <- 'darkgoldenrod2'}
  if (tissue_type == 'Club') {color <- 'pink'}
  color
}
choose_color <- Vectorize(choose_color)

quant_edges <- 0.8
test_data_tumor <- fread(paste('../results/high_values_concat.csv', sep = ""))


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
######################

all_graph <- quantile_data %>%
  dplyr::select(source_gene, target_gene) %>% 
  dplyr::group_by(source_gene, target_gene) %>%
  dplyr::summarize(weight = n()) %>%
  dplyr::filter(weight>500) %>%
  dplyr::select(-weight) %>%
  as.matrix() %>%
  graph_from_edgelist() %>%
  as.undirected()

all_graph <- quantile_data %>%
  group_by(cell_type_epi, patient_id) %>%
  dplyr::mutate(IpC = 1/(n()+10)) %>%
  dplyr::group_by(source_gene, target_gene) %>%
  dplyr::summarize(IpC = sum(IpC)) %>%
  dplyr::filter(IpC>1e-2) %>%
  dplyr::select(-IpC) %>%
  as.matrix() %>%
  graph_from_edgelist() %>%
  as.undirected()

#E(all_graph)$weight <-1
#all_graph <- simplify(all_graph, edge.attr.comb =list(weight='sum'))
#all_graph <- delete_edges(all_graph,E(all_graph)[E(all_graph)$weight<20])
E(all_graph)$weight <-1
#E(all_graph)$weight
#########
# number of edges berücksichtigen
#########
#######
set.seed(0)
modules <- cluster_louvain(all_graph)
print(modules)
modules$membership

degree(all_graph, V(all_graph))
degreeframe <- data.frame(name = V(all_graph)$name, degree = degree(all_graph, V(all_graph)), cluster = modules$membership) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(N = n()) %>%
  dplyr::ungroup() #%>%
  #dplyr::filter(N!=max(N))
cluster_sizes <-degreeframe %>% group_by(cluster) %>% dplyr::summarize(N = n()) %>% filter(N>1)
###


###
patient_ids <- quantile_data %>% dplyr::select(patient_id) %>% unique() 



################################################
##### long format to normalize per cluster
################################################

###
patient_ids <- quantile_data %>% dplyr::select(patient_id) %>% unique() 
get_relative_counts <- function(idx, cell_type) {
  
  gene_set <- degreeframe %>% filter(degree>-1, cluster==idx) %>% .$name
  
  ncells_perid <- quantile_data %>% 
    filter(cell_type_epi == cell_type) %>%
    dplyr::select(patient_id, sample_name) %>%
    unique() %>%
    group_by(patient_id) %>%
    dplyr::summarize('all_cellcount' = n())
  
  count_data <- quantile_data %>% 
    filter(cell_type_epi == cell_type) %>%
    filter(source_gene %in% gene_set, target_gene %in% gene_set) %>%
    dplyr::select(patient_id, sample_name) %>%
    group_by(patient_id) %>%
    dplyr::summarize('spec_cell_count' = n())
  
  intermediate <- left_join(ncells_perid, count_data) %>% 
    dplyr::mutate('ratio' = spec_cell_count/all_cellcount)
  
  data.frame('patient_id' = intermediate$patient_id, 'clustervalue' = intermediate$ratio, 'cluster' = idx)
}


get_all_clusters <- function(cell_type) {
  all_clusters <- rbindlist(lapply(seq(9), function(x) get_relative_counts(x, cell_type)))
  all_clusters$cell_type_epi <- cell_type
  all_clusters
}

get_all_clusters('Tumor')
get_relative_counts(6, 'Tumor')

cell_types <- quantile_data %>% .$cell_type_epi %>% unique()

all_data <- rbindlist(lapply(cell_types, get_all_clusters))
normalized_data <- all_data %>% dplyr::group_by(cluster) %>%
  dplyr::mutate(maxv = max(clustervalue, na.rm=T), minv = min(clustervalue, na.rm=T)) %>%
  dplyr::mutate(normv = (clustervalue)/(maxv))

patient_mutations <- data.frame(patient_id = c('p018', 'p019','p023','p024','p027','p030','p031','p032','p033','p034'),
                                'KRAS' = c(1,0,1,0,0,1,1,1,1,0),
                                'TP53' = c(0,0,1,0,1,0,0,0,0,0),
                                'PIK3CA' = c(0,0,0,0,0,0,1,0,0,0),
                                'EML4' = c(0,0,0,1,0,0,0,0,0,0))

patient_histology <- data.frame(patient_id = c('p018', 'p019','p023','p024','p027','p030','p031','p032','p033','p034'),
                                Histology = as.factor(c('acinar', 'acinar', 'solid', 'acinar', 'solid', 'mucinuous', 'acinar', 'lepidic', 'papillary', 'sarcomatoid')))

p_mat <- patient_mutations %>% dplyr::select(-patient_id) %>% as.matrix()
rownames(p_mat) <- patient_mutations$patient_id

ht_opt(legend_border = "black",
       heatmap_border = TRUE)
col_fun = colorRamp2(c(0, 1), c("white", "red"))
col_fun2 = colorRamp2(c(0, 1), c("white", "black"))
w <- unit(4, "cm")
h <- unit(6, "cm")

mutant_heatmap <- p_mat %>% 
  Heatmap(cluster_columns = F, width=0.5*w, height=h, name = 'Mutation', col = col_fun2, 
          show_heatmap_legend = F, column_names_gp = grid::gpar(fontsize = 15), column_title = 'Mutation')

h_mat <- patient_histology %>% dplyr::select(-patient_id) %>% as.matrix()
rownames(h_mat) <- patient_histology$patient_id
histology_heatmap <- h_mat %>% Heatmap(cluster_columns = F, width=0.3*w, height=h, name = 'Histology', col=(rainbow(14)[c(1,3,4,6,8,10)]), row_names_gp = gpar(fontsize=15),
                              column_names_gp = gpar(fontsize=15), heatmap_legend_param =  list(title_gp =  gpar(fontsize=13), labels_gp = gpar(fontsize = 13)))


prepare_matrix <- function(cell_type) {
  all_clusters_per_celltype <- normalized_data %>% dplyr::filter(cell_type_epi == cell_type) %>%
    dplyr::select(patient_id, cell_type_epi, normv, -cell_type_epi) %>%
    dplyr::mutate(cluster = paste0('C', cluster)) %>%
    pivot_wider(names_from = cluster, values_from = normv)
  
  all_clusters_per_celltype[is.na(all_clusters_per_celltype)] <- 0
  
  all_clusters_per_celltype <- all_clusters_per_celltype %>%
    merge(patient_ids, all=T)
  cell_type_matrix <- all_clusters_per_celltype %>% 
    dplyr::select(-patient_id) %>%
    as.matrix()
  rownames(cell_type_matrix) <- all_clusters_per_celltype$patient_id
  cell_type_matrix

}
prepare_matrix('AT1')


png('./figures/patient_heatmap.png', width=1000, height=250)
Heatmap(prepare_matrix('Tumor'), col = col_fun, cluster_columns = F, column_title ='Tumor', name = 'IpC (scaled)', column_names_gp = gpar(fontsize=15),
        heatmap_legend_param =  list(title_gp =  gpar(fontsize=13), labels_gp = gpar(fontsize = 13))) + 
  Heatmap(prepare_matrix('Club'), col = col_fun, cluster_columns = F, column_title ='Club',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15))+  
  Heatmap(prepare_matrix('AT1'), col = col_fun, cluster_columns = F, column_title ='AT1',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15)) + 
  Heatmap(prepare_matrix('AT2'), col = col_fun, cluster_columns = F, column_title ='AT2',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15)) + 
  Heatmap(prepare_matrix('Ciliated'), col = col_fun, cluster_columns = F, column_title ='Ciliated',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15)) + 
  mutant_heatmap+
  histology_heatmap
dev.off()
#############################
###visualize general networks
#############################

plot_general_clusters <- function(module_id) {
  df <- degreeframe %>% dplyr::filter(cluster == module_id) 

  spec_data <- quantile_data %>% dplyr::filter(source_gene %in% df$name, target_gene %in% df$name)
  
  cell_count <- quantile_data %>%
    dplyr::select(cell_id, cell_type_epi) %>%
    unique() %>%
    dplyr::group_by(cell_type_epi) %>%
    dplyr::summarize(N = n())
  
  general_specific_subedges <-  spec_data %>%
    dplyr::group_by(source_gene, target_gene, cell_type_epi) %>%
    dplyr::summarize(interactioncount = n()) %>%
    left_join(cell_count) %>%
    group_by(source_gene, target_gene) %>%
    ungroup() %>%
    dplyr::mutate(interactioncount = interactioncount/N) %>%
    dplyr::select(from =source_gene, to = target_gene,  tissue_type = cell_type_epi, weight = interactioncount)
  
  general_specific_subgraph <- graph_from_data_frame(general_specific_subedges, directed=F)
  #tumor_specific_graph <- graph_from_edgelist(tumor_specific_edges %>% dplyr::select(from, to), directed=F)
  
  E(general_specific_subgraph)$color <- choose_color(E(general_specific_subgraph)$tissue_type)
  E(general_specific_subgraph)$color
  E(general_specific_subgraph)$weight <- 0.01* E(general_specific_subgraph)$weight
  E(general_specific_subgraph)$width <- E(general_specific_subgraph)$weight*2000
  set.seed(0)
  l <- layout_with_fr(general_specific_subgraph, niter = 200)
  
  plot(general_specific_subgraph, vertex.size=0.2,
       edge.curved = curve_multiple(general_specific_subgraph, start = 0.05),
       edge.width = E(general_specific_subgraph)$width,
       vertex.label.color = 'black',
       vertex.color=NA,
       vertex.label = ifelse((degree(general_specific_subgraph) < 130) , V(general_specific_subgraph)$name, NA),
       vertex.label.cex= 3.5,
       layout=l)
  text(-1.2,1.0, labels = paste0('C', module_id), cex=5)
}
clusters <- degreeframe$cluster %>% unique()
#par(mfrow = c(1,length(clusters)))
#dev.off()
layout.matrix <- matrix(c(1,1,2,2,3,3,1,1,2,2,3,3,4,5,6,7,8,9), nrow= 6)
#png('./figures/general_specific_graph', width=600*3, height = 500*3 )
#par(bg='white', mfrow = c(3,3), mar= c(0,0,0,0))
#for (i in clusters) {
#  plot_general_clusters(i)
#}
#dev.off()

png('./figures/general_specific_graph', width=600*3, height = 500*6)
layout(mat = layout.matrix)
par(mar=c(0,0,0,0))
for (i in c(1,3,6,2,4,5,7,8,9)) {
  plot_general_clusters(i)
}
dev.off()
##############################  
## cancer networks
##############################

N1 <- c('ID2', 'INTS6', 'GADD45G', 'DDIT3', 'HEXIM1', 'KLF4')
N2 <- c('EP300', 'RBX1', 'PCM1', 'ST13', 'PACSIN2', 'XRCC6', 'CDC42EP1')
N3 <- c('STMN1', 'CKS2', 'CKS1B', 'HMGB1')
N4 <- c('HSP90B1', 'HSPA5', 'CALR', 'MANF')
own_clusters <- list(N1,N2,N3)

get_relative_counts <- function(idx, cell_type) {
  
  gene_set <- own_clusters[[idx]]
  
  ncells_perid <- quantile_data %>% 
    filter(cell_type_epi == cell_type) %>%
    dplyr::select(patient_id, sample_name) %>%
    unique() %>%
    group_by(patient_id) %>%
    dplyr::summarize('all_cellcount' = n())
  
  count_data <- quantile_data %>% 
    filter(cell_type_epi == cell_type) %>%
    filter(source_gene %in% gene_set, target_gene %in% gene_set) %>%
    dplyr::select(patient_id, sample_name) %>%
    group_by(patient_id) %>%
    dplyr::summarize('spec_cell_count' = n())
  
  intermediate <- left_join(ncells_perid, count_data) %>% 
    dplyr::mutate('ratio' = spec_cell_count/all_cellcount)
  
  data.frame('patient_id' = intermediate$patient_id, 'clustervalue' = intermediate$ratio, 'cluster' = idx)
}


get_all_clusters <- function(cell_type) {
  all_clusters <- rbindlist(lapply(seq(3), function(x) get_relative_counts(x, cell_type)))
  all_clusters$cell_type_epi <- cell_type
  all_clusters
}

get_all_clusters('Tumor')
get_relative_counts(3, 'Tumor')

cell_types <- quantile_data %>% .$cell_type_epi %>% unique()


all_data <- rbindlist(lapply(cell_types, get_all_clusters))
normalized_data <- all_data %>% dplyr::group_by(cluster) %>%
  dplyr::mutate(maxv = max(clustervalue, na.rm=T), minv = min(clustervalue, na.rm=T)) %>%
  dplyr::mutate(normv = (clustervalue)/(maxv))

prepare_matrix <- function(cell_type) {
  all_clusters_per_celltype <- normalized_data %>% dplyr::filter(cell_type_epi == cell_type) %>%
    dplyr::select(patient_id, cell_type_epi, normv, -cell_type_epi) %>%
    dplyr::mutate(cluster = paste0('T', cluster)) %>%
    pivot_wider(names_from = cluster, values_from = normv)
  
  all_clusters_per_celltype[is.na(all_clusters_per_celltype)] <- 0
  
  all_clusters_per_celltype <- all_clusters_per_celltype %>%
    merge(patient_ids, all=T)
  cell_type_matrix <- all_clusters_per_celltype %>% 
    dplyr::select(-patient_id) %>%
    as.matrix()
  rownames(cell_type_matrix) <- all_clusters_per_celltype$patient_id
  cell_type_matrix
  
}
prepare_matrix('AT1')

ht_opt(legend_border = "black",
       heatmap_border = TRUE)
col_fun = colorRamp2(c(0, 1), c("white", "red"))
col_fun2 = colorRamp2(c(0, 1), c("white", "black"))
#Heatmap(prepare_matrix('Tumor'))
w <- unit(4, "cm")
h <- unit(4, "cm")

mutant_heatmap <- p_mat %>% 
  Heatmap(cluster_columns = F, width=0.5*w, height=h, name = 'Mutation', col = col_fun2, 
          show_heatmap_legend = F, column_names_gp = grid::gpar(fontsize = 15), column_title = 'Mutation')

histology_heatmap <- h_mat %>% Heatmap(cluster_columns = F, width=0.3*w, height=h, name = 'Histology', col=(rainbow(14)[c(1,3,4,6,8,10)]), row_names_gp = gpar(fontsize=15),
                                       column_names_gp = gpar(fontsize=15), heatmap_legend_param =  list(title_gp =  gpar(fontsize=15), labels_gp = gpar(fontsize = 15)))


png('./figures/patient_heatmap_own_clusters.png', width=900, height=400)
par(mar = c(2,2,2,2), oma = c(2,2,2,2))
Heatmap(prepare_matrix('Tumor'), col = col_fun, cluster_columns = F,column_names_rot = 0,  column_title ='Tumor', name = 'IpC (scaled)', width =w, 
        column_names_gp = gpar(fontsize=15), heatmap_legend_param =  list(title_gp =  gpar(fontsize=13), labels_gp = gpar(fontsize = 13))) + 
  Heatmap(prepare_matrix('Club'), col = col_fun, cluster_columns = F,column_names_rot = 0, column_title ='Club',show_heatmap_legend = F,width =w, column_names_gp = gpar(fontsize=15))+  
  Heatmap(prepare_matrix('AT1'), col = col_fun, cluster_columns = F,column_names_rot = 0, column_title ='AT1',show_heatmap_legend = F,width =w, column_names_gp = gpar(fontsize=15)) + 
  Heatmap(prepare_matrix('AT2'), col = col_fun, cluster_columns = F, column_names_rot = 0, column_title ='AT2',show_heatmap_legend = F, width =w, column_names_gp = gpar(fontsize=15)) + 
  Heatmap(prepare_matrix('Ciliated'), col = col_fun, cluster_columns = F, column_names_rot = 0, column_title ='Ciliated',show_heatmap_legend = F, width =w, column_names_gp = gpar(fontsize=15)) + 
  mutant_heatmap+
  histology_heatmap
dev.off()
#######
#quantify the four networks from above
own_network <- function(network) {
  ncells <- quantile_data %>%
    dplyr::select(cell_type_epi, cell_id) %>%
    unique() %>%
    group_by(cell_type_epi) %>%
    dplyr::summarize('ncells' = n())
    
    
  data <- quantile_data %>%
    dplyr::filter(source_gene %in% network, target_gene %in% network) %>%
    dplyr::group_by(cell_type_epi) %>%
    dplyr::summarize('Interactions' = n()) %>%
    ungroup() %>%
    left_join(ncells) %>%
    dplyr::mutate('IpC' = Interactions/ncells)
  data
  }

own_network(N1)
alldata <- rbind(cbind(own_network(N1), 'network' = 'N1'), cbind(own_network(N2),'network' =  'N2'), 
                 cbind(own_network(N3),'network' =  'N3'))

#######
#visualize the four networks from above

choose_color <- function(tissue_type, mode='all') {
  if (tissue_type == 'Tumor') {color <- 'red'}
  if (tissue_type == 'AT1') {color <- 'green'}
  if (tissue_type == 'AT2') {color <- 'steelblue1'}
  if (tissue_type == 'Ciliated') {color <- 'yellow'}
  if (tissue_type == 'Club') {color <- 'pink'}
  color
}
choose_color <- Vectorize(choose_color)

own_network_visualize <- function(network) {
  ncells <- quantile_data %>%
    dplyr::select(cell_type_epi, cell_id) %>%
    unique() %>%
    group_by(cell_type_epi) %>%
    dplyr::summarize('ncells' = n())
  
  
  data <- quantile_data %>%
    dplyr::filter(source_gene %in% network, target_gene %in% network) %>%
    dplyr::group_by(cell_type_epi, source_gene, target_gene) %>%
    dplyr::summarize('Interactions' = n()) %>%
    ungroup() %>%
    left_join(ncells) %>%
    dplyr::mutate('IpC' = Interactions/ncells)
  data
}

network_data <- own_network_visualize(N3)

graph_all_patients <- network_data %>%
  dplyr::select(from =source_gene, to = target_gene,  tissue_type = cell_type_epi, weight = IpC) %>%
  graph_from_data_frame(directed=F)

E(graph_all_patients)$color <- choose_color(E(graph_all_patients)$tissue_type)
E(graph_all_patients)$color
E(graph_all_patients)$width <- E(graph_all_patients)$weight*50

plot_network_data <- function(network, sometitle) {
  network_data <- own_network_visualize(network)
  
  graph_all_patients <- network_data %>%
    dplyr::select(from =source_gene, to = target_gene,  tissue_type = cell_type_epi, weight = IpC) %>%
    graph_from_data_frame(directed=F)
  
  E(graph_all_patients)$color <- choose_color(E(graph_all_patients)$tissue_type)
  E(graph_all_patients)$color
  E(graph_all_patients)$width <- E(graph_all_patients)$weight*30
  set.seed(0)
  l <- layout_with_fr(graph_all_patients, niter=500)
  par(bg='black')
  plot(graph_all_patients, vertex.size=0.1,
       edge.curved = curve_multiple(graph_all_patients, start = 0.04),
       edge.width = E(graph_all_patients)$width,
       vertex.label.color = 'white',
       vertex.color=NA,
       vertex.label.cex= 4.0,
       main= list(sometitle, col = 'white', cex = 4),
       
       layout=l)
}

png('./figures/graph_subnetwork.png', width=500, height=1200)
#par(bg='black', mfrow = c(1,3), mar = c(9,9,9,9), oma = c(2,2,2,2))
layout(matrix(c(1,2,3), nc=1, byrow=T))
plot_network_data(N1, 'T1')
plot_network_data(N2, 'T2')
plot_network_data(N3, 'T3')
dev.off()

####################################
#quantify Communities (C1-9) for results
###################################
cluster<-4
community_data <- quantile_data %>% filter(source_gene %in% degreeframe[degreeframe$cluster==cluster,]$name, target_gene %in% degreeframe[degreeframe$cluster==cluster,]$name)
genes <- c(community_data$source_gene, community_data$target_gene) %>% unique()
write.csv(genes,quote=F,  './figures/genes.csv', row.names=F, col.names=F)
length(genes)

all_cells <- quantile_data$cell_id %>% unique %>% length
cells_by_celltype <- quantile_data %>% dplyr::select(cell_id, cell_type_epi) %>%
  unique() %>%
  group_by(cell_type_epi) %>%
  dplyr::summarize(cellcount = n())

cells_by_celltypeandpatient <- quantile_data %>% dplyr::select(cell_id, cell_type_epi, patient_id) %>%
  unique() %>%
  group_by(cell_type_epi, patient_id) %>%
  dplyr::summarize(cellcount = n())

interactions_overall <- community_data %>% 
  dplyr::select(source_gene,target_gene, cell_id) %>%
  group_by(source_gene, target_gene) %>%
  dplyr::summarize('IpC' = n()/all_cells)

interactions_by_celltype <- community_data %>% 
  dplyr::select(source_gene,target_gene, cell_id, cell_type_epi) %>%
  group_by(source_gene, target_gene, cell_type_epi) %>%
  dplyr::summarize('interaction_number' = n()) %>%
  left_join(cells_by_celltype) %>%
  dplyr::mutate(IpC = interaction_number/cellcount)

network_by_celltype <-  community_data %>% 
  dplyr::select(source_gene,target_gene, cell_id, cell_type_epi) %>%
  group_by(cell_type_epi) %>%
  dplyr::summarize('interaction_number' = n()) %>%
  left_join(cells_by_celltype) %>%
  dplyr::mutate(IpC = interaction_number/cellcount)
  

interactions_by_celltype_and_patient <- community_data %>% 
  dplyr::select(source_gene,target_gene, cell_id, cell_type_epi, patient_id) %>%
  group_by(source_gene, target_gene, cell_type_epi, patient_id) %>%
  dplyr::summarize('interaction_number' = n()) %>%
  left_join(cells_by_celltypeandpatient) %>%
  dplyr::mutate(IpC = interaction_number/cellcount) %>%
  filter(cellcount >10)

network_by_celltype_and_patient <-  community_data %>% 
  dplyr::select(source_gene,target_gene, cell_id, cell_type_epi, patient_id) %>%
  group_by(cell_type_epi, patient_id) %>%
  dplyr::summarize('interaction_number' = n()) %>%
  left_join(cells_by_celltypeandpatient) %>%
  dplyr::mutate(IpC = interaction_number/cellcount)

######
#compare network with TP53 mutation
data <- network_by_celltype_and_patient %>% dplyr::filter(cell_type_epi == 'Tumor') %>%
  left_join(patient_mutations)

wilcox.test(data$TP53, data$IpC)

wilcox.test(data$EML4, data$IpC)
ggplot(data, aes(x=as.factor(EML4), y= IpC)) + geom_point()
#####
####################################
#quantify T1-3 for results
###################################
cluster<-N3
community_data <- quantile_data %>% filter(source_gene %in% cluster, target_gene %in% cluster)
genes <- c(community_data$source_gene, community_data$target_gene) %>% unique()
length(genes)

#interactions_overall <- community_data %>% 
#  dplyr::select(source_gene,target_gene, cell_id) %>%
#  group_by(source_gene, target_gene) %>%
#  dplyr::summarize('IpC' = n()/all_cells)

#interactions_by_celltype <- community_data %>% 
#  dplyr::select(source_gene,target_gene, cell_id, cell_type_epi) %>%
#  group_by(source_gene, target_gene, cell_type_epi) %>%
#  dplyr::summarize('interaction_number' = n()) %>%
#  left_join(cells_by_celltype) %>%
#  dplyr::mutate(IpC = interaction_number/cellcount)

network_by_celltype <-  community_data %>% 
  dplyr::select(source_gene,target_gene, cell_id, cell_type_epi) %>%
  group_by(cell_type_epi) %>%
  dplyr::summarize('interaction_number' = n()) %>%
  left_join(cells_by_celltype) %>%
  dplyr::mutate(IpC = interaction_number/cellcount)


interactions_by_celltype_and_patient <- community_data %>% 
  dplyr::select(source_gene,target_gene, cell_id, cell_type_epi, patient_id) %>%
  group_by(source_gene, target_gene, cell_type_epi, patient_id) %>%
  dplyr::summarize('interaction_number' = n()) %>%
  left_join(cells_by_celltypeandpatient) %>%
  dplyr::mutate(IpC = interaction_number/cellcount) %>%
  filter(cellcount >10)

network_by_celltype_and_patient <-  community_data %>% 
  dplyr::select(source_gene,target_gene, cell_id, cell_type_epi, patient_id) %>%
  group_by(cell_type_epi, patient_id) %>%
  dplyr::summarize('interaction_number' = n()) %>%
  left_join(cells_by_celltypeandpatient) %>%
  dplyr::mutate(IpC = interaction_number/cellcount)

#see interactions for T2 when p032 is excluded
network_ungrouped <- dim(community_data)[1] / length(quantile_data %>% dplyr::filter(patient_id != 'p032') %>% .$cell_id %>% unique())

network_ungrouped
network_by_celltype
network_by_celltype_and_patient %>% filter(cell_type_epi == 'Tumor')

###########################
#### find objective tumor specific networks
###########################
quantile_data <-  quantile_data %>%
  dplyr::mutate('istumor' = ifelse(cell_type_epi =='Tumor', 'Tumor', 'NoTumor'))

cellcounts <- quantile_data %>% 
  dplyr::select(cell_id, istumor) %>%
  unique() %>%
  group_by(istumor) %>%
  dplyr::summarize(N = n())
  
tumor_specific <- quantile_data %>%
  dplyr::group_by(source_gene, target_gene, istumor) %>%
  dplyr::summarize(interactioncount = n()) %>%
  left_join(cellcounts) %>%
  group_by(source_gene, target_gene) %>%
  ungroup() %>%
  dplyr::mutate(interactioncount = interactioncount/N) %>%
  dplyr::group_by(source_gene, target_gene) %>%
  dplyr::mutate(maxIpC = max(interactioncount)) %>%
  ungroup() %>%
  dplyr::filter(maxIpC > 0.03) %>% #0.04
  dplyr::select(istumor, interactioncount, source_gene, target_gene) %>%
  pivot_wider(names_from=istumor, values_from=interactioncount) %>%
  dplyr::mutate(Tumor = ifelse(is.na(Tumor), 0, Tumor)) %>%
  dplyr::mutate(NoTumor = ifelse(is.na(NoTumor), 0, NoTumor)) %>%
  dplyr::mutate(ratio = ifelse(Tumor>1000*NoTumor,log(1000), log(Tumor/NoTumor))) %>%
  dplyr::filter(!is.na(ratio), ratio>3) # 4
  
write.csv(tumor_specific, './figures/tumor_edges.csv')

tumor_specific_frame <- quantile_data %>%
  dplyr::filter(paste0(source_gene,target_gene) %in% paste0(tumor_specific$source_gene, tumor_specific$target_gene))

tumor_specific_graph <- tumor_specific_frame  %>%
  dplyr::select(source_gene, target_gene) %>%
  as.matrix() %>%
    graph_from_edgelist() %>%
    as.undirected()

tumor_specific_graph <- quantile_data %>%
  dplyr::filter(paste0(source_gene,target_gene) %in% paste0(tumor_specific$source_gene, tumor_specific$target_gene)) %>%
  dplyr::select(source_gene, target_gene) %>%
  group_by(source_gene, target_gene) %>%
  dplyr::summarize(weight = 0.001 * n()) %>%
  graph_from_data_frame(directed=F)

E(tumor_specific_graph)$weight
E(tumor_specific_graph)$width <- E(tumor_specific_graph)$weight * 5
E(tumor_specific_graph)$color = 'black'

png('./figures/tumor_specific_graph', width=2000, height=2000)
par(bg='black')
plot(tumor_specific_graph, vertex.size=0.1, vertex.label.cex= 2, edge.color='red', vertex.label.color= 'white')
     #edge.width = E(graph_all_patients)$width,)
dev.off()
###################################
#####same for normal epithel
###################################

epit_specific <- quantile_data %>%
  dplyr::group_by(source_gene, target_gene, istumor) %>%
  dplyr::summarize(interactioncount = n()) %>%
  left_join(cellcounts) %>%
  group_by(source_gene, target_gene) %>%
  ungroup() %>%
  dplyr::mutate(interactioncount = interactioncount/N) %>%
  dplyr::group_by(source_gene, target_gene) %>%
  dplyr::mutate(maxIpC = max(interactioncount)) %>%
  dplyr::filter(maxIpC > 0.04) %>%
  ungroup()%>%
  dplyr::select(istumor, interactioncount, source_gene, target_gene) %>%
  pivot_wider(names_from=istumor, values_from=interactioncount) %>%
  dplyr::mutate(Tumor = ifelse(is.na(Tumor), 0, Tumor)) %>%
  dplyr::mutate(NoTumor = ifelse(is.na(NoTumor), 0, NoTumor)) %>%
  dplyr::mutate(ratio = ifelse(NoTumor>1000*Tumor,log(1000), log(NoTumor/Tumor))) %>%
  dplyr::filter(!is.na(ratio), ratio>4)


epit_specific_graph <- quantile_data %>%
  dplyr::filter(paste0(source_gene,target_gene) %in% paste0(epit_specific$source_gene, epit_specific$target_gene)) %>%
  dplyr::select(source_gene, target_gene) %>%
  group_by(source_gene, target_gene) %>%
  dplyr::summarize(weight = 0.001 * n()) %>%
  graph_from_data_frame(directed=F) #%>%
as.undirected()

#tumor_specific_graph <- simplify(tumor_specific_graph, edge.attr.comb =list(weight='sum'))
E(epit_specific_graph)$weight
E(epit_specific_graph)$width <- E(epit_specific_graph)$weight * 10
E(epit_specific_graph)$color = 'black'

png('./figures/epit_specific_graph', width=2000, height=2000)
par(bg='black')
plot(epit_specific_graph, vertex.size=0.1, vertex.label.cex= 2, edge.color='steelblue', vertex.label.color= 'white')
dev.off()
#########
# tumor specific modules and network visualization
#########
#######
set.seed(0)
tumor_specific_modules <- cluster_louvain(tumor_specific_graph)
print(tumor_specific_modules)
tumor_specific_modules$membership

degree(tumor_specific_graph, V(tumor_specific_graph))
tumor_specific_degreeframe <- data.frame(name = V(tumor_specific_graph)$name, degree = degree(tumor_specific_graph, V(tumor_specific_graph)), cluster = tumor_specific_modules$membership) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(N = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(N>2)

tumor_specific_degreeframe$clusterrank <- as.numeric(as.factor(tumor_specific_degreeframe$cluster))

cluster_sizes <-tumor_specific_degreeframe %>% group_by(cluster) %>% dplyr::summarize(N = n()) %>% filter(N>1)

plot_tumor_clusters <- function(module_id) {
  df <- tumor_specific_degreeframe %>% dplyr::filter(clusterrank == module_id) 
  tumor_specific_frame <- quantile_data %>%
    dplyr::filter(paste0(source_gene,target_gene) %in% paste0(tumor_specific$source_gene, tumor_specific$target_gene))
  spec_data <- tumor_specific_frame %>% dplyr::filter(source_gene %in% df$name, target_gene %in% df$name)
    
  cell_count <- quantile_data %>%
      dplyr::select(cell_id, cell_type_epi) %>%
      unique() %>%
      dplyr::group_by(cell_type_epi) %>%
      dplyr::summarize(N = n())
    
  tumor_specific_subedges <-  spec_data %>%
      dplyr::group_by(source_gene, target_gene, cell_type_epi) %>%
      dplyr::summarize(interactioncount = n()) %>%
      left_join(cell_count) %>%
      group_by(source_gene, target_gene) %>%
      ungroup() %>%
      dplyr::mutate(interactioncount = interactioncount/N) %>%
      dplyr::select(from =source_gene, to = target_gene,  tissue_type = cell_type_epi, weight = interactioncount)
      
  tumor_specific_subgraph <- graph_from_data_frame(tumor_specific_subedges, directed=F)
  #tumor_specific_graph <- graph_from_edgelist(tumor_specific_edges %>% dplyr::select(from, to), directed=F)
  
    E(tumor_specific_subgraph)$color <- choose_color(E(tumor_specific_subgraph)$tissue_type)
    E(tumor_specific_subgraph)$color
    E(tumor_specific_subgraph)$width <- E(tumor_specific_subgraph)$weight*100
    set.seed(0)
    l <- layout_with_fr(tumor_specific_subgraph, niter=500)
  
    plot(tumor_specific_subgraph, vertex.size=0.2,
         edge.curved = curve_multiple(tumor_specific_subgraph, start = 0.05),
         edge.width = E(tumor_specific_subgraph)$width,
         vertex.label.color = 'black',
         vertex.color=NA,
         vertex.label.cex= 3.5,
         layout=l)
    
    text(-1.2,1.0, labels = paste0('T', module_id), cex=5)
    
}
clusters <- tumor_specific_degreeframe$clusterrank %>% unique() %>% sort()
#par(mfrow = c(1,length(clusters)))
#dev.off()
png('./figures/tumor_specific_graph_3', width=600, height = 500*6 )
par(bg='white', mfrow = c(6,1))
for (i in clusters) {
plot_tumor_clusters(i)
}
dev.off()

#########c

# tumor specific modules as heatmap
#########

tumor_specific_degreeframe$clusterrank <- as.numeric(as.factor(tumor_specific_degreeframe$cluster))

get_relative_counts <- function(idx, cell_type) {
  
  gene_set <- tumor_specific_degreeframe %>% filter(degree>-1, cluster==idx) %>% .$name
  clusterrank <- tumor_specific_degreeframe %>% filter(degree>-1, cluster==idx)  %>% .$clusterrank %>% unique()
  tumor_specific_frame <- quantile_data %>%
    dplyr::filter(paste0(source_gene,target_gene) %in% paste0(tumor_specific$source_gene, tumor_specific$target_gene))
  spec_data <- tumor_specific_frame %>% dplyr::filter(source_gene %in% gene_set, target_gene %in% gene_set)
  
  ncells_perid <- quantile_data %>% 
    filter(cell_type_epi == cell_type) %>%
    dplyr::select(patient_id, sample_name) %>%
    unique() %>%
    group_by(patient_id) %>%
    dplyr::summarize('all_cellcount' = n())
  
  count_data <- spec_data %>% 
    filter(cell_type_epi == cell_type) %>%
    dplyr::select(patient_id, sample_name) %>%
    group_by(patient_id) %>%
    dplyr::summarize('spec_cell_count' = n())
  
  intermediate <- left_join(ncells_perid, count_data) %>% 
    dplyr::mutate('ratio' = spec_cell_count/all_cellcount)

  data.frame('patient_id' = intermediate$patient_id, 'clustervalue' = intermediate$ratio, 'cluster' = clusterrank)
}


get_all_clusters <- function(cell_type) {
  threegeneclusters <- tumor_specific_degreeframe %>% dplyr::filter(N>2) %>% .$cluster %>%unique() %>% sort()
  all_clusters <- rbindlist(lapply(threegeneclusters, function(x) get_relative_counts(x, cell_type)))
  all_clusters$cell_type_epi <- cell_type
  all_clusters
}

get_all_clusters('Tumor')
get_relative_counts(6, 'Tumor')

cell_types <- quantile_data %>% .$cell_type_epi %>% unique()

all_data <- rbindlist(lapply(cell_types, get_all_clusters))
normalized_data <- all_data %>% dplyr::group_by(cluster) %>%
  dplyr::mutate(maxv = max(clustervalue, na.rm=T), minv = min(clustervalue, na.rm=T)) %>%
  dplyr::mutate(normv = (clustervalue)/(maxv))

ht_opt(legend_border = "black",
       heatmap_border = TRUE)
col_fun = colorRamp2(c(0, 1), c("white", "red"))
col_fun2 = colorRamp2(c(0, 1), c("white", "black"))
w <- unit(4, "cm")
h <- unit(6, "cm")

mutant_heatmap <- p_mat %>% 
  Heatmap(cluster_columns = F, width=0.5*w, height=h, name = 'Mutation', col = col_fun2, 
          show_heatmap_legend = F, column_names_gp = grid::gpar(fontsize = 15), column_title = 'Mutation')

h_mat <- patient_histology %>% dplyr::select(-patient_id) %>% as.matrix()
rownames(h_mat) <- patient_histology$patient_id
histology_heatmap <- h_mat %>% Heatmap(cluster_columns = F, width=0.3*w, height=h, name = 'Histology', col=(rainbow(14)[c(1,3,4,6,8,10)]), row_names_gp = gpar(fontsize=15),
                                       column_names_gp = gpar(fontsize=15), heatmap_legend_param =  list(title_gp =  gpar(fontsize=13), labels_gp = gpar(fontsize = 13)))


prepare_matrix <- function(cell_type) {
  all_clusters_per_celltype <- normalized_data %>% dplyr::filter(cell_type_epi == cell_type) %>%
    dplyr::select(patient_id, cell_type_epi, normv, -cell_type_epi) %>%
    dplyr::mutate(cluster = paste0('T', cluster)) %>%
    pivot_wider(names_from = cluster, values_from = normv)
  
  all_clusters_per_celltype[is.na(all_clusters_per_celltype)] <- 0
  
  all_clusters_per_celltype <- all_clusters_per_celltype %>%
    merge(patient_ids, all=T)
  cell_type_matrix <- all_clusters_per_celltype %>% 
    dplyr::select(-patient_id) %>%
    as.matrix()
  rownames(cell_type_matrix) <- all_clusters_per_celltype$patient_id
  cell_type_matrix
  
}
prepare_matrix('AT1')


png('./figures/patient_heatmap_tumor_specific.png', width=1000, height=250)
Heatmap(prepare_matrix('Tumor'), col = col_fun, cluster_columns = F, column_title ='Tumor', name = 'IpC (scaled)', column_names_gp = gpar(fontsize=15),column_names_rot = 90,
        heatmap_legend_param =  list(title_gp =  gpar(fontsize=13), labels_gp = gpar(fontsize = 13))) + 
  Heatmap(prepare_matrix('Club'), col = col_fun, cluster_columns = F, column_title ='Club',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15),column_names_rot = 90)+  
  Heatmap(prepare_matrix('AT1'), col = col_fun, cluster_columns = F, column_title ='AT1',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15),column_names_rot = 90) + 
  Heatmap(prepare_matrix('AT2'), col = col_fun, cluster_columns = F, column_title ='AT2',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15),column_names_rot = 90) + 
  Heatmap(prepare_matrix('Ciliated'), col = col_fun, cluster_columns = F, column_title ='Ciliated',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15),column_names_rot = 90) + 
  mutant_heatmap+
  histology_heatmap
dev.off()

#####################
###exmine tumor interactions
######################
Tx <- get_all_clusters('Tumor')
idx = 6
gene_set <- tumor_specific_degreeframe %>% filter(degree>-1, clusterrank==idx) %>% .$name
#clusterrank <- tumor_specific_degreeframe %>% filter(degree>-1, cluster==idx)  %>% .$clusterrank %>% unique()
write.csv(gene_set, './figures/tumor_genes.csv',row.names=F, col.names=F,quote=F)

tumor_specific_frame <- quantile_data %>%
  dplyr::filter(paste0(source_gene,target_gene) %in% paste0(tumor_specific$source_gene, tumor_specific$target_gene))
spec_data <- tumor_specific_frame %>% dplyr::filter(source_gene %in% gene_set, target_gene %in% gene_set)

ncells_tumor <- quantile_data %>% filter(cell_type_epi == 'Tumor') %>%
  dplyr::select(patient_id, sample_name) %>% unique() %>% dim() %>% .[1]

#here IpC over all patients
interaction_IpC_average <- spec_data %>% 
  filter(cell_type_epi == 'Tumor') %>%
  dplyr::select(patient_id, sample_name, source_gene, target_gene) %>%
  group_by(source_gene, target_gene) %>%
  dplyr::summarize('spec_cell_count' = n()/ncells_tumor)

# now IpC stratified by patient
ncells_tumor_per_patient <- quantile_data %>% filter(cell_type_epi == 'Tumor') %>%
  dplyr::select(patient_id, sample_name) %>% unique() %>% group_by(patient_id) %>% dplyr::summarize(N = n())


interaction_IpC_per_patient <-spec_data %>% 
  filter(cell_type_epi == 'Tumor') %>%
  dplyr::select(patient_id, sample_name, source_gene, target_gene) %>%
  group_by(source_gene, target_gene,patient_id) %>%
  dplyr::summarize('spec_cell_count' = n()) %>%
  left_join(ncells_tumor_per_patient) %>%
  dplyr::mutate(IpC = spec_cell_count/N)

