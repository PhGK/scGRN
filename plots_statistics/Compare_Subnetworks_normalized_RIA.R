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
patient_ids <- quantile_data %>% dplyr::select(patient_id) %>% unique() 
cell_types <- quantile_data %>% .$cell_type_epi %>% unique()


####
ntissue <- quantile_data %>% dplyr::select(sample_name, cell_type_epi)%>%
  unique() %>%
  .$cell_type_epi %>%as.factor %>% summary()
ntissue
####
######################

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

all_graph_edges <- quantile_data %>%
  group_by(cell_type_epi, patient_id) %>%
  dplyr::mutate(IpC = 1/(n()+10)) %>%
  dplyr::group_by(source_gene, target_gene) %>%
  dplyr::summarize(IpC = sum(IpC)) %>%
  dplyr::filter(IpC>1e-2) %>%
  dplyr::select(-IpC)


E(all_graph)$weight <-1
#########
# number of edges berücksichtigen
#########
#######
set.seed(0)
modules <- cluster_louvain(all_graph)
print(modules)
modules$membership

degreeframe <- data.frame(name = V(all_graph)$name, degree = degree(all_graph, V(all_graph)), cluster = modules$membership) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(N = n()) %>%
  dplyr::ungroup() #%>%
  #dplyr::filter(N!=max(N))

cluster_sizes <-degreeframe %>% group_by(cluster) %>% dplyr::summarize(N = n()) %>% filter(N>1)

first_gene <- degreeframe %>% dplyr::select('first_gene' = name, 'cluster1' = cluster)
second_gene <- degreeframe %>% dplyr::select('second_gene' = name, 'cluster2'  =cluster)

edge_number_frame <- all_graph_edges %>%
  left_join(first_gene, by = c('source_gene' = 'first_gene')) %>%
  left_join(second_gene, by = c('target_gene' = 'second_gene')) %>%
  dplyr::filter(cluster1==cluster2) %>%
  dplyr::group_by(cluster1) %>%
  dplyr::summarize('N_edges' = n()) %>%
  dplyr::mutate(cluster=cluster1) %>%
  select(-cluster1)

###

################################################
##### long format to normalize per cluster
################################################

###
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
Heatmap(prepare_matrix('Tumor'), col = col_fun, cluster_columns = F, column_title ='Tumor', name = 'RIA (scaled)', column_names_gp = gpar(fontsize=15),
        heatmap_legend_param =  list(title_gp =  gpar(fontsize=13), labels_gp = gpar(fontsize = 13))) + 
  Heatmap(prepare_matrix('Club'), col = col_fun, cluster_columns = F, column_title ='Club',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15))+  
  Heatmap(prepare_matrix('AT1'), col = col_fun, cluster_columns = F, column_title ='AT1',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15)) + 
  Heatmap(prepare_matrix('AT2'), col = col_fun, cluster_columns = F, column_title ='AT2',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15)) + 
  Heatmap(prepare_matrix('Ciliated'), col = col_fun, cluster_columns = F, column_title ='Ciliated',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15)) + 
  mutant_heatmap+
  histology_heatmap
dev.off()

################################################
##### now normalize by edge number of blueprint
################################################
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

all_data <- rbindlist(lapply(cell_types, get_all_clusters))
normalized_data <- all_data %>% 
  left_join(edge_number_frame) %>%
  dplyr::mutate(normv = clustervalue/N_edges, 'cluster_' = cluster)

patient_mutations <- data.frame(patient_id = c('p018', 'p019','p023','p024','p027','p030','p031','p032','p033','p034'),
                                'KRAS' = c(1,0,1,0,0,1,1,1,1,0),
                                'TP53' = c(0,0,1,0,1,0,0,0,0,0),
                                'PIK3CA' = c(0,0,0,0,0,0,1,0,0,0),
                                'EML4' = c(0,0,0,1,0,0,0,0,0,0))

patient_histology <- data.frame(patient_id = c('p018', 'p019','p023','p024','p027','p030','p031','p032','p033','p034'),
                                Histology = as.factor(c('acinar', 'acinar', 'solid', 'acinar', 'solid', 'mucinuous', 'acinar', 'lepidic', 'papillary', 'sarcomatoid')))

p_mat <- patient_mutations %>% dplyr::select(-patient_id) %>% as.matrix()
rownames(p_mat) <- patient_mutations$patient_id

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
    dplyr::select(patient_id, cell_type_epi, normv, cluster,  -cell_type_epi) %>%
    dplyr::mutate(cluster = paste0('C', as.character(cluster))) %>%
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
prepare_matrix('AT2')


ht_opt(legend_border = "black",
       heatmap_border = TRUE)
col_fun = colorRamp2(c(0, 1), c("white","red"))
col_fun2 = colorRamp2(c(0, 1), c("white", "black"))


png('./figures/patient_heatmap_byblueprint.png', width=1000, height=250)
Heatmap(prepare_matrix('Tumor'), col = col_fun, cluster_columns = F, column_title ='Tumor', name = 'RIA', column_names_gp = gpar(fontsize=15),
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
layout.matrix <- matrix(c(1,1,2,2,3,3,1,1,2,2,3,3,4,5,6,7,8,9), nrow= 6)


png('./figures/general_specific_graph', width=600*3, height = 500*6)
layout(mat = layout.matrix)
par(mar=c(0,0,0,0))
for (i in c(1,3,6,2,4,5,7,8,9)) {
  plot_general_clusters(i)
}
dev.off()


####################################
#quantify Communities (C1-9) for results
###################################
cluster<-9
sel_cluster <- cluster
community_data <- quantile_data %>% 
  filter(source_gene %in% degreeframe[degreeframe$cluster==cluster,]$name, target_gene %in% degreeframe[degreeframe$cluster==cluster,]$name) %>%
  dplyr::mutate('n_edges' = edge_number_frame %>% filter(cluster==sel_cluster) %>% .$N_edges)
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
  dplyr::select(source_gene,target_gene, cell_id, cell_type_epi, n_edges) %>%
  group_by(cell_type_epi) %>%
  dplyr::summarize('interaction_number' = n(), n_edges = first(n_edges)) %>%
  left_join(cells_by_celltype) %>%
  dplyr::mutate(IpC = interaction_number/cellcount) %>%
  dplyr::mutate(normalized_IpC = IpC/n_edges)
  
interactions_by_celltype_and_patient <- community_data %>% 
  dplyr::select(source_gene,target_gene, cell_id, cell_type_epi, patient_id) %>%
  group_by(source_gene, target_gene, cell_type_epi, patient_id) %>%
  dplyr::summarize('interaction_number' = n()) %>%
  left_join(cells_by_celltypeandpatient) %>%
  dplyr::mutate(IpC = interaction_number/cellcount) %>%
  filter(cellcount >10)

network_by_celltype_and_patient <- community_data %>% 
  dplyr::select(source_gene,target_gene, cell_id, cell_type_epi, patient_id, n_edges) %>%
  group_by(cell_type_epi, patient_id) %>%
  dplyr::summarize('interaction_number' = n(), n_edges = first(n_edges)) %>%
  left_join(cells_by_celltypeandpatient) %>%
  dplyr::mutate(IpC = interaction_number/cellcount)  %>%
  dplyr::mutate(normalized_IpC = IpC/n_edges)

######
#compare network with TP53 mutation
data <- network_by_celltype_and_patient %>% dplyr::filter(cell_type_epi == 'Tumor') %>%
  left_join(patient_mutations)

wilcox.test(data$TP53, data$IpC)

wilcox.test(data$EML4, data$IpC)
ggplot(data, aes(x=as.factor(EML4), y= IpC)) + geom_point()
#####

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
  dplyr::filter(maxIpC > 0.04) %>% #0.04
  dplyr::select(istumor, interactioncount, source_gene, target_gene) %>%
  pivot_wider(names_from=istumor, values_from=interactioncount) %>%
  dplyr::mutate(Tumor = ifelse(is.na(Tumor), 0, Tumor)) %>%
  dplyr::mutate(NoTumor = ifelse(is.na(NoTumor), 0, NoTumor)) %>%
  dplyr::mutate(ratio = ifelse(Tumor>1000*NoTumor,log(1000), log(Tumor/NoTumor))) %>%
  dplyr::filter(!is.na(ratio), ratio>4) # 4
  
write.csv(tumor_specific, './figures/tumor_edges.csv')

tumor_specific_frame <- quantile_data %>%
  dplyr::filter(paste0(source_gene,target_gene) %in% paste0(tumor_specific$source_gene, tumor_specific$target_gene))

#tumor_specific_graph <- tumor_specific_frame  %>%
#  dplyr::select(source_gene, target_gene) %>%
#  as.matrix() %>%
#    graph_from_edgelist() %>%
#    as.undirected()

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
  graph_from_data_frame(directed=F) %>%
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

first_tumorgene <- tumor_specific_degreeframe %>% dplyr::select('first_gene' = name, 'clusterrank1' = clusterrank)
second_tumorgene <- tumor_specific_degreeframe %>% dplyr::select('second_gene' = name, 'clusterrank2'  =clusterrank)

#count edges for RIA normalization
tumor_specific_edge_number_frame <- tumor_specific_frame %>%
  dplyr::select(source_gene, target_gene) %>%
  unique() %>%
  left_join(first_tumorgene, by = c('source_gene' = 'first_gene')) %>%
  left_join(second_tumorgene, by = c('target_gene' = 'second_gene')) %>%
  dplyr::filter(clusterrank1==clusterrank2) %>%
  dplyr::group_by(clusterrank1) %>%
  dplyr::summarize('N_edges' = n()) %>%
  dplyr::mutate(clusterrank=clusterrank1) %>%
  select(-clusterrank1)


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

#########
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
Heatmap(prepare_matrix('Tumor'), col = col_fun, cluster_columns = F, column_title ='Tumor', name = 'RIA (scaled)', column_names_gp = gpar(fontsize=15),column_names_rot = 90,
        heatmap_legend_param =  list(title_gp =  gpar(fontsize=13), labels_gp = gpar(fontsize = 13))) + 
  Heatmap(prepare_matrix('Club'), col = col_fun, cluster_columns = F, column_title ='Club',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15),column_names_rot = 90)+  
  Heatmap(prepare_matrix('AT1'), col = col_fun, cluster_columns = F, column_title ='AT1',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15),column_names_rot = 90) + 
  Heatmap(prepare_matrix('AT2'), col = col_fun, cluster_columns = F, column_title ='AT2',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15),column_names_rot = 90) + 
  Heatmap(prepare_matrix('Ciliated'), col = col_fun, cluster_columns = F, column_title ='Ciliated',show_heatmap_legend = F, column_names_gp = gpar(fontsize=15),column_names_rot = 90) + 
  mutant_heatmap+
  histology_heatmap
dev.off()

#####################
###examine tumor interactions
######################
Tx <- get_all_clusters('Tumor')
idx = 4
gene_set <- tumor_specific_degreeframe %>% filter(degree>-1, clusterrank==idx) %>% .$name
#clusterrank <- tumor_specific_degreeframe %>% filter(degree>-1, cluster==idx)  %>% .$clusterrank %>% unique()
write.csv(gene_set, './figures/tumor_genes.csv',row.names=F, col.names=F,quote=F)

tumor_specific_frame <- quantile_data %>%
  dplyr::filter(paste0(source_gene,target_gene) %in% paste0(tumor_specific$source_gene, tumor_specific$target_gene))
spec_data <- tumor_specific_frame %>% dplyr::filter(source_gene %in% gene_set, target_gene %in% gene_set) %>%
  dplyr::mutate(clusterrank = idx)

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

network_IpC_per_patient <- spec_data %>% 
  filter(cell_type_epi == 'Tumor') %>%
  dplyr::select(patient_id, sample_name, source_gene, target_gene, clusterrank) %>%
  group_by(patient_id) %>%
  dplyr::summarize('spec_cell_count' = n(), clusterrank = first(clusterrank)) %>%
  left_join(ncells_tumor_per_patient) %>%
  dplyr::mutate(IpC = spec_cell_count/N) %>%
  left_join(tumor_specific_edge_number_frame) %>%
  dplyr::mutate(normalized_IpC = IpC / N_edges)
  
