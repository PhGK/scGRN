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

quant_edges <- 0.8
test_data_tumor <- fread(paste('../data/high_values_concat.csv', sep = ""))

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
  unique() #%>%
  dplyr::group_by(patient_id) %>%
  dplyr::summarize('counts' = n())

#return_sample('Neuroendocrine', cells, ncells_per_tissue)
cell_counts <-quantile_filtered_test_data_all_patients %>% dplyr::select(sample_name,cell_type_epi) %>%
  unique() %>%
  group_by(cell_type_epi) %>%
  dplyr::summarize('counts' = n())

######################
######################
######################
######################
#find modules in qftdap
all_graph <- quantile_filtered_test_data_all_patients %>%
  dplyr::select(source_gene, target_gene) %>% 
  as.matrix() %>%
  graph_from_edgelist() %>%
  as.undirected()

E(all_graph)$weight <-1
all_graph <- simplify(all_graph)#, edge.attr.comb =list(weight='sum'))
#all_graph <- delete_edges(all_graph,E(all_graph)[E(all_graph)$weight<2])
#E(all_graph)$weight <-1

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
  dplyr::ungroup() %>%
  dplyr::filter(N!=max(N))
cluster_sizes <-degreeframe %>% group_by(cluster) %>% dplyr::summarize(N = n()) %>% filter(N>1)
###

##########
#first find prior distribution

gene_expression <- fread('../data/epi_top2000.csv')
get_all_interactions <- function(gene1, gene2, cluster) {
  expr_data <- gene_expression %>% dplyr::select(gene1, gene2, cell_type_epi, patient_id) %>%
    dplyr::mutate(inter_exists = (.[,1]!=0)&(.[,2]!=0)) %>%
    dplyr::group_by(cell_type_epi, patient_id) %>%
    dplyr::summarize('n_interactions' = sum(inter_exists))
  expr_data$gene1 <- gene1
  expr_data$gene2 <- gene2
  expr_data$cluster <- cluster
  expr_data
  #print(expr_data)
}

get_priors_per_cluster <- function(cluster_id){ 
  current_frame <- degreeframe %>% filter(N>1, cluster==cluster_id) 
  if (dim(current_frame[1])==0) return(NULL) 
  
  all_genes <- list()
  i = 1
  for (gene1 in current_frame$name) {
    for (gene2 in current_frame$name) {
      #get_all_interactions(gene1, gene2)
      all_genes[[i]] <- get_all_interactions(gene1, gene2, cluster = cluster_id)
      i <- i + 1
      if (i>1000) break
      
    }
  }
  all_genes_frame <- rbindlist(all_genes) %>% dplyr::filter(gene1>gene2) %>%
    group_by(cell_type_epi, patient_id) %>%
    dplyr::summarize('n_interactions' = sum(n_interactions), 'cluster' = mean(cluster))
all_genes_frame
}


all_priors <- rbindlist(lapply(seq(10), get_priors_per_cluster)) %>% filter(cell_type_epi != 'Neuroendocrine')

##########
##########
###
patient_ids <- quantile_filtered_test_data_all_patients %>% dplyr::select(patient_id) %>% unique() 
get_relative_counts <- function(idx, cell_type) {
  
  gene_set <- degreeframe %>% filter(degree>10, cluster==idx) %>% .$name
  
  ncells_perid <- quantile_filtered_test_data_all_patients %>% 
    filter(cell_type_epi == cell_type) %>%
    dplyr::select(patient_id, sample_name) %>%
    unique() %>%
    group_by(patient_id) %>%
    dplyr::summarize('all_cellcount' = n())
  
  
  count_data <- quantile_filtered_test_data_all_patients %>% 
    filter(cell_type_epi == cell_type) %>%
    filter(source_gene %in% gene_set, target_gene %in% gene_set) %>%
    dplyr::select(patient_id, sample_name) %>%
    #unique() %>%
    group_by(patient_id) %>%
    dplyr::summarize('spec_cell_count' = n())
    
    intermediate <- left_join(ncells_perid, count_data) %>% 
      dplyr::mutate('ratio' = spec_cell_count/all_cellcount)
    
    data.frame('patient_id' = intermediate$patient_id, 'clustervalue' = intermediate$ratio)
}


get_all_clusters <- function(cell_type) {
  all_clusters <- get_relative_counts(1, cell_type)
  colnames(all_clusters)[2] <- paste0('C', 1)
  for (i in seq(2,13)) {
    all_clusters <- left_join(all_clusters, get_relative_counts(i, cell_type))
    colnames(all_clusters)[i+1] <- paste0('C', i)
  }
  mask <- apply(all_clusters,2, function(x) all(is.na(x)))
  all_clusters <- all_clusters[,!mask]
  all_clusters[is.na(all_clusters)] <- 0
  full_join(all_clusters, patient_ids)
}

get_all_clusters('Tumor')
get_relative_counts(6, 'Tumor')


patient_mutations <- data.frame(patient_id = c('p018', 'p019','p023','p024','p027','p030','p031','p032','p033','p034'),
                                'KRAS' = c(1,0,1,0,0,1,1,1,1,0),
                                'TP53' = c(0,0,1,0,1,0,0,0,0,0),
                                'PIK3CA' = c(0,0,0,0,0,0,1,0,0,0),
                                'EML4' = c(0,0,0,1,0,0,0,0,0,0))

patient_histology <- data.frame(patient_id = c('p018', 'p019','p023','p024','p027','p030','p031','p032','p033','p034'),
                                histology = as.factor(c('acinar', 'acinar', 'solid', 'acinar', 'solid', 'mucinuous', 'acinar', 'lepidic', 'papillary', 'sarcomatoid')))

p_mat <- patient_mutations %>% dplyr::select(-patient_id) %>% as.matrix()
rownames(p_mat) <- patient_mutations$patient_id

h_mat <- patient_histology %>% dplyr::select(-patient_id) %>% as.matrix()
rownames(h_mat) <- patient_histology$patient_id

prepare_matrix <- function(cell_type) {
  all_clusters_patient <- get_all_clusters(cell_type)
  rownames(all_clusters_patient) <- all_clusters_patient$patient_id
  heatmap_matrix <- all_clusters_patient %>% dplyr::select(-patient_id) %>% as.matrix()
  #heatmap_matrix <- log(1+heatmap_matrix)
  #heatmap_matrix <- apply(heatmap_matrix, 2, function(x) (x-min(x))/max(x-min(x)))
  heatmap_matrix
}
prepare_matrix('AT2') %>% dim()

ht_opt(legend_border = "black",
       heatmap_border = TRUE)
col_fun = colorRamp2(c(0, 30), c("white", "red"))
#Heatmap(prepare_matrix('Tumor'))
w <- unit(4, "cm")
h <- unit(4, "cm")
mutant_heatmap <- p_mat %>% Heatmap(cluster_columns = F, width=0.3*w, height=h, name = 'Mutation')
histology_heatmap <- h_mat %>% Heatmap(cluster_columns = F, width=0.3*w, height=h, name = 'Histology', col=(rainbow(14)[c(1,3,4,6,8,10)]))
histology_heatmap
png('./figures/patient_heatmap.png')
Heatmap(prepare_matrix('Tumor'), col = col_fun, cluster_columns = F, width=w, height=h) + 
  Heatmap(prepare_matrix('Club'), col = col_fun, cluster_columns = F, width=w, height=h)+  
  Heatmap(prepare_matrix('AT1'), col = col_fun, cluster_columns = F, width=w, height=h) + 
  Heatmap(prepare_matrix('AT2'), col = col_fun, cluster_columns = F, width=w, height=h) + 
  Heatmap(prepare_matrix('Ciliated'), col = col_fun, cluster_columns = F, width=w, height=h) + 
  mutant_heatmap
dev.off()

png('./figures/patient_heatmap.png', width=1000, height=200)
Heatmap(prepare_matrix('Tumor'), col = col_fun, cluster_columns = F, column_title ='Tumor') + 
  Heatmap(prepare_matrix('Club'), col = col_fun, cluster_columns = F, column_title ='Club')+  
  Heatmap(prepare_matrix('AT1'), col = col_fun, cluster_columns = F, column_title ='AT1') + 
  Heatmap(prepare_matrix('AT2'), col = col_fun, cluster_columns = F, column_title ='AT2') + 
  Heatmap(prepare_matrix('Ciliated'), col = col_fun, cluster_columns = F, column_title ='Ciliated') + 
  mutant_heatmap
dev.off()


################################################
##### long format to normalize per cluster
################################################

###
patient_ids <- quantile_filtered_test_data_all_patients %>% dplyr::select(patient_id) %>% unique() 
get_relative_counts <- function(idx, cell_type) {
  
  gene_set <- degreeframe %>% filter(degree>10, cluster==idx) %>% .$name
  
  ncells_perid <- quantile_filtered_test_data_all_patients %>% 
    filter(cell_type_epi == cell_type) %>%
    dplyr::select(patient_id, sample_name) %>%
    unique() %>%
    group_by(patient_id) %>%
    dplyr::summarize('all_cellcount' = n())
  
  
  count_data <- quantile_filtered_test_data_all_patients %>% 
    filter(cell_type_epi == cell_type) %>%
    filter(source_gene %in% gene_set, target_gene %in% gene_set) %>%
    dplyr::select(patient_id, sample_name) %>%
    #unique() %>%
    group_by(patient_id) %>%
    dplyr::summarize('spec_cell_count' = n())
  
  intermediate <- left_join(ncells_perid, count_data) %>% 
    dplyr::mutate('ratio' = spec_cell_count/all_cellcount)
  
  data.frame('patient_id' = intermediate$patient_id, 'clustervalue' = intermediate$ratio, 'cluster' = idx)
}


get_all_clusters <- function(cell_type) {
  all_clusters <- rbindlist(lapply(seq(10), function(x) get_relative_counts(x, cell_type)))
  all_clusters$cell_type_epi <- cell_type
  all_clusters
}

get_all_clusters('Tumor')
get_relative_counts(6, 'Tumor')

cell_types <- quantile_filtered_test_data_all_patients %>% .$cell_type_epi %>% unique()


all_data <- rbindlist(lapply(cell_types, get_all_clusters))
normalized_data <- all_data %>% dplyr::group_by(cluster) %>%
  dplyr::mutate(maxv = max(clustervalue, na.rm=T), minv = min(clustervalue, na.rm=T)) %>%
  dplyr::mutate(normv = (clustervalue)/(maxv))

patient_mutations <- data.frame(patient_id = c('p018', 'p019','p023','p024','p027','p030','p031','p032','p033','p034'),
                                'KRAS' = c(1,0,1,0,0,1,1,1,1,0),
                                'TP53' = c(0,0,1,0,1,0,0,0,0,0),
                                'PIK3CA' = c(0,0,0,0,0,0,1,0,0,0),
                                'EML4' = c(0,0,0,1,0,0,0,0,0,0))


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

ht_opt(legend_border = "black",
       heatmap_border = TRUE)
col_fun = colorRamp2(c(0, 1), c("white", "red"))
col_fun2 = colorRamp2(c(0, 1), c("white", "black"))
#Heatmap(prepare_matrix('Tumor'))
w <- unit(4, "cm")
h <- unit(4, "cm")
mutant_heatmap <- p_mat %>% 
  Heatmap(cluster_columns = F, width=0.3*w, height=h, name = 'Mutation', col = col_fun2, 
          show_heatmap_legend = F, column_names_gp = grid::gpar(fontsize = 9), column_title = 'Mutation')

png('./figures/patient_heatmap.png')
Heatmap(prepare_matrix('Tumor'), col = col_fun, cluster_columns = F, width=w, height=h) + 
  Heatmap(prepare_matrix('Club'), col = col_fun, cluster_columns = F, width=w, height=h)+  
  Heatmap(prepare_matrix('AT1'), col = col_fun, cluster_columns = F, width=w, height=h) + 
  Heatmap(prepare_matrix('AT2'), col = col_fun, cluster_columns = F, width=w, height=h) + 
  Heatmap(prepare_matrix('Ciliated'), col = col_fun, cluster_columns = F, width=w, height=h) + 
  mutant_heatmap
dev.off()

png('./figures/patient_heatmap.png', width=1000, height=200)
Heatmap(prepare_matrix('Tumor'), col = col_fun, cluster_columns = F, column_title ='Tumor', name = 'IpC (scaled)') + 
  Heatmap(prepare_matrix('Club'), col = col_fun, cluster_columns = F, column_title ='Club',show_heatmap_legend = F)+  
  Heatmap(prepare_matrix('AT1'), col = col_fun, cluster_columns = F, column_title ='AT1',show_heatmap_legend = F) + 
  Heatmap(prepare_matrix('AT2'), col = col_fun, cluster_columns = F, column_title ='AT2',show_heatmap_legend = F) + 
  Heatmap(prepare_matrix('Ciliated'), col = col_fun, cluster_columns = F, column_title ='Ciliated',show_heatmap_legend = F) + 
  mutant_heatmap+
  histology_heatmap
dev.off()


  

