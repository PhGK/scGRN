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

# - build modules
# - compare modules between cell types and patients (Heatmap and chi2 test)
# - show overlap between cell types



quant_edges <- 0.8
test_data_tumor <- fread(paste('../usedata/high_values_concat.csv', sep = ""))

filtered_test_data <- test_data_tumor %>% filter(target_gene<source_gene) %>% dplyr::select(-c(V1, cell_id)) %>%
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
all_graph <- resampled_data %>% #quantile_filtered_test_data_all_patients %>%
  dplyr::select(source_gene, target_gene) %>% 
  as.matrix() %>%
  graph_from_edgelist() %>%
  as.undirected()

E(all_graph)$weight <-1
all_graph <- simplify(all_graph, edge.attr.comb =list(weight='sum'))
all_graph <- delete_edges(all_graph,E(all_graph)[E(all_graph)$weight<100])
plot(all_graph, vertex.label=NA, vertex.size=0.1)
E(all_graph)$weight <-1


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

get_relative_counts <- function(idx, cell_type) {
    gene_set <- degreeframe %>% filter(N>3, cluster==idx) %>% .$name
  
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
    unique() %>%
    group_by(patient_id) %>%
    dplyr::summarize('spec_cell_count' = n())

    intermediate <- left_join(ncells_perid, count_data) %>% 
      dplyr::mutate('ratio' = spec_cell_count/all_cellcount)
    
    data.frame('patient_id' = intermediate$patient_id, 'clustervalue' = intermediate$ratio)
}

get_all_clusters <- function(cell_type) {
all_clusters <- get_relative_counts(1, cell_type)
colnames(all_clusters)[2] <- paste0('cluster_', 1)
for (i in seq(2,13)) {
  all_clusters <- left_join(all_clusters, get_relative_counts(i, cell_type))
  colnames(all_clusters)[i+1] <- paste0('cluster_', i)
}
all_clusters
}



patient_mutations <- data.frame(patient_id = c('p018', 'p019','p023','p024','p027','p030','p031','p032','p033','p034'),
                                'KRAS' = c(1,0,1,0,0,1,1,1,1,0),
                                'TP53' = c(0,0,1,0,1,0,0,0,0,0),
                                'PIK3CA' = c(0,0,0,0,0,0,1,0,0,0),
                                'EML4' = c(0,0,0,1,0,0,0,0,0,0))

p_mat <- patient_mutations %>% dplyr::select(-patient_id) %>% as.matrix()
rownames(p_mat) <- patient_mutations$patient_id
mutant_heatmap <- p_mat %>% Heatmap()

prepare_matrix <- function(cell_type) {
  all_clusters_patient <- get_all_clusters(cell_type)
  rownames(all_clusters_patient) <- all_clusters_patient$patient_id
  heatmap_matrix <- all_clusters_patient %>% dplyr::select(-patient_id) %>% as.matrix()
  heatmap_matrix <- log(1+heatmap_matrix)
  #heatmap_matrix <- apply(heatmap_matrix, 2, function(x) (x-min(x))/max(x-min(x)))
  mask <- apply(heatmap_matrix,2, function(x) all(is.na(x)))
  heatmap_matrix <- heatmap_matrix[,!mask]
  heatmap_matrix
}

ht_opt(legend_border = "black",
       heatmap_border = TRUE)
col_fun = colorRamp2(c(0, 1), c("white", "red"))
#Heatmap(prepare_matrix('Tumor'))
Heatmap(prepare_matrix('Tumor'), col = col_fun, cluster_columns = F) + 
  Heatmap(prepare_matrix('Club'), col = col_fun, cluster_columns = F)+  
  Heatmap(prepare_matrix('AT2'), col = col_fun, cluster_columns = F) + 
  mutant_heatmap


####
#Heatmap according to celltype
###

get_relative_counts_celltype <- function(idx) {
  
  #gene_set <- V(all_graph)$name[modules$membership==idx]
  gene_set <- degreeframe %>% filter(N>1, cluster==idx) %>% .$name
  
  ncells_perid <- quantile_filtered_test_data_all_patients %>% 
    dplyr::select(cell_type_epi, sample_name) %>%
    unique() %>%
    group_by(cell_type_epi) %>%
    dplyr::summarize('all_cellcount' = n())
  
  count_data <- quantile_filtered_test_data_all_patients %>% 
    filter(source_gene %in% gene_set, target_gene %in% gene_set) %>%
    dplyr::select(cell_type_epi, sample_name) %>%
    unique() %>%
    group_by(cell_type_epi) %>%
    dplyr::summarize('spec_cell_count' = n())
  
  intermediate <- left_join(ncells_perid, count_data) %>% 
    dplyr::mutate('ratio' = spec_cell_count/all_cellcount)
  
  data.frame('cell_type_epi' = intermediate$cell_type_epi, 'clustervalue' = intermediate$ratio)
}

all_clusters <- get_relative_counts_celltype(1)
colnames(all_clusters)[2] <- paste0('C', 1)
for (i in seq(2, dim(cluster_sizes)[1])) {
  all_clusters <- left_join(all_clusters, get_relative_counts_celltype(i))
  colnames(all_clusters)[i+1] <- paste0('C', i)
}
rownames(all_clusters)<- all_clusters$cell_type_epi
heatmap_matrix <- all_clusters %>% dplyr::select(-cell_type_epi) %>% as.matrix 
heatmap_matrix <- apply(heatmap_matrix, 2, function(x) (x-min(x, na.rm=T))/max(x-min(x, na.rm=T), na.rm=T))
#heatmap_matrix <- log(1+heatmap_matrix)

mask <- apply(heatmap_matrix,2,function(x) !all(is.na(x)))
heatmap_matrix <- heatmap_matrix[,mask]

col_fun = colorRamp2(c(0, 1), c("white", "red"))

Heatmap(heatmap_matrix, col=col_fun, name = 'IpC (scaled)')

png('./figures/celltype_heatmap.png', width=1000, height=600, res = 200)
Heatmap(heatmap_matrix, col=col_fun, name = 'IpC (scaled)')
dev.off()

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


cell_counts2 <- quantile_filtered_test_data_all_patients$sample_name %>% unique %>% length()
quantification_all <- quantile_filtered_test_data_all_patients %>% 
  group_by(source_gene, target_gene) %>%
  dplyr::summarize('interactioncount' = n()) %>%
  dplyr::mutate('relative_count' = interactioncount /cell_counts2)

####
#cluster
current_cluster <- degreeframe %>% filter(degree>10, cluster==10)
quantification_cluster_all <- quantification_all %>% dplyr::filter((source_gene %in% current_cluster$name)&  (target_gene %in% current_cluster$name))
quantification_cluster <- quantification %>% dplyr::filter((source_gene %in% current_cluster$name)&  (target_gene %in% current_cluster$name))

quantification_cluster_average <- quantile_filtered_test_data_all_patients %>% 
  group_by(cell_type_epi) %>%
  dplyr::filter((source_gene %in% current_cluster$name) &  (target_gene %in% current_cluster$name)) %>%
  dplyr::summarize('interactioncount' = n()) %>%
  left_join(cell_counts) %>%
  dplyr::mutate('relativecount' = interactioncount/cell_count)

write.csv(current_cluster$name %>% unique() %>% as.character(), '../usedata/cluster_genes.csv', 
          row.names = F, quote = F)

numberofgenes <- c(quantification_cluster$source_gene, quantification_cluster$target_gene) %>% unique()
length(numberofgenes)
######################
#################
#compute statistical comparison for pair of genes 
#C1
#first find prior distribution

gene_expression <- fread('../data/epi_top2000.csv')
get_all_interactions <- function(gene1, gene2) {
  expr_data <- gene_expression %>% dplyr::select(gene1, gene2, cell_type_epi) %>%
    dplyr::mutate(inter_exists = (.[,1]!=0)&(.[,2]!=0)) %>%
    dplyr::group_by(cell_type_epi) %>%
    dplyr::summarize('n_interactions' = sum(inter_exists))
  expr_data$gene1 <- gene1
  expr_data$gene2 <- gene2
  expr_data
}

current_cluster <- degreeframe %>% filter(degree>10, cluster==1)

all_genes <- list()
i = 1
for (gene1 in current_cluster$name) {
  for (gene2 in current_cluster$name) {
    get_all_interactions(gene1, gene2)
    all_genes[[i]] <- get_all_interactions(gene1, gene2)
    i <- i + 1
  }
}

all_genes_frame <- rbindlist(all_genes) %>% filter(gene1>gene2) %>%
  group_by(cell_type_epi) %>%
  dplyr::summarize('n_interactions' = sum(n_interactions))

for_chi2 <- quantile_filtered_test_data_all_patients %>% 
  group_by(cell_type_epi) %>%
  dplyr::filter((source_gene %in% current_cluster$name) &  (target_gene %in% current_cluster$name)) %>%
  dplyr::summarize('interactioncount' = n()) %>%
  left_join(all_genes_frame)

for_chi2_matrix <- for_chi2 %>% dplyr::select(-cell_type_epi)%>% as.matrix()
rownames(for_chi2_matrix) <- for_chi2$cell_type_epi


TumorClub_vs_AT1AT2Ciliated <- for_chi2 %>%
  dplyr::mutate('cell_type_group' = ifelse(cell_type_epi %in% c('Tumor', 'Club'),'TumorClub', 'AT1AT2Ciliated')) %>%
  group_by(cell_type_group) %>%
  dplyr::summarize('interactioncount' = sum(interactioncount), n_interactions = sum(n_interactions))


TumorClub_vs_AT1AT2Ciliated_t <- TumorClub_vs_AT1AT2Ciliated %>% dplyr::select(-cell_type_group) %>% t() 
colnames(TumorClub_vs_AT1AT2Ciliated_t) <- TumorClub_vs_AT1AT2Ciliated$cell_type_group
chisq.test(TumorClub_vs_AT1AT2Ciliated_t)$p.value

chisq.test(for_chi2_matrix)$p.value
chisq.test(for_chi2_matrix)$observed

TumorClubAT2_vs_AT1Ciliated <- for_chi2 %>%
  dplyr::mutate('cell_type_group' = ifelse(cell_type_epi %in% c('Tumor', 'Club', 'AT2'),'TumorClub', 'AT1AT2Ciliated')) %>%
  group_by(cell_type_group) %>%
  dplyr::summarize('interactioncount' = sum(interactioncount), n_interactions = sum(n_interactions))

TumorClubAT2_vs_AT1Ciliated_t <- TumorClubAT2_vs_AT1Ciliated %>% dplyr::select(-cell_type_group) %>% t() 
colnames(TumorClubAT2_vs_AT1Ciliated_t) <- TumorClubAT2_vs_AT1Ciliated$cell_type_group
chisq.test(TumorClubAT2_vs_AT1Ciliated_t)$p.value

Ciliated_vs_rest <- for_chi2 %>%
  dplyr::mutate('cell_type_group' = ifelse(cell_type_epi %in% c('Ciliated'),'Ciliated', 'rest')) %>%
  group_by(cell_type_group) %>%
  dplyr::summarize('interactioncount' = sum(interactioncount), n_interactions = sum(n_interactions))

Ciliated_vs_rest_t <- Ciliated_vs_rest %>% dplyr::select(-cell_type_group) %>% t() 
colnames(Ciliated_vs_rest_t) <- Ciliated_vs_rest$cell_type_group
chisq.test(Ciliated_vs_rest_t)$p.value

Tumor_vs_rest <- for_chi2 %>%
  dplyr::mutate('cell_type_group' = ifelse(cell_type_epi %in% c('Tumor'),'Tumor', 'AT1AT2Ciliated')) %>%
  group_by(cell_type_group) %>%
  dplyr::summarize('interactioncount' = sum(interactioncount), n_interactions = sum(n_interactions))


Tumor_vs_rest_t <- Tumor_vs_rest %>% dplyr::select(-cell_type_group) %>% t() 
colnames(Tumor_vs_rest_t) <- Tumor_vs_rest$cell_type_group
chisq.test(Tumor_vs_rest_t)$p.value

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

######################################################
#cell counts, validity
unique_cells <- test_data_tumor %>% 
  group_by(cell_id, cell_type_epi) %>% 
  dplyr::summarize('meanLRP' = mean(LRP))

test_data_tumor %>% select(-c(LRP, source_gene, target_gene, V1, 'Unnamed: 0')) %>% unique() %>% dim()

(test_data_tumor %>% 
    dplyr::select(-c(LRP, source_gene, target_gene)) %>% 
    unique %>% dim())[1] /500 

count_data <- fread('../data/epi_top2000.csv')[,-c(5:2000)]
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

