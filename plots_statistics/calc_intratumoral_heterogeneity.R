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
library(Matrix)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#select highest 100 (out of highest 500) LRP scores
quant_edges <- 0.8 # 0.8
quant_edges2 <- 1.0
test_data_tumor <- fread(paste('../results/high_values_concat.csv', sep = ""))

############
quantile_data <- test_data_tumor %>% group_by(sample_name) %>% 
  do(data.frame(t(quantile(.$LRP, c(quant_edges, quant_edges2)))))

colnames(quantile_data)[2] <- 'quant'
colnames(quantile_data)[3] <- 'quant2'

quantile_data <- test_data_tumor %>% left_join(quantile_data) %>% filter(LRP>quant, LRP<quant2) %>%
  dplyr::mutate('interaction' = paste0(source_gene, target_gene))

all_interactions <- quantile_data %>% dplyr::select(interaction) %>% unique()
patient_ids <- quantile_data$patient_id %>% unique()
celltypes <- quantile_data$cell_type_epi %>% unique()

##
#######
# prepare interaction vectors (in matrix for all cells) specific for tissue type and patient idx 
prepare_matrix <- function(celltype, idx) {
  data_one_celltype <- quantile_data %>%
    dplyr::mutate('one' = 1) %>%
    dplyr::filter(cell_type_epi == celltype, patient_id == idx)  %>%
    dplyr::select(interaction, sample_name, one) %>%
    pivot_wider(names_from = sample_name, values_from =one)
  data_one_celltype <- left_join(all_interactions, data_one_celltype)
  data_matrix <- data_one_celltype %>%
    dplyr::select(-interaction) %>%
    as.matrix()
  rownames(data_matrix) <- data_one_celltype$interaction
  data_matrix[is.na(data_matrix)] <- 0
  Matrix(data_matrix, sparse = T)
  }


#function to compute intrapatient homogeneity 
get_homogeneity <- function(celltype, idx) {
  data_matrix <- prepare_matrix(celltype, idx)
  if (dim(data_matrix)[2]<5) {return(NULL)}
  product <- t(data_matrix) %*% data_matrix
  diag(product) <- NA
  data.frame('homogeneity' = mean(product, na.rm=T), 'cell_type_epi' = celltype, 'patient_id' = idx)
}

get_all_celltypes <- function(patient_id) {
  print(patient_id)
  rbindlist(lapply(celltypes, function(celltype) get_homogeneity(celltype, patient_id)))
  }


all_results <- rbindlist(lapply(patient_ids, get_all_celltypes))
all_results$comparison <- 'intrapatient'

png('./figures/intratumoral_heterogeneity.png')
ggplot(all_results, aes(x = cell_type_epi, y = homogeneity)) + 
  geom_boxplot()  +
  xlab('Cell Type') + 
  ylab('Homogeneity') +
  theme_bw()
dev.off()


####
#between patients
####

between_patients <- function(celltype) {
  prepare_matrix <- function(idx) {
    data_one_celltype <- quantile_data %>%
      dplyr::mutate('one' = 1) %>%
      dplyr::filter(cell_type_epi == celltype, patient_id == idx)  %>%
      dplyr::select(interaction, sample_name, one) %>%
      pivot_wider(names_from = sample_name, values_from =one)
    data_one_celltype <- left_join(all_interactions, data_one_celltype)
    data_matrix <- data_one_celltype %>%
      dplyr::select(-interaction) %>%
      as.matrix()
    rownames(data_matrix) <- data_one_celltype$interaction
    data_matrix[is.na(data_matrix)] <- 0
    Matrix(data_matrix, sparse=T)
  }
  
  
  get_homogeneity <- function(idx1, idx2) {
    if (idx1<=idx2) return(NULL)
    data_matrix1 <- prepare_matrix(idx1) %>% t()
    data_matrix2 <- prepare_matrix(idx2)
    print(dim(data_matrix1))
    print(dim(data_matrix2))
    if( (dim(data_matrix1)[1]<5) | (dim(data_matrix2)[2]<5)) {return(NULL)}
    product <- data_matrix1 %*% data_matrix2
    #diag(product) <- 0
    data.frame('homogeneity' = mean(product), 'idx1' = idx1, 'idx2' = idx2)
  }
  
  patient_ids <- quantile_data$patient_id %>% unique()
  
  get_all_celltypes <- function(idx2) {
    print(idx2)
    rbindlist(lapply(patient_ids, function(idx1) get_homogeneity(idx1, idx2)))
  }
  
  
  inter_patient_results <- rbindlist(lapply(patient_ids, get_all_celltypes))
  inter_patient_results$celltype_epi = celltype
  inter_patient_results
  
}

inter_patient_results <- rbindlist(lapply(celltypes, between_patients))
inter_patient_results$comparison <- 'interpatient'
inter_patient_results$patient_id <- paste0(inter_patient_results$idx1, inter_patient_results$idx2)
inter_patient_results$cell_type_epi <- inter_patient_results$celltype_epi

png('./figures/interpatient_heterogeneity.png')
ggplot(inter_patient_results, aes(y = homogeneity,x=celltype_epi)) + 
  geom_boxplot() +
  xlab('Cell Type') + 
  ylab('Homogeneity') + 
  theme_bw()
dev.off()
        
#####
intra_and_inter <- rbind(all_results, inter_patient_results[,-c('idx1', 'idx2', 'celltype_epi')])
intra_and_inter$comparison <- as.factor(intra_and_inter$comparison)
#levels(intra_and_inter$comparison) <- c('Between patients', 'Within patients')
png('./figures/inter_and_intra_heterogeneity.png', width=2000, height=2000,  res=300)
ggplot(intra_and_inter, aes(y = homogeneity,x=cell_type_epi)) +
  facet_wrap(~comparison) + 
  geom_boxplot(width = 0.2, outlier.size = 1.0) +
  xlab('Cell Type') + 
  ylab('Network similarity') + 
  theme_minimal() + 
  theme(panel.spacing = unit(2, "lines"),
        text = element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1))
  
dev.off()

median_and_IQR <- intra_and_inter  %>% 
  #filter(homogeneity!=0) %>% 
  group_by(comparison, cell_type_epi) %>%
  dplyr::summarize(med = median(homogeneity), iqr = IQR(homogeneity))

tumor_inter<- inter_patient_results %>% dplyr::filter(cell_type_epi == 'Tumor')
normal_inter<- inter_patient_results %>% dplyr::filter(cell_type_epi != 'Tumor')
wilcox.test(tumor_inter$homogeneity, normal_inter$homogeneity)

tumor_intra<- all_results %>% dplyr::filter(cell_type_epi == 'Tumor')
normal_intra<- all_results %>% dplyr::filter(cell_type_epi != 'Tumor')
wilcox.test(tumor_intra$homogeneity, normal_intra$homogeneity)
median(normal_intra$homogeneity)
median(tumor_intra$homogeneity)

#############
#compute shared interactions between Tumor cells and other cells
#############

compare_to_tumor <- function(idx) {
  prepare_matrix <- function(celltype) {
    
    data_one_celltype <- quantile_data %>%
      dplyr::mutate('one' = 1) %>%
      dplyr::filter(cell_type_epi == celltype, patient_id == idx)  %>%
      dplyr::select(interaction, sample_name, one) %>%
      pivot_wider(names_from = sample_name, values_from =one)
    
    data_one_celltype <- left_join(all_interactions, data_one_celltype)
    
    data_matrix <- data_one_celltype %>%
      dplyr::select(-interaction) %>%
      as.matrix()
    rownames(data_matrix) <- data_one_celltype$interaction
    data_matrix[is.na(data_matrix)] <- 0
    #data_matrix <- as.integer(data_matrix)
    mode(data_matrix) <- 'integer'
    print(class(data_matrix))
    Matrix(data_matrix, sparse = T)
  }
  
  data_matrix1 <- prepare_matrix('Tumor') %>% t() 
  
  get_homogeneity <- function(cellt) {
    data_matrix2 <- prepare_matrix(cellt)
    print(dim(data_matrix1))
    print(dim(data_matrix2))
    if( (dim(data_matrix1)[1]<5) | (dim(data_matrix2)[2]<5)) {return(NULL)}
    product <- data_matrix1 %*% data_matrix2
    #diag(product) <- 0
    df <- data.frame('homogeneity' = mean(product), 'idx' = idx, 'celltype' = cellt)
    print(df)
    df
    
  }
  
  celltypes <- quantile_data$cell_type_epi %>% unique()
  celltypes_wo_tumor <- celltypes[celltypes != 'Tumor']
  print(celltypes_wo_tumor)
    
  
  get_all_celltypes <- function(idx) {
    print(idx)
    rbindlist(lapply(celltypes_wo_tumor, get_homogeneity))
  }
  
  
  inter_patient_results <- get_all_celltypes(idx)
  inter_patient_results
  
}

comparison_to_tumor <- rbindlist(lapply(patient_ids,compare_to_tumor))

png('./figures/compare_to_tumor_interactions.png', width=800, height=800,  res=150)
ggplot(comparison_to_tumor, aes(y = homogeneity,x=celltype)) + 
  geom_boxplot() +
  xlab('Cell Type') + 
  ylab('Average shared interactions') + 
  theme_bw()
dev.off()           

comparison_to_tumor %>%
  group_by(celltype) %>%
  dplyr::summarize('Median' = median(homogeneity), 'IQR' = IQR(homogeneity), 'Mean' = mean(homogeneity))

patient_mutations <- data.frame(patient_id = c('p018', 'p019','p023','p024','p027','p030','p031','p032','p033','p034'),
                                'KRAS' = c(1,0,1,0,0,1,1,1,1,0),
                                'TP53' = c(0,0,1,0,1,0,0,0,0,0),
                                'PIK3CA' = c(0,0,0,0,0,0,1,0,0,0),
                                'EML4' = c(0,0,0,1,0,0,0,0,0,0))

patient_histology <- data.frame(patient_id = c('p018', 'p019','p023','p024','p027','p030','p031','p032','p033','p034'),
                                Histology = as.factor(c('acinar', 'acinar', 'solid', 'acinar', 'solid', 'mucinuous', 'acinar', 'lepidic', 'papillary', 'sarcomatoid')))



comparison_to_tumor <- comparison_to_tumor %>% left_join(patient_histology, by =c('idx' = 'patient_id')) %>%
  left_join(patient_mutations, by =c('idx' = 'patient_id'))

ggplot(comparison_to_tumor, aes(y = homogeneity,x=celltype, color = KRAS)) + 
  geom_point() +
  geom_line(aes(group = idx)) + 
  xlab('Cell Type') + 
  ylab('Average shared interactions') + 
  theme_bw()

