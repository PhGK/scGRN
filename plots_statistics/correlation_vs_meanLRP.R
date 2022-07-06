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

quant_edges <- 0.8
average_LRP_scores <- fread(paste('../results/mean_interactions_all_types.csv', sep = "")) %>% filter(source_gene>target_gene)

############

######################################################
expression_data <- fread('../data/epi_top2000.csv')
boolean <- ((colnames(expression_data) %in% average_LRP_scores$source_gene) | (colnames(expression_data) %in% average_LRP_scores$target_gene))
length(boolean)

selected_gene_data <- expression_data %>% dplyr::select(colnames(expression_data)[boolean])
selected_matrix <- selected_gene_data %>% as.matrix()

##
corrvalues <- cor(selected_matrix)#, use = "pairwise.complete.obs")
dim(corrvalues)
corrvalues %>% is.na %>% sum()
corrvalues_frame <- as.data.frame(corrvalues)
corrvalues_frame$source_gene <- rownames(corrvalues)

corr_long <- corrvalues_frame %>% 
  as.data.frame() %>%
  pivot_longer(!source_gene, names_to = 'target_gene', values_to = 'corr_value') %>%
  filter(target_gene < source_gene) %>%
  dplyr::mutate('corr_value' = abs(corr_value))
#################
#partial correlation
library(ppcor)
pcorrvalues <- pcor(selected_matrix)$estimate
dim(pcorrvalues)
rownames(pcorrvalues) <- colnames(selected_matrix)
colnames(pcorrvalues) <- colnames(selected_matrix)

pcorrvalues %>% is.na %>% sum()
pcorrvalues_frame <- as.data.frame(pcorrvalues)
pcorrvalues_frame$source_gene <- rownames(pcorrvalues)

pcorr_long <- pcorrvalues_frame %>% 
  as.data.frame() %>%
  pivot_longer(!source_gene, names_to = 'target_gene', values_to = 'pcorr_value') %>%
  filter(target_gene < source_gene) %>%
  dplyr::mutate('pcorr_value' = abs(pcorr_value))

##################
all_data <- inner_join(average_LRP_scores, corr_long, by=c('source_gene', 'target_gene')) %>% 
  inner_join(pcorr_long, by=c('source_gene', 'target_gene')) %>%
  dplyr::select(-c('V1'))

result <- cor.test(all_data$LRP, all_data$corr_value, method =c('spearman'))
result$p.value
result$estimate
result <- cor.test(all_data$LRP, all_data$pcorr_value, method =c('spearman'))
result$p.value
result$estimate
######################################################

LRP_frame_for_graph <- all_data %>% dplyr::select(LRP, source_gene, target_gene) %>%
  arrange(LRP, decreasing=T) %>% .[1:500,] %>%
  dplyr::select('from' = source_gene, 'to' = target_gene) %>% as.matrix()

LRP_graph <- LRP_frame_for_graph %>% graph_from_edgelist(directed=F)
plot(LRP_graph, vertex.size=1, vertex.label.cex = 0.8)

Corr_frame_for_graph <- all_data %>% dplyr::select(corr_value, source_gene, target_gene) %>%
  arrange(corr_value, decreasing=T) %>% .[1:500,] %>%
  dplyr::select('from' = source_gene, 'to' = target_gene) %>% as.matrix()

Corr_graph <- Corr_frame_for_graph %>% graph_from_edgelist(directed=F)
plot(Corr_graph, vertex.size=1, vertex.label.cex = 0.8)

png('./figures/correlation_vs_LRP.png', height=1000, width=1000, res=100)
ggplot(all_data, aes(x=LRP, y = corr_value)) + geom_point(size=0.2)+theme_bw()
dev.off()

#####################################
#comparison with reactome
gene_names <- c(all_data$source_gene, all_data$target_gene) %>% unique()
gene_names %>% length()
#gene_names %>% write.csv('../data/gene_names.csv', col.names=FALSE, row.names=F)

Reactome_data <- fread('../data/Reactome_edgelist.sif', header = F) #%>% as.character()

Reactome_data <- Reactome_data %>% dplyr::filter(!is.na(str_match(V1, '-')))
Reactome_new <- str_split_fixed(Reactome_data$V1, "\t-\t", 2)
Reactome_undirected <- rbind(Reactome_new, Reactome_new[, c(2,1)]) %>% as.data.frame()
colnames(Reactome_undirected) <- c('source_gene', 'target_gene')
Reactome_undirected$source_gene <- as.character(Reactome_undirected$source_gene)
Reactome_undirected$target_gene <- as.character(Reactome_undirected$target_gene)
Reactome_undirected <- Reactome_undirected %>% dplyr::mutate('Reactome_value'=1) 


#############
library(pROC)
all_data_ROC <- all_data %>% left_join(Reactome_undirected) %>%
  dplyr::mutate('corrLRP' = LRP * corr_value)
    


all_data_ROC$Reactome_value[is.na(all_data_ROC$Reactome_value)] <- 0

AUCvalueLRP <- roc(all_data_ROC$Reactome_value, all_data_ROC$LRP)
AUCvalueLRP$auc
plot( roc(all_data_ROC$Reactome_value, all_data_ROC$LRP))


AUCvaluecorr <- roc(all_data_ROC$Reactome_value, all_data_ROC$corr_value)
AUCvaluecorr$auc
plot(roc(all_data_ROC$Reactome_value, all_data_ROC$corr_value))

AUCvaluepcorr <- roc(all_data_ROC$Reactome_value, all_data_ROC$pcorr_value)
AUCvaluepcorr$auc
plot(roc(all_data_ROC$Reactome_value, all_data_ROC$pcorr_value))

AUCvaluecorr <- roc(all_data_ROC$Reactome_value, all_data_ROC$corrLRP)
AUCvaluecorr$auc
plot(roc(all_data_ROC$Reactome_value, all_data_ROC$corrLRP))
#############################
#############################

estimate_AUC <- function(counts_) {
  all_data_ROC <- all_data %>% left_join(Reactome_undirected) %>%
    filter(counts>counts_) %>%
    dplyr::filter(source_gene != 'meanpercell', target_gene !='meanpercell') 
  
  all_data_ROC$Reactome_value[is.na(all_data_ROC$Reactome_value)] <- 0
  
  AUCvalueLRP <- roc(all_data_ROC$Reactome_value, all_data_ROC$meanLRP)
  data.frame('AUC' = AUCvalueLRP$auc,'counts' = counts_, 'ninteractions'=dim(all_data_ROC)[1])
}


##########
# see n highest LRP values
##########
n = 100
all_data_highest_LRP <- all_data_ROC %>% dplyr::arrange(desc(LRP)) %>%.[1:n,]
correct_LRP <- sum(all_data_highest_LRP$Reactome_value)
false_LRP <- n-correct_LRP
correct_all <- sum(all_data_ROC$Reactome_value)
false_all <-  sum(all_data_ROC$Reactome_value==0)
LRP_mat <- matrix(c(correct_LRP, false_LRP,correct_all, false_all), nrow=2, byrow=T)
chisq.test(LRP_mat)

##########
all_data_highest_correlation <- all_data_ROC %>% dplyr::arrange(desc(corr_value)) %>%.[1:n,]
correct_corr <- sum(all_data_highest_correlation$Reactome_value)
false_corr<- n-correct_corr
corr_mat <- matrix(c(correct_corr, false_corr,correct_all, false_all), nrow=2, byrow=T)
chisq.test(corr_mat)

##
all_data_highest_pcorrelation <- all_data_ROC %>% dplyr::arrange(desc(pcorr_value)) %>%.[1:n,]
correct_pcorr <- sum(all_data_highest_pcorrelation$Reactome_value)
false_pcorr<- n-correct_pcorr

chisq.test(all_data_highest_LRP$LRP, all_data_highest_LRP$Reactome_value)

##
correct_LRP
correct_corr
correct_pcorr

