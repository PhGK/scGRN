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
test_data_tumor <- fread(paste('../results/mean_interactions_all_types.csv', sep = "")) %>% filter(source_gene>target_gene)

############
quantile_data <- test_data_tumor# %>% group_by(sample_name) %>% do(data.frame(t(quantile(.$LRP, c(quant_edges)))))
colnames(quantile_data)[2] <- 'quant'
quantile_data <- test_data_tumor %>% left_join(quantile_data) %>% filter(LRP>quant)


######################################################
expression_data <- fread('../data/epi_top2000.csv')
boolean <- ((colnames(expression_data) %in% quantile_data$source_gene) | (colnames(expression_data) %in% quantile_data$target_gene))
length(boolean)

selected_gene_data <- expression_data %>% dplyr::select(colnames(expression_data)[boolean])
selected_matrix <- selected_gene_data %>% as.matrix()

##
selected_matrix[selected_matrix==0] <- NA
corrvalues <- cor(selected_matrix)#, use = "pairwise.complete.obs")
##
#corrvalues <- cor(selected_matrix, method = 'spearman')
#corrvalues[is.na(corrvalues)] <- 0
corrvalues %>% is.na %>% sum()
corrvalues_frame <- as.data.frame(corrvalues)
corrvalues_frame$source_gene <- rownames(corrvalues)

corr_long <- corrvalues_frame %>% 
  as.data.frame() %>%
  pivot_longer(!source_gene, names_to = 'target_gene', values_to = 'corr_value') %>%
  filter(target_gene < source_gene) %>%
  dplyr::mutate('corr_value' = abs(corr_value))

##################
all_data <- inner_join(quantile_data, corr_long, by=c('source_gene', 'target_gene')) %>% 
  dplyr::select(-c('V1')) #%>%
  dplyr::filter(cell_type_epi== 'Tumor')

result <- cor.test(all_data$meanLRP, all_data$corr_value, method =c('pearson'))
result$p.value
result$estimate
######################################################

LRP_frame_for_graph <- all_data %>% dplyr::select(meanLRP, source_gene, target_gene) %>%
  arrange(meanLRP, decreasing=T) %>% .[1:500,] %>%
  dplyr::select('from' = source_gene, 'to' = target_gene) %>% as.matrix()

LRP_graph <- LRP_frame_for_graph %>% graph_from_edgelist(directed=F)
plot(LRP_graph, vertex.size=1, vertex.label.cex = 0.8)

Corr_frame_for_graph <- all_data %>% dplyr::select(corr_value, source_gene, target_gene) %>%
  arrange(corr_value, decreasing=T) %>% .[1:500,] %>%
  dplyr::select('from' = source_gene, 'to' = target_gene) %>% as.matrix()

Corr_graph <- Corr_frame_for_graph %>% graph_from_edgelist(directed=F)
plot(Corr_graph, vertex.size=1, vertex.label.cex = 0.8)

ggplot(all_data, aes(x=meanLRP, y = corr_value)) + geom_point()


#####################################
#comparison with reactome
gene_names <- c(all_data$source_gene, all_data$target_gene) %>% unique()
gene_names %>% length()
gene_names %>% write.csv('/mnt/scratch2/mlprot/Projekte/singlecell/data/gene_names.csv', col.names=FALSE, row.names=F)

Reactome_data <- fread('/mnt/scratch2/mlprot/Projekte/singlecell/data/Reactome_edgelist.sif', header = F) #%>% as.character()

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
  filter((counts>14000)&(counts<17000)) %>%
  dplyr::mutate('corrLRP' = meanLRP * corr_value) %>%
  #dplyr::filter(source_gene != 'MALAT1', target_gene !='MALAT1') %>%
  dplyr::filter(source_gene != 'meanpercell', target_gene !='meanpercell') 
    


all_data_ROC$Reactome_value[is.na(all_data_ROC$Reactome_value)] <- 0

AUCvalueLRP <- roc(all_data_ROC$Reactome_value, all_data_ROC$meanLRP)
AUCvalueLRP$auc
plot( roc(all_data_ROC$Reactome_value, all_data_ROC$meanLRP))


AUCvaluecorr <- roc(all_data_ROC$Reactome_value, all_data_ROC$corr_value)
AUCvaluecorr$auc
plot(roc(all_data_ROC$Reactome_value, all_data_ROC$corr_value))

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

AUC_results <- rbindlist(lapply(seq(1,164000, by=300), estimate_AUC))
library(scales)
log10_rev_trans <- trans_new(
  "log10_rev",
  function(x) log10(rev(x)),
  function(x) rev(10 ^ (x)),
  log_breaks(10),
  domain = c(1e-100, Inf)
)

png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/AUC_REACTOME.png', width=1000, height=1000, res=200)
ggplot(AUC_results, aes(x=counts, y= AUC))+ 
  geom_point()  + 
  #scale_x_continuous(trans=log10_rev_trans) +
  xlab('Minimum Cell Count per Interaction')+
  theme_bw()
dev.off()

AUC_results_long <- AUC_results %>%
  dplyr::mutate(AUC = as.double(AUC), 'ninteractions' = log10(ninteractions)) %>%
  pivot_longer(!counts, values_to='values', names_to = 'ys')


png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/AUC_REACTOME.png', width=1000, height=1000, res=200)
ggplot(AUC_results_long, aes(x=counts, y= values))+ 
  geom_point()  + 
  facet_wrap(~ ys , scale="free", nrow=2)+
  #scale_x_continuous(trans=log10_rev_trans) +
  xlab('Minimum Cell Count per Interaction')+
  ylim(0.5,1)+
  theme_bw()
dev.off()
######
##look for relevant genes
######
all_data_frequent_network <- all_data %>% left_join(Reactome_undirected) %>%
  filter(counts>10000) %>%
  dplyr::filter(source_gene != 'meanpercell', target_gene !='meanpercell') %>%
  dplyr::mutate('interaction' = paste0(source_gene,target_gene))

all_data_frequent_network$Reactome_value[is.na(all_data_frequent_network$Reactome_value)] <- 0

(roc(all_data_frequent_network$Reactome_value, all_data_frequent_network$meanLRP))
genes <- c(all_data_frequent_network$source_gene, all_data_frequent_network$target_gene) %>% unique()
write.csv(genes, '/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/use_data/frequentnetwork.csv', row.names = F)
###############################################
###############################################
frequent_network_LRP <- fread('/mnt/scratch2/mlprot/Projekte/singlecell/data/frequent_network_LRP.csv')
summary_network <- frequent_network_LRP %>%
  group_by(source_gene, target_gene) %>% 
  dplyr::mutate('interaction' = paste0(source_gene,target_gene)) %>%
  dplyr::filter(inpv!=0, tinpv!=0, interaction %in% all_data_frequent_network$interaction) %>%
  dplyr::summarize('meanLRP' = mean(LRP))

summary_REACTOME <- summary_network %>% 
  dplyr::select(meanLRP, source_gene, target_gene) %>%
  left_join(Reactome_undirected)

summary_REACTOME$Reactome_value[is.na(summary_REACTOME$Reactome_value)] <-0

roc(summary_REACTOME$Reactome_value, summary_REACTOME$meanLRP)

get_single_AUCs <- function(n) {
  #cell_type = 'Tumor'
sample_names <- frequent_network_LRP %>% 
  #dplyr::filter(cell_type_epi != cell_type) %>% 
  .$cell_id %>% unique()
sampled_names <- sample(sample_names, n, replace=F)

one_network <-  frequent_network_LRP %>%
  #dplyr::filter(cell_type_epi != cell_type) %>%
  group_by(source_gene, target_gene) %>%
  dplyr::mutate('interaction' = paste0(source_gene,target_gene)) %>%
  dplyr::filter(inpv!=0, tinpv!=0, interaction %in% all_data_frequent_network$interaction, sample_name %in% sampled_names) %>%
  dplyr::summarize('meanLRP' = mean(LRP))


one_REACTOME <- one_network %>% 
  dplyr::select(meanLRP, source_gene, target_gene) %>%
  left_join(Reactome_undirected)

one_REACTOME$Reactome_value[is.na(one_REACTOME$Reactome_value)] <-0

data.frame('AUC' = as.numeric(roc(one_REACTOME$Reactome_value, one_REACTOME$meanLRP)$auc, direction = '<'))
}

AUCs1 <- lapply(rep(1,10), get_single_AUCs)
AUCs10 <- lapply(rep(10,10), get_single_AUCs)
AUCs100 <- lapply(rep(100,10), get_single_AUCs)
AUCs1000 <- lapply(rep(1000,10), get_single_AUCs)
AUCs10000 <- lapply(rep(10000,10), get_single_AUCs)

a1 <-cbind(rbindlist(AUCs1),1)
a2 <-cbind(rbindlist(AUCs10),10)
a3 <-cbind(rbindlist(AUCs100),100)
a4 <-cbind(rbindlist(AUCs1000),1000)
a5 <-cbind(rbindlist(AUCs10000),10000)

frame <- data.frame(rbind(a1,a2,a3,a4,a5))
png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/AUC_frequent_network.png', width=1000, height=1000, res=200)
ggplot(frame, aes(x = V2, y=AUC, group=as.factor(V2))) + 
  geom_boxplot() + 
  scale_x_continuous(trans='log10')+
  xlab('Number of predicted interactions') +
  ylim(0.5,0.8)+
  theme_bw()
  
dev.off()
###############################################
###############################################
k = 100

ncorrectall <- all_data_ROC$Reactome_value %>% sum()
nfalseall <- (1-all_data_ROC$Reactome_value) %>% sum()

phyper(meanLRP_correct-1, ncorrectall, nfalseall, k, lower.tail = F, log.p = FALSE)
phyper(corr_correct-1, ncorrectall, nfalseall, k, lower.tail = F, log.p = FALSE)
phyper(LRP_corr_correct-1, ncorrectall, nfalseall, k, lower.tail = F, log.p = FALSE)

####################################
####################################
single_cell <- fread('/mnt/scratch2/mlprot/Projekte/singlecell/data/high_values_concat.csv') %>%
  filter(source_gene != 'meanpercell', target_gene != ',meanpercell', source_gene > target_gene)


########
#correct for dropouts
expr_data <- fread('/mnt/scratch2/mlprot/Projekte/singlecell/data/epi_top2000.csv')[,3:2002]
cell_type_epi <- fread('/mnt/scratch2/mlprot/Projekte/singlecell/data/epi_top2000.csv')$cell_type_epi
#zerorate_source <- data.frame('source_gene' = colnames(expr_data), 'zerorate_source' = colMeans(expr_data==0))
#zerorate_target <- data.frame('target_gene' = colnames(expr_data), 'zerorate_target' = colMeans(expr_data==0))

cell_frame <- data.frame('cell_type_epi' = cell_type_epi, expr_data)

zerorate_source <-  aggregate(cell_frame[,-(1)]==0, list(cell_frame$cell_type_epi), mean) %>%
  pivot_longer(!Group.1, names_to ='source_gene', values_to = 'zerorate_source')

colnames(zerorate_source)[1] <- 'cell_type_epi'
zerorate_target <-zerorate_source
colnames(zerorate_target)[2:3] <- c('target_gene', 'zerorate_target') 

resampled_data <- left_join(single_cell, zerorate_source) %>%
  left_join(zerorate_target) %>%
  mutate(random = runif(dim(.)[1]), zerorate = zerorate_target * zerorate_source) %>%
  mutate(selected_rows = zerorate>(random-0.01))
resampled_data <- resampled_data %>% filter(selected_rows)

########

resampled_data_counts <-resampled_data%>% 
  dplyr::filter(as.character(cell_type_epi) == 'Tumor') %>%
  dplyr::select(source_gene, target_gene) %>%
  dplyr::group_by(source_gene, target_gene) %>%
  dplyr::summarize('counts' = n())

single_cell_ROC <- resampled_data_counts %>% left_join(Reactome_undirected)
single_cell_ROC$Reactome_value[is.na(single_cell_ROC$Reactome_value)] <- 0

single_cell_ROC$Reactome_value %>% mean()

