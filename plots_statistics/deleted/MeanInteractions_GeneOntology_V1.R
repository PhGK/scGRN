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

test_data_all <- fread(paste('../results/mean_interactions_all_types.csv', sep = ""))%>%
  dplyr::mutate('cell_type_epi' = 'all')

filtered_test_data <- test_data_all %>% filter(target_gene<source_gene)

###########################
#GO TERMS
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

getGO <- function(go_id) {
  genelist <- getBM(attributes = "hgnc_symbol",
                    filters = "go_parent_term",
                    values = go_id,
                    mart=mart)
  
  print(length(genelist))
  geneframe <- as.data.frame(genelist)
  print(dim(geneframe))
  geneframe$intermediate <-1
  colnames(geneframe)[2]<-as.character(go_id)
  geneframe
}

GOTERMS <- data.frame("description" = c("regulation of apoptotic process",
                                        #"positive regulation of apoptotic process",
                                        #"negative regulation of apoptotic process",
                                        "cell population proliferation",
                                        #"positive cell population proliferation",
                                        #"negative cell population proliferation",
                                        "cellular metabolic process",
                                        "immune system process",
                                        #"cell death",
                                        #"cell division",
                                        #"cell cycle process",
                                        "cell cycle",
                                        #"cell_growth",
                                        "intercellular transport"),
                      "GO_ID" = c("GO:0042981",
                                  #"GO:0043065",
                                  #"GO:0043066",
                                  "GO:0008283",
                                  #"GO:0008284",
                                  #"GO:0008285",
                                  "GO:0044237",
                                  "GO:0002376",
                                  #"GO:0008129",
                                  #"GO:0051301",
                                  #"GO:0022402",
                                  "GO:0007049",
                                  #"GO:0048591",
                                  "GO:0010496"))
#####################################


genelists <- lapply(GOTERMS$GO_ID, getGO)

comb <- full_join(genelists[[1]], genelists[[2]])
comb <- full_join(comb, genelists[[3]])
comb <- full_join(comb, genelists[[4]])
comb <- full_join(comb, genelists[[5]])
comb <- full_join(comb, genelists[[6]])

for (i in seq(3,length(genelists))) {
  comb <- full_join(comb, genelists[[i]])
  }

apply(is.na(comb),1,sum) %>% as.factor() %>% summary()


new_comb <- comb
new_comb[is.na(new_comb)]<-0

which_main_pie <- function(gene_name) { 
  gene_name <- as.character(gene_name)
  print(gene_name)
  color <- c(0,
             new_comb$'GO:0042981'[new_comb$hgnc_symbol==gene_name],
             new_comb$'GO:0008283'[new_comb$hgnc_symbol==gene_name],
             new_comb$'GO:0044237'[new_comb$hgnc_symbol==gene_name],
             new_comb$'GO:0002376'[new_comb$hgnc_symbol==gene_name],
             new_comb$'GO:0007049'[new_comb$hgnc_symbol==gene_name],
             new_comb$'GO:0010496'[new_comb$hgnc_symbol==gene_name])
  if (length(color) == 0) {color <- c(1, 0, 0, 0, 0, 0, 0)}
  if (sum(color)==0 ) {color <- c(1, 0, 0, 0, 0, 0, 0)}
  color
}

which_main_pie <- function(gene_name) { 
  gene_name <- as.character(gene_name)
  print(gene_name)
  color <- c(0,
    new_comb$'GO:0042981'[new_comb$hgnc_symbol==gene_name], 
    new_comb$'GO:0022402'[new_comb$hgnc_symbol==gene_name],
    new_comb$'GO:0044237'[new_comb$hgnc_symbol==gene_name],
    new_comb$'GO:0002376'[new_comb$hgnc_symbol==gene_name],
    new_comb$'GO:0051301'[new_comb$hgnc_symbol==gene_name],
    new_comb$'GO:0007049'[new_comb$hgnc_symbol==gene_name],
    new_comb$'GO:0008283'[new_comb$hgnc_symbol==gene_name],
    new_comb$'GO:0010496'[new_comb$hgnc_symbol==gene_name])
    if (length(color) == 0) {color <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)}
    if (sum(color)==0 ) {color <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)}
    color
}



which_main_pie <- Vectorize(which_main_pie)

##########################
epi_tissue <- filtered_test_data$cell_type_epi %>% unique()
filtered_test_data <- filtered_test_data %>% mutate(meanLRP = LRP) %>% 
  filter(counts>1000) %>%
  filter(source_gene != 'meanpercell', target_gene !='meanpercell')


data_now <-filtered_test_data %>% 
  dplyr::arrange(desc(LRP)) 

nodes <- c(data_now$source_gene, data_now$target_gene) %>% unique()
selflooplist <- cbind("source_gene" = nodes, "target_gene"=nodes)
cutoff <- data_now$LRP[500]
data_now <- data_now %>% dplyr::mutate(edge = ifelse(meanLRP<cutoff,0,1))
edge_frame <- data_now %>% dplyr::select(source_gene, target_gene, edge, LRP)# %>% .[1:1000000,]
edge_list <- edge_frame %>% dplyr::filter(edge==1) %>% dplyr::select(source_gene, target_gene, LRP) #%>% rbind(selflooplist)
colnames(edge_list)[3] <- "weight"
graph <- graph_from_data_frame(edge_list)%>% as.undirected() %>% simplify() 

components(graph)$membership %>% as.factor() %>% summary()
large_c <- which(components(graph)$membership < 10)
large_component = induced_subgraph(graph, large_c)
E(large_component)$value <- E(large_component)$weight
E(large_component)$weight <- 0.1
set.seed(0)
l <- layout_nicely(large_component)
#l <- layout_with_dh(large_component)

a <- V(large_component)
b <- lapply(V(large_component)$name, which_main_pie)
b

border <- lapply(seq(length(b)), function(x) c(1,1,1,1,1,1,1))
border <- list(c('white','white','white','white','white','white','white'))
png(paste('./figures/singlecell_mean_pie_', epi_now,'.png', sep=''), width=1000, height=1000)
plot(large_component, 
     vertex.shape = "pie", #vertex.label=NA, 
     vertex.size=6, vertex.pie= b, vertex.pie.color = list(adjustcolor(c('white', 'red',  'pink', 'green', 'yellow', 'blue', 'brown'), alpha.f=0.4)),
     edge.color='black', layout=l, 
     vertex.pie.border= 'white',#list(adjustcolor(c('blue', 'red', 'green', 'yellow', 'darkblue', 'pink', 'orange'), alpha.f=0.4)),
     #vertex.pie.density=0,
     vertex.pie.lty = 'blank',
     edge.width= 5*log(1.1+E(large_component)$value*1.4),
    edge.color=adjustcolor(1, alpha.f = 1.0))#, xlim=c(0.6,1.2), ylim=c(-0.4,0.4)) 

#legend
dev.off()

################################################
#add bounding box
png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_mean_boxa.png', width=1000, height=1000)
par(mar = c(0,0,0,0), oma = c(0,0,0,0))


plot(large_component, #vertex.label=NA, 
     vertex.size=2, vertex.color= adjustcolor(which_main(V(large_component)$name),alpha.f=0.4),
     edge.color='black', layout=l, 
     edge.width= 5*log(1.1+E(large_component)$value*1.4),
     vertex.frame.color = adjustcolor("white", alpha.f = 0), edge.color=adjustcolor(1, alpha.f = 1.0), xlim=c(x1a,x2a), ylim=c(y1a,y2a)) 
rect(x1a,y1a,x2a,y2a)

dev.off()

png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_mean_boxb.png', width=1000, height=1000)
par(mar = c(0,0,0,0), oma = c(0,0,0,0))

plot(large_component, #vertex.label=NA, 
     vertex.size=2, vertex.color= adjustcolor(which_main(V(large_component)$name),alpha.f=0.4),
     edge.color='black', layout=l, 
     edge.width= 5*log(1.1+E(large_component)$value*1.4),
     vertex.frame.color = adjustcolor("white", alpha.f = 0), edge.color=adjustcolor(1, alpha.f = 1.0), xlim=c(x1b,x2b), ylim=c(y1b,y2b)) 
rect(x1b,y1b,x2b,y2b)

dev.off()

##############################################################
filtered_test_data$cell_type_epi %>% unique()
all_sub_types <- filtered_test_data %>% filter(cell_type_epi!='all') #%>%
  #dplyr::group_by(cell_type_epi, source_gene, target_gene) %>%
  #dplyr::summarize('meanLRP'=mean(mean))

quantile <- all_sub_types %>% group_by(cell_type_epi) %>% do(data.frame(t(quantile(.$mean, c(0.9998)))))
colnames(quantile)[2] <- 'quant'
quantile_filtered_test_data <- all_sub_types %>% left_join(quantile) %>% filter(mean>quant)

#quantile_filtered_test_data <- all_sub_types %>% left_join(quantile) %>% filter(mean>quantile(mean, 0.9999))

quantile_filtered_test_data$cell_type_epi <- as.factor(quantile_filtered_test_data$cell_type_epi)
quantile_filtered_test_data$group <- as.numeric(quantile_filtered_test_data$cell_type_epi)

all_graph <- quantile_filtered_test_data %>% dplyr::select(from =source_gene, to = target_gene, cluster = group) %>%
  graph_from_data_frame(directed=F)


mycolors <- c('pink', 'green', 'darkblue','steelblue',  'orange', 'red', 'yellow')

choose_color <- function(cluster) {
  colors <- mycolors
  colors[cluster]
  
}
choose_color(quantile_filtered_test_data$cell_type_epi[100])

choose_color(E(all_graph)$cluster)
E(all_graph)$color %>% length()
#choose_color <- Vectorize(choose_color)
E(all_graph)$color <- choose_color(E(all_graph)$cluster)
E(all_graph)$color 
E(all_graph)$cluster

png('/mnt/scratch2/mlprot/Projekte/singlecell/plots_statistics/figures/singlecell_all_celltypeepi_all_colored_small.png', width=1500, height=1500)
plot(all_graph, vertex.size=1, edge.width=1, 
     vertex.label=NA,
     vertex.label.cex= 1.0)
legend("topleft",bty = "n",
       legend=levels(quantile_filtered_test_data$cell_type_epi),
       fill=mycolors, border=NA,
       cex = 2.0)

dev.off()
########
for (second_epi_type in quantile_filtered_test_data$cell_type_epi %>% unique()) { 
one <- quantile_filtered_test_data %>% filter(cell_type_epi == c('Tumor')) 
two <- quantile_filtered_test_data %>% filter(cell_type_epi == c(second_epi_type)) 
print(second_epi_type)
print(one %>% filter(source_gene %in% two$source_gene, target_gene %in% two$target_gene) %>% dim() %>% .[1])
}
