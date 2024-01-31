rm(list = ls())

library(gplots)
library(RColorBrewer)
library(ape)
library(stats)
library(tidyverse)
library(reshape)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)

#####################################################
######## Pairwise SNP clustering #########
#####################################################

### load DNA alignment
all_ali <- read.dna('input/gubbins_Staph_all.filtered_polymorphic_sites.fasta', format='fasta')
### Calculate pairwise distance btw SNPS
pw_dist_all <- dist.gene(all_ali)
m_all <- as.matrix(pw_dist_all)
m_all <- melt(m_all)[melt(upper.tri(m_all))$value,]
names(m_all) <- c("c1", "c2", "distance")
m_all
#write.csv(m_all, "results/pairwise_distances_all.csv", row.names = F)

#cluster analysis
all_cluster <- hclust(pw_dist_all) #perform hierarchical clustering
all_clus14<-as.data.frame(cutree(all_cluster, h=14))  %>% tibble::rownames_to_column()
colnames(all_clus14)<- c("Isolate","SNP14")
all_SNP14_counts <- all_clus14 %>% dplyr::count(SNP14, sort=TRUE) %>% filter(n>1)
#write.csv(all_SNP14_counts, "results/transmission_cluster_sizes_SNP14.csv", row.names = F)
#list of isolates with transmission cluster id
all_SNP14_isolates<-all_clus14 %>% filter(SNP14 %in% c("4","6","15","17","18","19","22","23","34", "39", "40","45","50","51","58", "63"))

