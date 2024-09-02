setwd("D:/master_project_networks/workflows/otu_table")

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

Pvals <- read.table(file = "network_output_fastspar/wl_nw_fp10_output/pvalues.tsv", sep="\t", header=T, row.names=1,  comment.char = "")
colnames(Pvals) <- gsub("^X", "", colnames(Pvals))
Pvals1 <- as.matrix(Pvals)
Pvals1[upper.tri(Pvals, diag=TRUE)]<- NA
Pvals2 <- as.data.frame(as.table(Pvals1))
colnames(Pvals2)<-c("Node1","Node2","Pvalue")

head(Pvals2)

Cors <- read.table("network_output_fastspar/wl_nw_fp10_output/median_correlation.tsv", sep="\t", header=T, row.names=1 , comment.char = "")
colnames(Cors) <- gsub("^X", "", colnames(Cors))
Cors1 <- as.matrix(Cors)
Cors1[upper.tri(Cors1, diag=TRUE)]<-NA
Cors2 <- as.data.frame(as.table(Cors1))
colnames(Cors2)<-c("Node1","Node2","Cor")
head(Cors2)

#bind pvalue and correlation tables

Edge_table <- cbind(Pvals2, Cors2$Cor, deparse.level=2)
head(Edge_table)

#remove NA lines

Edge_table_final <- Edge_table[!is.na(Edge_table$Pvalue),]

colnames(Edge_table_final) <- c("Node1","Node2","Pvalue","Cor")

#select correlations with another value that 0
Edge_table_final <- subset(Edge_table_final,abs(Edge_table_final$Cor)!=0)
##head(Edge_table_final)

# add a new column and give it value according to correlation values
Edge_table_final$CorrType <- ifelse(Edge_table_final$Cor >=0, "Positive", "Negative")
head(Edge_table_final)
dim(Edge_table_final)

#remove pvalues lower than 0.05
Edge_table_final <- Edge_table_final[which(Edge_table_final$Pvalue<=0.05),]

# number of positive and negative interactions
table(Edge_table_final$CorrType)

#save the file

write.table(Edge_table_final, file="network_output_fastspar/wl_nw_fp10_output/Pvalue_correlation_no_filter.txt", sep="\t", row.names=FALSE , quote = F)

