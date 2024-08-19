setwd("D:/master_project_networks/workflows/otu_table")

library("FactoMineR")
library('corrr')
library(ggcorrplot)
library(vegan)
library(ade4)
library(ggplot2)
library(ellipse)
library(tidyverse)
library(tibble)
library(dplyr)
library(car)
library(gridExtra)

## Read in the dataframes and clean them for further analysis ####################################################################

whole_leaf <- read.table('raw_otu_and_tax/Approach6/_all_otu_table_normalized_relAB_wl_ap6.tsv', sep="\t", h=T)
endophytes <- read.table('raw_otu_and_tax/Approach6/_all_otu_table_normalized_relAB_endo_ap6.tsv', sep="\t", h=T)
otu_table <- merge(whole_leaf, endophytes, by = "Feature.ID", all = TRUE)


clean_tables_for_analysis <- function(table_input){
  rownames(table_input) <- table_input$Feature.ID
  table_input <- table_input[, -1]
  table_input <- t(table_input)
  table_input[is.na(table_input)] <- 0
  return(table_input)
}

whole_leaf <- clean_tables_for_analysis(whole_leaf)
endophytes <- clean_tables_for_analysis(endophytes)
otu_table <- clean_tables_for_analysis(otu_table)

otu_rows <- row.names(otu_table)

metadata <- read.table("metadata_collapsed.csv", sep=";", h=T) %>% 
  rename(Description = Identifier) %>% 
  mutate(Country = as.factor(Country)) %>% 
  mutate(Compartment = tolower(Compartment))%>%
  mutate(Compartment = as.factor(Compartment)) %>% 
  filter(Description %in% otu_rows) 
  
metadata_wl <- metadata %>% filter(Compartment == "wl")
metadata_endo <- metadata %>% filter(Compartment == "endo")

###########################################################################################################################
## PCOA ##################################################################################################################
###########################################################################################################################

# For reproducability
set.seed(19980607)

calculate_bray_curtis_distance_matrix <- function(otu_table_in){
  # construction of distance matrices
  distmat.pca_bray_otu_table_con <- vegdist(otu_table_in, method="bray", binary=FALSE)
  pca_bray_otu_table_con <- dudi.pco(vegdist(otu_table_in, method="bray", binary=FALSE), scannf=F)
  
  # Extract eigenvalues from the dudi.pco result to Calculate the proportion of variance explained
  eigenvalues_con <- pca_bray_otu_table_con$eig
  variance_proportion_con <- eigenvalues_con / sum(eigenvalues_con)
  
  plot(1:length(variance_proportion_con), cumsum(variance_proportion_con), type = "b", 
       main = "Cumulative Proportion of Variance Explained - Bray-Curtis PCoA",
       xlab = "Number of Principal Coordinates", ylab = "Cumulative Proportion")
  
  for (i in 1:4) {
    cat("Axis", i, "explains", round((variance_proportion_con[[i]]*100), digits = 2), "% of the variance.\n")
  }
  
  distmat.pca_bray_otu_table_con <- as.matrix(distmat.pca_bray_otu_table_con)
  
  return(list(dudi_obj = pca_bray_otu_table_con, 
              distmat = distmat.pca_bray_otu_table_con, 
              axis1 = round((variance_proportion_con[[1]]*100), digits = 2),
              axis2 = round((variance_proportion_con[[2]]*100), digits = 2)))
}

plot_bray_curtis_distance_matrix_country <- function(dist_mat_in, metadata_in, Axis1, Axis2, title_in, pval_in){
  metadata_in$Compartment <- factor(metadata_in$Compartment, levels = c("wl", "endo"))
  
  plot_data <- as.data.frame(dist_mat_in$li) %>% 
    rownames_to_column(var = "Description") %>%
    inner_join(metadata_in, by = "Description") %>% 
    mutate(row.names = Description) %>%
    select(-Description) %>%
    column_to_rownames(var = "row.names") %>% 
    rename(x = A1, y = A2, group = Country, shapes = Compartment)
  
  colors <- c('DE' = '#1f77b4','FR' = '#9467bd','SE' = '#bcbd22', 'ES' = '#17becf')
  # colors <- c('DE' = 'deeppink','FR' = 'darkorchid','SE' = 'dodgerblue', 'ES' = 'blue')
  
  plot_pcoa_bray_con <- ggplot(plot_data, aes(x = x, y = y, color = group)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
    geom_point(size = 1, aes(shape = shapes)) +  
    stat_ellipse(aes(fill = group), geom = "polygon", alpha = 0.1, show.legend = FALSE, size = 0.3) +  # Ellipses
    theme_minimal() +
    labs(title = paste("PCoA"), color = "Countries", shape = 'Compartment', 
         subtitle = paste("Adonis2 p-value", round(pval_in, 6))) +
    xlab(paste0('PCoA Axis 1 (', Axis1, '%)')) +
    ylab(paste0('PCoA Axis 2 (', Axis2, '%)')) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 11), 
          plot.subtitle = element_text(hjust = 0.5, size = 8), legend.box="vertical", 
          legend.margin=margin(-5,-5,-5,-5), legend.box.margin=margin(3,3,3,3),
          legend.spacing.x = unit(0.1, "cm"),
          axis.title = element_text(size = 8)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_shape_discrete(drop = FALSE)
  
  return(plot_pcoa_bray_con)
}

plot_bray_curtis_distance_matrix_country_kingdom <- function(dist_mat_in, metadata_in, Axis1, Axis2, title_in, pval_in, kingdom_in){
  metadata_in$Compartment <- factor(metadata_in$Compartment, levels = c("wl", "endo"))
  
  plot_data <- as.data.frame(dist_mat_in$li) %>% 
    rownames_to_column(var = "Description") %>%
    inner_join(metadata_in, by = "Description") %>% 
    mutate(row.names = Description) %>%
    select(-Description) %>%
    column_to_rownames(var = "row.names") %>% 
    rename(x = A1, y = A2, group = Country, shapes = Compartment)
  
  # colors <- c('DE' = 'deeppink','FR' = 'darkorchid','SE' = 'dodgerblue', 'ES' = 'blue')
  colors <- c('DE' = '#1f77b4','FR' = '#9467bd','SE' = '#bcbd22', 'ES' = '#17becf')
  
  plot_pcoa_bray_con <- ggplot(plot_data, aes(x = x, y = y, color = group)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
    geom_point(size = 1, aes(shape = shapes)) +  
    stat_ellipse(aes(fill = group), geom = "polygon", alpha = 0.1, show.legend = FALSE, size = 0.3) +  # Ellipses
    theme_minimal() +
    labs(title = paste("PCoA: ", title_in, ' - ', kingdom_in), color = "Countries", shape = 'Compartment',
         subtitle = paste("Adonis2 p-value", round(pval_in, 5))) +
    xlab(paste0('PCoA Axis 1 (', Axis1, '%)')) +
    ylab(paste0('PCoA Axis 2 (', Axis2, '%)')) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 11), 
          plot.subtitle = element_text(hjust = 0.5), legend.box="vertical", 
          legend.margin=margin(-5,-5,-5,-5), legend.box.margin=margin(3,3,3,3),
          legend.spacing.x = unit(0.1, "cm"),
          axis.title = element_text(size = 8)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_shape_discrete(drop = FALSE)
  
  return(plot_pcoa_bray_con)
}

statistic_bray_curtis_distance_matrix_country <- function(dist_mat_in, metadata_in){
  dist_df_bray <- as.data.frame(as.matrix(dist_mat_in))
  dist_df_bray$Description <- rownames(dist_df_bray)
  
  meta_dist_df_bray <- inner_join(metadata_in, dist_df_bray, by = 'Description')
  
  permutations <- 10000
  permanova_test_bray <- adonis(dist_mat_in ~ Country, data = meta_dist_df_bray,  
                                permutations = permutations, method = "bray")
  cat(permanova_test_bray$aov.tab$`Pr(>F)`)
  return(list(all = permanova_test_bray, pval = permanova_test_bray$aov.tab$`Pr(>F)`[[1]]))
}

analyze_pcoa <- function(otu_table_in, metadata_in, title_in){
  otu_bray_curtis_distance <- calculate_bray_curtis_distance_matrix(otu_table_in)
  otu_bray_curtis_distance_dudi <- otu_bray_curtis_distance$dudi_obj
  otu_bray_curtis_distance_dist <- otu_bray_curtis_distance$distmat
  otu_bray_curtis_distance_statistic <- statistic_bray_curtis_distance_matrix_country(otu_bray_curtis_distance_dist, 
                                                                                      metadata_in)
  otu_bray_curtis_distance_plot <- plot_bray_curtis_distance_matrix_country(otu_bray_curtis_distance_dudi, 
                                                                            metadata_in, 
                                                                            otu_bray_curtis_distance$axis1, 
                                                                            otu_bray_curtis_distance$axis2,
                                                                            title_in, otu_bray_curtis_distance_statistic$pval
                                                                            )
  return(list(stat = otu_bray_curtis_distance_statistic, 
              plot = otu_bray_curtis_distance_plot))
}

analyze_pcoa_by_kingdom <- function(otu_table_in, metadata_in, title_in, kingdom_in){
  kingdom_cols <- which(grepl(kingdom_in, colnames(otu_table_in)))
  otu_table_in <- otu_table_in[, kingdom_cols, drop = FALSE]
  
  otu_bray_curtis_distance <- calculate_bray_curtis_distance_matrix(otu_table_in)
  otu_bray_curtis_distance_dudi <- otu_bray_curtis_distance$dudi_obj
  otu_bray_curtis_distance_dist <- otu_bray_curtis_distance$distmat
  otu_bray_curtis_distance_statistic <- statistic_bray_curtis_distance_matrix_country(otu_bray_curtis_distance_dist, 
                                                                                      metadata_in)
  
  otu_bray_curtis_distance_plot <- plot_bray_curtis_distance_matrix_country_kingdom(otu_bray_curtis_distance_dudi, 
                                                                            metadata_in, 
                                                                            otu_bray_curtis_distance$axis1, 
                                                                            otu_bray_curtis_distance$axis2,
                                                                            title_in, 
                                                                            otu_bray_curtis_distance_statistic$pval,
                                                                            kingdom_in)
  return(list(stat = otu_bray_curtis_distance_statistic, 
              plot = otu_bray_curtis_distance_plot))
}

##########################################
otu_all_bray_curtis_distance <- analyze_pcoa(otu_table, metadata, 'All Samples')
otu_all_bray_curtis_distance$plot

otu_wl_bray_curtis_distance <- analyze_pcoa(whole_leaf, metadata_wl, 'Whole leaf')
otu_wl_bray_curtis_distance$plot

otu_endo_bray_curtis_distance <- analyze_pcoa(endophytes, metadata_endo, 'Endophytes')
otu_endo_bray_curtis_distance$plot


pdf(file = 'D:/master_project_networks/workflows/otu_table/div_plots/beta-div-plot-otu-wl.pdf',
    width = 2.5,
    height = 2.5)

otu_wl_bray_curtis_distance$plot + theme(legend.position="none")

dev.off()