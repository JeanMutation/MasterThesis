setwd("D:/master_project_networks/workflows/otu_table")

library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")     # Needed for converting column to row names
library(ggpubr)
library(ggtree)
library(gtable)
library("gridExtra")
library(tidyr)
library(cowplot)

## reading the data and splitting wl and endo ##################################################

whole_leaf <- read.table('raw_otu_and_tax/Approach6/_all_otu_table_normalized_relAB_wl_ap6.tsv', sep="\t", h=T)
whole_leaf <- whole_leaf %>% 
  tibble::column_to_rownames("Feature.ID")

endophytes <- read.table('raw_otu_and_tax/Approach6/_all_otu_table_normalized_relAB_endo_ap6.tsv', sep="\t", h=T)
endophytes <- endophytes %>% 
  tibble::column_to_rownames("Feature.ID")

taxonomy <- read.table("raw_otu_and_tax/Approach6/taxonomy.tsv", sep="\t", h=T)
rownames(taxonomy) <- taxonomy$Feature.ID
taxonomy$sample_kingdom <- sub("^[^_]+_(.*)", "\\1", taxonomy$Feature.ID)
taxonomy <- subset(taxonomy, select = -c(Feature.ID))
taxonomy$Genus <- trimws(taxonomy$Genus, "left")
taxonomy <- as.matrix(taxonomy)


metadata <- read.table("metadata_collapsed.csv", sep=";", h=T)
rownames(metadata) <- metadata$Identifier
metadata <- subset(metadata, select = -c(Identifier))

metadata_wl <- metadata[grep("wl", metadata$Compartment), ]
metadata_endo <- metadata[grep("endo", metadata$Compartment), ]

## Create phyloseq objects #############################################

OTU_wl = otu_table(whole_leaf, taxa_are_rows = TRUE)
OTU_endo = otu_table(endophytes, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
samples_wl = sample_data(metadata_wl)
samples_endo = sample_data(metadata_endo)

obj_wl <- phyloseq(OTU_wl, TAX, samples_wl)
obj_endo <- phyloseq(OTU_endo, TAX, samples_endo)

sample_names(obj_wl)
rank_names(obj_wl)
sample_variables(obj_wl)

sample_names(obj_endo)
rank_names(obj_endo)
sample_variables(obj_endo)

## Create the barplots #################################################################################################
# merge by Country
obj_wl_merged_site <- merge_samples(obj_wl, 'Country')
obj_endo_merged_site <- merge_samples(obj_endo, 'Country')

# rewrite the metadata
columns_to_delete <- c("InfectionStatus", 'Year', 'YearTotal', 'ReplicateNr', 'Run')
metadata_wl <- sample_data(obj_wl_merged_site)
metadata_wl[, columns_to_delete] <- NULL
metadata_wl$Compartment <- c('whole-leaf', 'whole-leaf', 'whole-leaf', 'whole-leaf')
metadata_wl$Country <- c('Germany', 'Spain', 'France', 'Sweden')
obj_wl_merged_site <- phyloseq(otu_table(obj_wl_merged_site), tax_table(obj_wl_merged_site), sample_data(metadata_wl))

metadata_endo <- sample_data(obj_endo_merged_site)
metadata_endo[, columns_to_delete] <- NULL
metadata_endo$Compartment <- c('endophyte', 'endophyte', 'endophyte', 'endophyte')
metadata_endo$Country <- c('Germany', 'Spain', 'France', 'Sweden')
obj_endo_merged_site <- phyloseq(otu_table(obj_endo_merged_site), tax_table(obj_endo_merged_site), sample_data(metadata_endo))

# subset the data into the different kingdoms
obj_wl_bac <- subset_taxa(obj_wl_merged_site, sample_kingdom == 'bacteria')
obj_wl_fun <- subset_taxa(obj_wl_merged_site, sample_kingdom == 'fungi')
obj_wl_oom <- subset_taxa(obj_wl_merged_site, sample_kingdom == 'oomycete')
obj_wl_euk <- subset_taxa(obj_wl_merged_site, sample_kingdom == 'eukaryote')

obj_endo_bac <- subset_taxa(obj_endo_merged_site, sample_kingdom == 'bacteria')
obj_endo_fun <- subset_taxa(obj_endo_merged_site, sample_kingdom == 'fungi')
obj_endo_oom <- subset_taxa(obj_endo_merged_site, sample_kingdom == 'oomycete')
obj_endo_euk <- subset_taxa(obj_endo_merged_site, sample_kingdom == 'eukaryote')

# Normalize the data again
standf <- function(x) (x / sum(x))


norm_obj_wl_bac <- transform_sample_counts(obj_wl_bac, standf)
norm_obj_wl_fun <- transform_sample_counts(obj_wl_fun, standf)
norm_obj_wl_oom <- transform_sample_counts(obj_wl_oom, standf)
norm_obj_wl_euk <- transform_sample_counts(obj_wl_euk, standf)

norm_obj_endo_bac <- transform_sample_counts(obj_endo_bac, standf)
norm_obj_endo_fun <- transform_sample_counts(obj_endo_fun, standf)
norm_obj_endo_oom <- transform_sample_counts(obj_endo_oom, standf)
norm_obj_endo_euk <- transform_sample_counts(obj_endo_euk, standf)

plot_bar(norm_obj_wl_oom, 'Country','Abundance', 'Genus')
plot_bar(norm_obj_endo_euk, 'Country','Abundance', 'Class')

# trim to the 10 most abundant classes

top_10_by_class <- function(norm_obj) {
  TaxaSums <- data.frame(Class = tax_table(norm_obj)[,"Class"],
                         taxa_sums = taxa_sums(norm_obj)) %>%
    group_by(Class) %>%
    summarize(taxa_sums = sum(taxa_sums)) %>% 
    arrange(desc(taxa_sums))
  
  top10_level <- head(TaxaSums$Class, 10)
  
  taxa_sum_top_10 <- data.frame(Class = tax_table(norm_obj)[,"Class"],
                                taxa_sums = taxa_sums(norm_obj)) %>%
    arrange(desc(taxa_sums)) 
  
  taxa_sum_top_10 <- subset(taxa_sum_top_10, Class %in% top10_level)
  
  top10 <- rownames(taxa_sum_top_10)
  
  y <- prune_taxa(top10, norm_obj)
  
  df1 <- data.frame(ID = c(taxa_names(y), "Other"),
                    Phylum = c(tax_table(y)[,"Class"], "Other"))
  
  df2 <- t(cbind(otu_table(y), data.frame(low_abundant = 1 - sample_sums(y))))
  
  df <- cbind(df1, df2) %>%
    pivot_longer(-c(ID, Phylum), names_to = "SampleType", values_to = "Abundance") %>%
    as.data.frame
  
  df$Phylum = as.factor(df$Phylum)
  
  df <- df %>% 
    group_by(Phylum, SampleType) %>%
    summarize(Abundance = sum(Abundance)) %>%
    ungroup() %>%
    group_by(Phylum) %>% 
    filter(Abundance != 0) %>% 
    ungroup() %>%
    mutate(across(Phylum, as.character)) %>%
    mutate_all(~ifelse(. == "", 'unclassified', .)) %>%
    mutate(across(Phylum, as.factor))
  
  return(df)
}

top_10_by_genus <- function(norm_obj) {
  TaxaSums <- data.frame(Genus = tax_table(norm_obj)[,"Genus"],
                         taxa_sums = taxa_sums(norm_obj)) %>%
    group_by(Genus) %>%
    summarize(taxa_sums = sum(taxa_sums)) %>% 
    arrange(desc(taxa_sums))
  
  top10_level <- head(TaxaSums$Genus, 10)
  
  taxa_sum_top_10 <- data.frame(Genus = tax_table(norm_obj)[,"Genus"],
                                taxa_sums = taxa_sums(norm_obj)) %>%
    arrange(desc(taxa_sums)) 
  
  taxa_sum_top_10 <- subset(taxa_sum_top_10, Genus %in% top10_level)
  
  top10 <- rownames(taxa_sum_top_10)
  
  y <- prune_taxa(top10, norm_obj)
  
  df1 <- data.frame(ID = c(taxa_names(y), "Other"),
                    Phylum = c(tax_table(y)[,"Genus"], "Other"))
  
  df2 <- t(cbind(otu_table(y), data.frame(Other = 1 - sample_sums(y))))
  
  df <- cbind(df1, df2) %>%
    pivot_longer(-c(ID, Phylum), names_to = "SampleType", values_to = "Abundance") %>%
    as.data.frame
  
  df$Phylum = as.factor(df$Phylum)
  
  df <- df %>% 
    group_by(Phylum, SampleType) %>%
    summarize(Abundance = sum(Abundance)) %>%
    ungroup() %>%
    group_by(Phylum) %>% 
    filter(Abundance != 0) %>% 
    ungroup() %>%
    mutate(across(Phylum, as.character)) %>%
    mutate_all(~ifelse(. == "", 'unclassified', .)) %>%
    mutate(across(Phylum, as.factor))
  
  return(df)
}

top_10_wl_bac <- top_10_by_class(norm_obj_wl_bac)
top_10_wl_fun <- top_10_by_class(norm_obj_wl_fun)
top_10_wl_oom <- top_10_by_genus(norm_obj_wl_oom)
top_10_wl_euk <- top_10_by_class(norm_obj_wl_euk)

top_10_endo_bac <- top_10_by_class(norm_obj_endo_bac)
top_10_endo_fun <- top_10_by_class(norm_obj_endo_fun)
top_10_endo_oom <- top_10_by_genus(norm_obj_endo_oom)
top_10_endo_euk <- top_10_by_class(norm_obj_endo_euk)

average_abundance <- top_10_endo_euk %>%
  group_by(Phylum) %>%
  summarize(Average_Abundance = mean(Abundance)) %>%
  arrange(desc(Average_Abundance))

print(average_abundance)

# plot the data
plot_barplot_class <- function(top_10_obj_wl, top_10_obj_endo, Kingdom){
  top_10_obj_wl$Dataset <- "whole-leaf"
  top_10_obj_endo$Dataset <- "endophytes"
  combined_data <- bind_rows(top_10_obj_wl, top_10_obj_endo)
  
  ggplot(combined_data, aes(SampleType, Abundance, fill = Phylum, group = Dataset)) +
    geom_col(color = NA) + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme(axis.title = element_text(size = 8),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          legend.key.size = unit(0.15, 'cm')) +
    labs(x = " ", y = "Average relative abundance", fill = paste0(Kingdom, ' - Class')) +
    font("xlab", size = 8) +
    facet_wrap(~Dataset, scales = "free_y")
}

plot_barplot_genus <- function(top_10_obj_wl, top_10_obj_endo, Kingdom){
  top_10_obj_wl$Dataset <- "whole-leaf"
  top_10_obj_endo$Dataset <- "endophytes"
  combined_data <- bind_rows(top_10_obj_wl, top_10_obj_endo)
  
  ggplot(combined_data, aes(SampleType, Abundance, fill = Phylum, group = Dataset)) +
    geom_col(color = NA) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme(axis.title = element_text(size = 8),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          legend.key.size = unit(0.15, 'cm')) +
    labs(x = " ", y = "Average relative abundance", fill = paste0(Kingdom, ' - Genus')) +
    font("xlab", size = 8) +
    facet_wrap(~Dataset, scales = "free_y")
}

barplot_bac <- plot_barplot_class(top_10_wl_bac, top_10_endo_bac, 'Bacteria')
barplot_bac
barplot_fun <- plot_barplot_class(top_10_wl_fun, top_10_endo_fun, 'Fungi')
barplot_fun
barplot_oom <- plot_barplot_genus(top_10_wl_oom, top_10_endo_oom, 'Oomycetes')
barplot_oom
barplot_euk <- plot_barplot_class(top_10_wl_euk, top_10_endo_euk, 'Eukaryotes')
barplot_euk

# SAVE THE DATA ###########################################################################################

pdf(file = 'D:/master_project_networks/workflows/otu_table/div_plots_final/abundance_top10_otus.pdf',
    width = 6.7,
    height = 9)

plot_grid(barplot_bac,
          barplot_fun,
          barplot_oom,
          barplot_euk,
          align = 'v',
          nrow = 4, ncol = 1)

dev.off()

pdf(file = 'D:/master_project_networks/workflows/otu_table/div_plots_final/abundance_top10_otus_wo_legend.pdf',
    width = 6.7,
    height = 9)

plot_grid(barplot_bac + theme(legend.position="none"),
          barplot_fun + theme(legend.position="none"),
          barplot_oom + theme(legend.position="none"),
          barplot_euk + theme(legend.position="none"),
          align = 'v',
          nrow = 4, ncol = 1)

dev.off()

#########################################################################################################
