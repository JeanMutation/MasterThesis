setwd("D:/master_project_networks/workflows/otu_table")

library(gridExtra)
library(vegan)
library(dplyr)
library(ggplot2)
library(agricolae)
library(tibble)
library(tidyr)
library(reshape2)
library(VennDiagram)

whole_leaf <- read.table('raw_otu_and_tax/Approach6/_all_otu_table_normalized_relAB_wl_ap6.tsv', sep="\t", h=T) 
endophytes <- read.table('raw_otu_and_tax/Approach6/_all_otu_table_normalized_relAB_endo_ap6.tsv', sep="\t", h=T) 
otu_table <- merge(whole_leaf, endophytes, by = "Feature.ID", all = TRUE)

make_otu_readable <- function(otu_table_in){
  rownames(otu_table_in) <- otu_table_in$Feature.ID
  otu_table_in <- otu_table_in[, -1]
  otu_table_in <- t(otu_table_in)
  otu_table_in[is.na(otu_table_in)] <- 0
  otu_table_in <- as.data.frame(otu_table_in)
  return(otu_table_in)
}

whole_leaf <- make_otu_readable(whole_leaf)
endophytes <- make_otu_readable(endophytes)
otu_table <- make_otu_readable(otu_table)

otu_rows <- row.names(otu_table)

metadata <- read.table("metadata_collapsed.csv", sep=";", h=T) %>% 
  rename(Description = Identifier) %>% 
  mutate(Country = as.factor(Country)) %>% 
  mutate(Compartment = tolower(Compartment))%>% 
  mutate(Site = as.factor(Site)) %>%
  filter(Description %in% otu_rows)

metadata_wl <- metadata %>% filter(Compartment == "wl")
metadata_endo <- metadata %>% filter(Compartment == "endo")

## Generate the diversity indices per sample #############################################################

calculate_diversity <- function(otu_table_in, metadata_in) {
  species_richness <- rowSums(otu_table_in > 0)
  shannon_index <- diversity(otu_table_in, index = "shannon")
  simpson_index <- diversity(otu_table_in, index = "simpson")
  
  diversity_df <- data.frame(
    Description = rownames(otu_table_in),
    Species_richness = species_richness,
    Shannon_index = shannon_index,
    Simpson_index = simpson_index
  )
  
  merged_df <- merge(metadata_in, diversity_df, by = "Description", all.x = TRUE)
  
  return(merged_df)
}

metadata <- calculate_diversity(otu_table, metadata)
metadata_wl <- calculate_diversity(whole_leaf, metadata_wl)
metadata_endo <- calculate_diversity(endophytes, metadata_endo)

######################################################################################################################################
## plot the alpha diversity ##########################################################################################################

wilcox_test_statistic <- function(metadata_in, div_index){
  metadata_in$div_index <- metadata_in[[div_index]]
  pairwise_test <- pairwise.wilcox.test(metadata_in$div_index, metadata_wl$Country, p.adjust.method = "hochberg", exact = FALSE)
  pairwise_test$p.value
}

kruskal_test_statistic <- function(metadata_in, div_index){
  metadata_in$div_index <- metadata_in[[div_index]]
  kruskal_test <- kruskal.test(div_index ~ Country, data = metadata_in)
  kruskal_test$p.value
}

plot_alpha_div_country <- function(metadata_in, div_index, title_in, limit_y, text_y_shift, custom_letters, kruskal_p_value) {
  
  colors <- c('DE' = '#1f77b4','FR' = '#9467bd','SE' = '#bcbd22', 'ES' = '#17becf')
  
  metadata_in$div_index <- metadata_in[[div_index]]
  
  sig.letters <- custom_letters
  
  plot <- ggplot(metadata_in, aes(x = Country, y = div_index, fill = Country)) +
    geom_boxplot(outlier.size = 1) +
    labs(title = title_in, 
         subtitle = paste("Kruskal-Wallis p-value:", formatC(kruskal_p_value, digits = 2)),
         x = '', y = title_in) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = limit_y) +
    geom_text(data = sig.letters, aes(x = Country, y = text_y_shift, label = groups), vjust = -0.5) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.key.size = unit(6, 'mm'),
          legend.text = element_text(size = 8),
          legend.margin = margin(-5, -5, -5, -5),
          legend.box.margin = margin(3, 3, 3, 3),
          plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          axis.title = element_text(size = 8))
  
  return(plot)
}

#####

wl_species_richness_stat <- wilcox_test_statistic(metadata_wl, 'Species_richness')
wl_species_richness_stat

letters <- data.frame(
  Country = c("DE", "FR", "SE", "ES"),
  groups = c("a", "a", "a", "a")
)

kruskal_p_value <- kruskal_test_statistic(metadata_wl, 'Species_richness')
wl_species_richness_plot <- plot_alpha_div_country(metadata_wl, 'Species_richness', 'Species Richness', 
                                                   c(0,360), 330, letters, kruskal_p_value)
wl_species_richness_plot

##

wl_shannon_index_stat <- wilcox_test_statistic(metadata_wl, 'Shannon_index')
wl_shannon_index_stat

letters <- data.frame(
  Country = c("DE", "FR", "SE", "ES"),
  groups = c("a", "b", "ab", "b")
)

kruskal_p_value <- kruskal_test_statistic(metadata_wl, 'Shannon_index')
wl_shannon_index_plot <- plot_alpha_div_country(metadata_wl, 'Shannon_index', 'Shannon index', 
                                                c(1.5,4.7), 4.4, letters, kruskal_p_value)
wl_shannon_index_plot

##
wl_simpson_index_stat <- wilcox_test_statistic(metadata_wl, 'Simpson_index')
wl_simpson_index_stat

letters <- data.frame(
  Country = c("DE", "FR", "SE", "ES"),
  groups = c("a", "b", "ab", "b")
)

kruskal_p_value <- kruskal_test_statistic(metadata_wl, 'Simpson_index')
wl_simpson_index_plot <- plot_alpha_div_country(metadata_wl, 'Simpson_index', 'Simpson index', 
                                                c(0.749, 1.03), 1, letters, kruskal_p_value)
wl_simpson_index_plot

#####
endo_species_richness_stat <- wilcox_test_statistic(metadata_endo, 'Species_richness')
endo_species_richness_stat

letters <- data.frame(
  Country = c("DE", "FR", "SE", "ES"),
  groups = c("a", "a", "b", "a")
)

kruskal_p_value <- kruskal_test_statistic(metadata_endo, 'Species_richness')
endo_species_richness_plot <- plot_alpha_div_country(metadata_endo, 'Species_richness', 'Species Richness', 
                                                     c(0,280), 250, letters, kruskal_p_value)
endo_species_richness_plot

##
endo_shannon_index_stat <- wilcox_test_statistic(metadata_endo, 'Shannon_index')
endo_shannon_index_stat

letters <- data.frame(
  Country = c("DE", "FR", "SE", "ES"),
  groups = c("a", "a", "a", "a")
)

kruskal_p_value <- kruskal_test_statistic(metadata_endo, 'Shannon_index')
endo_shannon_index_plot <- plot_alpha_div_country(metadata_endo, 'Shannon_index', 'Shannon index', 
                                                  c(1.5,4.4), 4, letters, kruskal_p_value)
endo_shannon_index_plot

##
endo_simpson_index_stat <- wilcox_test_statistic(metadata_endo, 'Simpson_index')
endo_simpson_index_stat

letters <- data.frame(
  Country = c("DE", "FR", "SE", "ES"),
  groups = c("a", "a", "a", "a")
)

kruskal_p_value <- kruskal_test_statistic(metadata_endo, 'Simpson_index')
endo_simpson_index_plot <- plot_alpha_div_country(metadata_endo, 'Simpson_index', 'Simpson index', 
                                                  c(0.749, 0.989), 0.96, letters, kruskal_p_value)
endo_simpson_index_plot

# SAVE THE PLOTS ###############################################
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <-g_legend(endo_simpson_index_plot)

title1 <- textGrob("Alpha Diversity in the whole leaf compartment", gp = gpar(fontsize = 12, fontface = "bold"))
title2 <- textGrob("Alpha Diversity in the endophytic compartment", gp = gpar(fontsize = 12, fontface = "bold"))


pdf(file = 'D:/master_project_networks/workflows/otu_table/div_plots_final/alpha-div-wl.pdf',
    width = 6.7,
    height = 2.3)

grid.arrange(title1,
             arrangeGrob(wl_species_richness_plot + theme(legend.position="none"), 
                         wl_shannon_index_plot + theme(legend.position="none"), 
                         wl_simpson_index_plot + theme(legend.position="none"), 
                         ncol = 3, nrow = 1),
             mylegend, nrow=3,heights=c(1, 10, 1))


dev.off()

pdf(file = 'D:/master_project_networks/workflows/otu_table/div_plots_final/alpha-div-endo.pdf',
    width = 6.7,
    height = 2.5)

grid.arrange(title2,
             arrangeGrob(endo_species_richness_plot + theme(legend.position="none"), 
                         endo_shannon_index_plot + theme(legend.position="none"), 
                         endo_simpson_index_plot + theme(legend.position="none"), 
                         ncol = 3, nrow = 1),
             mylegend, nrow=3,heights=c(1, 10, 1))


dev.off()


#######################################################################################################################
# Random Subsampling ##################################################################################################
######################################################################################################################
n_iterations <- 100

p_values_df <- data.frame(
  Iteration = 1:n_iterations,
  Species_richness_p = NA,
  Shannon_index_p = NA,
  Simpson_index_p = NA
)

perform_subsampling <- function(data, min_count) {
  subsampled_data <- data %>%
    group_by(Country) %>%
    sample_n(min_count) %>%
    ungroup()
  return(subsampled_data)
}

perform_kruskal_test <- function(subsampled_data, div_index) {
  kruskal_test <- kruskal.test(subsampled_data[[div_index]] ~ subsampled_data$Country)
  return(kruskal_test$p.value)
}

for (i in 1:n_iterations) {
  set.seed(400*i)
  country_counts <- metadata_endo %>%
    count(Country)
  
  min_count <- min(country_counts$n)
  
  subsampled_data <- perform_subsampling(metadata_endo, min_count)
  
  species_richness_p <- perform_kruskal_test(subsampled_data, 'Species_richness')
  shannon_index_p <- perform_kruskal_test(subsampled_data, 'Shannon_index')
  simpson_index_p <- perform_kruskal_test(subsampled_data, 'Simpson_index')
  
  p_values_df$Species_richness_p[i] <- species_richness_p
  p_values_df$Shannon_index_p[i] <- shannon_index_p
  p_values_df$Simpson_index_p[i] <- simpson_index_p
}

print(p_values_df)


count_below_0_05 <- colSums(p_values_df[-1] < 0.05)
percentage_below_0_05 <- (count_below_0_05 / n_iterations) * 100

results_df <- data.frame(
  Column = names(count_below_0_05),
  Count_below_0_05 = count_below_0_05,
  Percentage_below_0_05 = percentage_below_0_05
)

print(results_df)

