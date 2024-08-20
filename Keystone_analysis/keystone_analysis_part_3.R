

library(ggplot2)
library(dplyr)
library(ggrepel)
library(reshape2) 
library(viridis)
library("ggpubr")
library(gridExtra)
library(vegan)
library(ade4)

setwd("D:/master_keystone_project/workflows/keystone_analysis_run_5_data")

real_comparison <- read.csv('Results/Comparison_full_and_dropout_nws.tsv', sep = '\t')

taxonomy <- read.table("taxonomy.tsv", sep="\t", h=T) %>%
 select(c(Feature.ID, Order, Family, Genus))

core_wl <- read.table("../Important_otus_wl.csv", sep=",", h=T) 
core_wl_list <- as.list(core_wl$X)

core_endo <- read.table("../Important_otus_endo.csv", sep=",", h=T) 
core_endo_list <- as.list(core_endo$X)

# taxonomy <- read.table("taxonomy.tsv", sep="\t", h=T) %>% 
#     select(c(OTUid, Order, Family, Genus))

clean_comparison_df <- function(comparison_df){
  comparison_df <- comparison_df %>% 
    select(-c(X, List_Network)) %>% 
    mutate(Shared_Edges_percent = Shared_Edges / Theoretical_edge_num,
           Shared_Edges_with_attributes_percent = Shared_Edges_with_attributes / Theoretical_edge_num,
           Kingdom = sub("^[^_]+_", "", comparison_df$Dict_Network))
  
  return(comparison_df)
}

real_comparison <- clean_comparison_df(real_comparison)

dropout_degree_density <- read.csv('Results/Dropout_nw_properties_1_degree_density.csv')
real_comparison <- inner_join(real_comparison, dropout_degree_density, by = c("Dict_Network" = "OTU_Name"))

dropout_modularity <- read.csv('Results/Dropout_nw_properties_2_modularity.tsv', sep='\t')
real_comparison <- inner_join(real_comparison, dropout_modularity, by = c("Dict_Network" = "OTU_Name"))

full_degree_density <- read.csv('Results/Full_nw_properties_1_degree_density.tsv', sep='\t')
real_comparison <- inner_join(real_comparison, full_degree_density, by = c("Dict_Network" = "X"))

full_centrality <- read.csv('Results/Full_nw_properties_2_centrality.tsv', sep='\t')
real_comparison <- inner_join(real_comparison, full_centrality, by = c("Dict_Network" = "X"))

full_module_nr <- read.csv('Results/Full_nw_properties_3_module_numbers.tsv', sep='\t') %>% select(c(Node, Module))
real_comparison <- inner_join(real_comparison, full_module_nr, by = c("Dict_Network" = "Node"))

real_comparison$Keystone_candidate <- ifelse(
  real_comparison$Shared_Edges_with_attributes_percent < 0.70, 'yes','no')

real_comparison <- inner_join(real_comparison, taxonomy, by = c("Dict_Network" = 'Feature.ID'))

real_comparison$core_microbe <- ifelse(real_comparison$Dict_Network %in% core_endo_list, "yes", "no")

write.csv(real_comparison, 'Results/Real_Comparison.csv')

####

random_comparison <- read.csv('Results/Result_comparison_random_networks.tsv', sep='\t')

random_comparison <- random_comparison %>%
  group_by(Dict_Network) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  mutate(Shared_Edges_percent = Shared_Edges / Theoretical_edge_num,
         Shared_Edges_with_attributes_percent = Shared_Edges_with_attributes / Theoretical_edge_num)

random_comparison <- random_comparison %>%
  mutate(Kingdom = sub("^[^_]+_", "", random_comparison$Dict_Network))

## Plot a histogram
min_degree <- 0 - 0.55
max_degree <- max(real_comparison$Degree) + 0.55

plot_histogram_shared_edge_perc <- function(comparison_df_short) {
  ggplot(comparison_df_short, aes(x = Shared_Edges_with_attributes_percent)) +
    geom_histogram(binwidth = 0.01, fill = "slateblue", color = "black", size = 0.1) +
    labs(#title = "Common Edges between \nfull and dropout networks", 
         x = "Shared Edges with Attributes [%]", 
         y = "Number of dropout nodes") +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 10)) +
    xlim(-0.05, 1.05)
}

plot_histogram_Degree <- function(comparison_df_short) {
  ggplot(comparison_df_short, aes(x = Degree)) +
    geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
    labs(title = "Degree of the Nodes of\nthe full network", 
         x = "Degree", 
         y = "Frequency") +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 12)) +
    xlim(min_degree, max_degree)
}

plot_histogram_Degree <- function(comparison_df_short) {
  ggplot(comparison_df_short, aes(x = Degree, fill = Keystone_candidate)) +
    geom_bar(binwidth = 1, color = "black", , size = 0.1, position = "stack") +
    scale_fill_manual(values = c("skyblue", "slateblue")) +  # You can adjust the colors as needed
    labs(title = "Degree of the Nodes of\nthe full network - endophytes", 
         x = "Degree", 
         y = "Frequency",
         fill = "Keystone Candidate") +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 12),
          legend.position="bottom") +
    scale_x_continuous(breaks = seq(min(comparison_df_short$Degree, na.rm = TRUE), 
                                    max(comparison_df_short$Degree, na.rm = TRUE), 
                                    by = 1),
                       limits = c(min(comparison_df_short$Degree, na.rm = TRUE) - 0.55, 
                                  max(comparison_df_short$Degree, na.rm = TRUE) +0.55))
}


hist_shared_edges <- plot_histogram_shared_edge_perc(real_comparison)
hist_degree_wl <- plot_histogram_Degree(real_comparison)

hist_shared_edges
hist_degree_wl

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <-g_legend(hist_degree_wl)

pdf(file = '../final_plots/Hist_SharedEdges_real_endo_fp10p.pdf', width = 3.5, height = 2.5)
hist_shared_edges + labs(title = "endophytes - real networks")
dev.off()

pdf(file = '../final_plots/Hist_Degree_real_wl_with_legend.pdf', width = 3, height = 2)
hist_degree_wl + labs(title = "whole leaf") + theme(legend.position="none")
dev.off()

# pdf(file = 'Results/Hist_SharedEdges_Degree_real.pdf', width = 6.7, height = 3)
# grid.arrange(hist_shared_edges, hist_degree, ncol=2)
# dev.off()

hist_shared_edges_rand <- plot_histogram_shared_edge_perc(random_comparison)
hist_degree_rand <- plot_histogram_Degree(random_comparison)

hist_shared_edges_rand
hist_degree_rand

pdf(file = '../final_plots/Hist_SharedEdges_random_endo_fp10p.pdf', width = 2.5, height = 2.5)
hist_shared_edges_rand + labs(title = "endophytes - random networks")
dev.off()


pdf(file = 'Results/Hist_SharedEdges_Degree_random.pdf', width = 6.7, height = 3)
grid.arrange(hist_shared_edges_rand, hist_degree_rand, ncol=2)
dev.off()


## Hus analysis
n <- 5

top5_BC_wl <- real_comparison[real_comparison$Betweenness > 
                                      quantile(real_comparison$Betweenness, prob=1-n/100), ]$Dict_Network

top5_BC_wl_min <- min(real_comparison[real_comparison$Betweenness > 
                                              quantile(real_comparison$Betweenness, prob=1-n/100), ]$Betweenness)

top5_CC_wl <- real_comparison[real_comparison$Closeness > 
                                      quantile(real_comparison$Closeness, prob=1-n/100), ]$Dict_Network

top5_CC_wl_min <- min(real_comparison[real_comparison$Closeness > 
                                              quantile(real_comparison$Closeness, prob=1-n/100), ]$Closeness)


p5_Intersection_BC_CC_wl <- top5_BC_wl[which((top5_BC_wl %in% top5_CC_wl) == TRUE)]
real_comparison$hub_BC_CC <- ifelse(real_comparison$Dict_Network %in% p5_Intersection_BC_CC_wl,'yes','no')



wl_BC_CC_keystone <- ggplot(real_comparison, aes(x = Closeness, y= Betweenness , label = Family, shape = Kingdom)) +
  geom_point(aes(colour = Keystone_candidate) ,stat= "identity",position="identity",alpha= 1 ) +
  scale_colour_manual(values =  c('skyblue', 'slateblue')) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +  # Adjust shapes accordingly
  geom_vline(xintercept = top5_CC_wl_min, linetype="dashed", linewidth = 0.3) + 
  geom_hline(yintercept = top5_BC_wl_min, linetype= "dashed", linewidth = 0.3) + 
  labs(title = "Hubs - endophytic samples") +
  theme_minimal() + 
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing = unit(0.5, "cm")) +
  geom_text_repel(data = subset(real_comparison, Keystone_candidate == "yes"),
                  size = 3,
                  box.padding = unit(0.5, "lines"),    # Adjust box padding
                  point.padding = unit(0.3, "lines"),  # Adjust point padding
                  force = 100,                        # Increase force to spread labels
                  segment.size = 0.2,
                  segment.color = "darkgreen") +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) 

wl_BC_CC_keystone

wl_BC_CC_core <- ggplot(real_comparison, aes(x = Closeness, y= Betweenness , label = Family, shape = Kingdom)) +
  geom_point(aes(colour = core_microbe) ,stat= "identity",position="identity",alpha= 1 ) +
  scale_colour_manual(values =  c('skyblue', '#A02B93')) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +  # Adjust shapes accordingly
  labs(title = "Hubs - endophytic samples") +
  theme_minimal() + 
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing = unit(0.5, "cm"))

wl_BC_CC_core

wl_BC_CC_hub <- ggplot(real_comparison, aes(x = Closeness, y= Betweenness , label = Family, shape = Kingdom)) +
  geom_point(aes(colour = Kingdom) ,stat= "identity",position="identity",alpha= 1 ) +
  scale_colour_manual(values =  c('#8c564b', '#7f7f7f', '#ff7f0e', '#e377c2' )) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +  # Adjust shapes accordingly
  geom_vline(xintercept = top5_CC_wl_min, linetype="dashed", linewidth = 0.3) + 
  geom_hline(yintercept = top5_BC_wl_min, linetype= "dashed", linewidth = 0.3) + 
  labs(title = "Hubs - whole leaf samples") +
  theme_minimal() + 
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom") +
  geom_text_repel(data=subset(real_comparison, hub_BC_CC =="yes")
                  ,size = 3,
                  box.padding   = 1,
                  point.padding = 0.2,
                  force         = 50,
                  segment.size  = 0.2,
                  segment.color = "darkgreen") 

wl_BC_CC_hub

pdf(file = '../final_plots/keystone_endo_hubs_keystone_label.pdf', width = 6.7, height = 3.5)
wl_BC_CC_keystone
dev.off()

pdf(file = '../final_plots/keystone_wl_hubs_hubs_label.pdf', width = 6.7, height = 3)
wl_BC_CC_hub
dev.off()


pdf(file = '../final_plots/Hubs.pdf', width = 6.7, height = 4)

grid.arrange(wl_BC_CC_hub, wl_BC_CC_keystone + theme(legend.position="none"), ncol=2)

dev.off()


## Test correlations
correlation_matrix <- cor(real_comparison[c("Shared_Edges_with_attributes_percent", 
                                                  "Density", "Num_Subnetworks", "Modularity_Coefficient", 
                                                  "Degree", "Betweenness", "Closeness")])
melted_correlation <- melt(correlation_matrix)

# Plot heatmap
correlation_factors <- ggplot(melted_correlation, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis(name = "Correlation") +
  theme_minimal() +
  labs(title = "Pairwise Correlation Heatmap") +
  theme(axis.title.x = element_text(size = 5),
        axis.title.y = element_text(size = 5))

pdf(file = 'Results/Correlation_matrix.pdf', width = 6.7, height = 4)
correlation_factors
dev.off()

pdf(file = 'Results/Correlation_matrix_diagramms.pdf', width = 6.7, height = 6.7)
pairs(real_comparison[c("Shared_Edges_with_attributes_percent", 
                        "Density", "Num_Subnetworks", "Modularity_Coefficient", 
                        "Degree", "Betweenness", "Closeness")])
dev.off()

###

## ALPHA DIVERSITY ##############################################################################################

otu_table <- read.csv('Dataframes/all_otus.tsv', sep='\t') %>% 
  column_to_rownames(var = "Unnamed..0") %>% 
  setNames(gsub("^X", "", names(.))) 

keystone_list <- real_comparison %>%
  filter(Keystone_candidate == "yes") %>%
  pull(Dict_Network) %>%
  as.list()

sanitize_colnames <- function(colnames) {
  make.names(colnames, unique = TRUE)
}
calculate_diversity <- function(otu_table_in, keystone_list) {
  species_richness <- rowSums(otu_table_in > 0)
  shannon_index <- diversity(otu_table_in, index = "shannon")
  simpson_index <- diversity(otu_table_in, index = "simpson")
  
  diversity_df <- data.frame(
    Description = rownames(otu_table_in),
    Species_richness = species_richness,
    Shannon_index = shannon_index,
    Simpson_index = simpson_index
  )
  
  sanitized_otu_table_in <- otu_table_in
  colnames(sanitized_otu_table_in) <- sanitize_colnames(colnames(otu_table_in))
  keystone_list_sanitized <- sanitize_colnames(keystone_list)
  
  for (keystone in keystone_list_sanitized) {
    if (keystone %in% colnames(sanitized_otu_table_in)) {
      diversity_df[[keystone]] <- ifelse(sanitized_otu_table_in[[keystone]] > 0, "yes", "no")
    } else {
      diversity_df[[keystone]] <- "no"
    }
  }
  
  return(diversity_df)
}

# Calculate diversity and add keystone columns
diversity_df <- calculate_diversity(otu_table, keystone_list) %>% 
  setNames(gsub("^X", "", names(.))) 
write.csv(diversity_df, 'Results/Alpha_Diversity_df.csv')

# Generate and display plots
sanitize_colnames <- function(colnames) {
  colnames <- gsub(" ", "_", colnames)
  colnames <- gsub("[^[:alnum:]_]", "", colnames)
  colnames <- make.names(colnames, unique = TRUE)
  return(colnames)
}

plot_alpha_diversity <- function(diversity_df, keystone_list) {
  sanitized_keystones <- sanitize_colnames(keystone_list)
  sanitized_colnames <- sanitize_colnames(colnames(diversity_df))
  
  missing_cols <- setdiff(sanitized_keystones, sanitized_colnames)
  if (length(missing_cols) > 0) {
    stop("The following keystone columns are missing in the dataframe: ", paste(missing_cols, collapse = ", "))
  }
  
  names(diversity_df) <- sanitize_colnames(names(diversity_df))
  
  diversity_long <- melt(diversity_df, id.vars = c("Description", sanitized_keystones),
                         variable.name = "Diversity_Metric", value.name = "Value")
  
  plots <- lapply(sanitized_keystones, function(keystone) {
    ggplot(diversity_long, aes_string(x = keystone, y = "Value", fill = keystone)) +
      geom_boxplot() +
      facet_wrap(~Diversity_Metric, scales = "free_y") +
      labs(title = paste("Alpha Diversity for\n", keystone),
           x = keystone,
           y = "Diversity Index") +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 10))
  })
  
  for (i in seq_along(plots)) {
    keystone <- sanitized_keystones[i]
    plot_name <- paste0("alpha_div_", gsub(" ", "_", keystone), ".pdf")
    ggsave(paste0("Results/", plot_name), plots[[i]], device = "pdf")
  }
}

plot_alpha_diversity(diversity_df, keystone_list)

# Calculate p values
calculate_p_values <- function(diversity_df, keystone, diversity_metric) {
  # Filter data for the specified keystone and diversity metric
  filtered_data <- diversity_df %>%
    select(Description, !!sym(keystone), !!sym(diversity_metric))
  
  # Filter data for 'yes' and 'no' values in the keystone column
  x <- filtered_data %>%
    filter(!!sym(keystone) == "yes") %>%
    pull(!!sym(diversity_metric))
  
  y <- filtered_data %>%
    filter(!!sym(keystone) == "no") %>%
    pull(!!sym(diversity_metric))
  
  # Perform Wilcoxon test to calculate p-value
  if (length(x) > 0 && length(y) > 0) {
    p_value <- wilcox.test(x, y, paired = FALSE)$p.value
  } else {
    p_value <- NA
  }
  
  # Create a data frame to store the result
  p_value_df <- data.frame(
    Keystone = paste("X-", keystone, sep = ""),
    Diversity_Metric = diversity_metric,
    P_Value = p_value,
    with_keystone = length(x), 
    without_keystone = length(y)
  )
  
  return(p_value_df)
}

p_value_df <- data.frame(Keystone = character(), Diversity_Metric = character(), P_Value = numeric(),
                         with_keystone = numeric(), without_keystone = numeric())
keystone_list_sanitized <- sanitize_colnames(keystone_list)

diversity_metrics <- c("Species_richness", "Shannon_index", "Simpson_index")

for (keystone in keystone_list) {
  for (diversity_metric in diversity_metrics) {
    # Calculate p-value for the current keystone and diversity metric
    p_value <- calculate_p_values(diversity_df, keystone, diversity_metric)
    
    # Append the result to the p_value_df
    p_value_df <- rbind(p_value_df, p_value)
  }
}

p_value_df <- p_value_df %>%
  mutate(Significant = ifelse(P_Value < 0.05, "yes", "no"))
write.csv(p_value_df, 'Results/Alpha_Diversity_df_p_val.csv')

p_value_df

for (div_metric in diversity_metrics){
  
  clean_df <- p_value_df %>% filter(!is.na(P_Value), Diversity_Metric == div_metric)  %>%
    mutate(log_P_Value = -log10(P_Value))
  
  correlation_result <- cor.test(clean_df$with_keystone, clean_df$P_Value, method = "spearman")
  
  plot <- ggplot(clean_df, aes(x = with_keystone, y = log_P_Value)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") + # Adding a linear trend line
    labs(title = paste("Keystone Sample Nr. vs. -log10(P_Value) for\n", div_metric),
         subtitle = paste0("Spearman Correlation p-value: ", format(correlation_result$p.value, digits = 3)),
         x = "With Keystone",
         y = "-log10(P_Value)") +
    theme_minimal() 
  
  ggsave(filename = paste0("Results/scatter_plot_p_val_keystone_sample_nr_", gsub(" ", "_", div_metric), ".pdf"), plot = plot)
  
}


## BETA DIVERSITY ##############################################################################################
# For reproducability
set.seed(19980607)

names(diversity_df) <- sanitize_colnames(names(diversity_df))
keystone_list_sanitized <- names(diversity_df)[5:14]

p_value_df_beta <- data.frame(keystone = character(), pval = numeric(), stringsAsFactors = FALSE)

for (keystone_species in keystone_list_sanitized) {
  
  distmat_dist <- vegdist(otu_table, method = "bray", binary = FALSE)
  distmat_dudi <- dudi.pco(distmat_dist, scannf = FALSE)
  
  eigenvalues_con <- distmat_dudi$eig
  variance_proportion_con <- eigenvalues_con / sum(eigenvalues_con)
  
  axis1 <- round((variance_proportion_con[[1]] * 100), digits = 2)
  axis2 <- round((variance_proportion_con[[2]] * 100), digits = 2)
  
  # numeric_cols <- sapply(distmat_dist, is.numeric)
  # numeric_data <- distmat_dist[, numeric_cols]
  dist_df_bray <- as.data.frame(as.matrix(distmat_dist))
  
  dist_df_bray$Description <- rownames(dist_df_bray)
  
  meta_dist_df_bray <- inner_join(diversity_df, dist_df_bray, by = 'Description')

  # Create a new column with a standard name for the keystone species
  meta_dist_df_bray$keystone_temp <- factor(meta_dist_df_bray[[keystone_species]], levels = c("yes", "no"))
  
  dist_df_bray <- dist_df_bray %>% select(-Description)
  
  permutations <- 10000
  permanova_test_bray <- adonis2(dist_df_bray ~ keystone_temp, data = meta_dist_df_bray,  
                                 permutations = permutations, method = "bray")
  pval <- permanova_test_bray$'Pr(>F)'[1]
  p_value_df_beta <- rbind(p_value_df_beta, data.frame(keystone = keystone_species, pval = pval, stringsAsFactors = FALSE))
  
  plot_data <- as.data.frame(distmat_dudi$li) %>% 
    rownames_to_column(var = "Description") %>%
    inner_join(diversity_df, by = "Description") %>% 
    mutate(row.names = Description) %>%
    select(-Description) %>%
    column_to_rownames(var = "row.names") %>% 
    rename(x = A1, y = A2, group = keystone_species)
  
  plot_pcoa_bray_con <- ggplot(plot_data, aes(x = x, y = y, color = group)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
    geom_point(size = 1) +  
    stat_ellipse(aes(fill = group), geom = "polygon", alpha = 0.1, show.legend = FALSE, size = 0.3) +  # Ellipses
    theme_minimal() +
    labs(title = paste("PCoA"), color = "Keystone in samples", 
         subtitle = paste("Adonis2 p-value", round(pval, 6))) +
    xlab(paste0('PCoA Axis 1 (', axis1, '%)')) +
    ylab(paste0('PCoA Axis 2 (', axis2, '%)')) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 11), 
          plot.subtitle = element_text(hjust = 0.5, size = 8), legend.box = "vertical", 
          legend.margin = margin(-5, -5, -5, -5), legend.box.margin = margin(3, 3, 3, 3),
          legend.spacing.x = unit(0.1, "cm"),
          axis.title = element_text(size = 8)) +
    scale_shape_discrete(drop = FALSE)
  
  ggsave(filename = paste0("Results/beta_diversity_", gsub(" ", "_", keystone_species), ".pdf"), plot = plot_pcoa_bray_con)
}

p_value_df_beta <- p_value_df_beta %>%
  mutate(Significant = ifelse(pval < 0.05, "yes", "no")) %>% 
  mutate(keystone = str_replace(keystone, "^X", ""))

p_value_df_beta <- inner_join(p_value_df_beta, taxonomy, by = c('keystone' = 'Feature.ID'))
write.csv(p_value_df_beta, 'Results/Beta_Diversity_df_p_val.csv')



## Compare Alpha-Diversity-Analysis to random samples ####################################################################
# Function to sanitize column names
sanitize_colnames <- function(colnames) {
  colnames <- gsub(" ", "_", colnames)
  colnames <- gsub("[^[:alnum:]_]", "", colnames)
  colnames <- make.names(colnames, unique = TRUE)
  return(colnames)
}

# Function to calculate diversity metrics
calculate_diversity <- function(otu_table_in, otus_list) {
  species_richness <- rowSums(otu_table_in > 0)
  shannon_index <- diversity(otu_table_in, index = "shannon")
  simpson_index <- diversity(otu_table_in, index = "simpson")
  
  diversity_df <- data.frame(
    Description = rownames(otu_table_in),
    Species_richness = species_richness,
    Shannon_index = shannon_index,
    Simpson_index = simpson_index
  )
  
  sanitized_otu_table_in <- otu_table_in
  colnames(sanitized_otu_table_in) <- sanitize_colnames(colnames(otu_table_in))
  otus_list_sanitized <- sanitize_colnames(otus_list)
  
  for (otu in otus_list_sanitized) {
    if (otu %in% colnames(sanitized_otu_table_in)) {
      diversity_df[[otu]] <- ifelse(sanitized_otu_table_in[[otu]] > 0, "yes", "no")
    } else {
      diversity_df[[otu]] <- "no"
    }
  }
  
  return(diversity_df)
}

# Function to calculate p-values
calculate_p_values <- function(diversity_df, keystone_list, diversity_metric) {
  p_value_df <- data.frame(Keystone = character(), Diversity_Metric = character(), P_Value = numeric(),
                           with_keystone = numeric(), without_keystone = numeric(), stringsAsFactors = FALSE)
  
  for (keystone in keystone_list) {
    # Filter data for the specified keystone and diversity metric
    filtered_data <- diversity_df %>%
      select(Description, !!sym(keystone), !!sym(diversity_metric))
    
    # Filter data for 'yes' and 'no' values in the keystone column
    x <- filtered_data %>%
      filter(!!sym(keystone) == "yes") %>%
      pull(!!sym(diversity_metric))
    
    y <- filtered_data %>%
      filter(!!sym(keystone) == "no") %>%
      pull(!!sym(diversity_metric))
    
    # Perform Wilcoxon test to calculate p-value
    if (length(x) > 0 && length(y) > 0) {
      p_value <- wilcox.test(x, y, paired = FALSE)$p.value
    } else {
      p_value <- NA
    }
    
    # Create a data frame to store the result
    p_value_row <- data.frame(
      Keystone = paste("X-", keystone, sep = ""),
      Diversity_Metric = diversity_metric,
      P_Value = p_value,
      with_keystone = length(x), 
      without_keystone = length(y),
      stringsAsFactors = FALSE
    )
    
    # Append the result to the p_value_df
    p_value_df <- rbind(p_value_df, p_value_row)
  }
  
  p_value_df$Significant <- ifelse(p_value_df$P_Value < 0.05, "yes", "no")
  
  return(p_value_df)
}

# Main code block
set.seed(123) 
# Define the number of runs
num_runs <- 10
sample_size <- length(keystone_list)

p_value_df_random_all <- data.frame(Keystone = character(), Diversity_Metric = character(), 
                                    P_Value = numeric(), with_keystone = numeric(), 
                                    without_keystone = numeric(), Origin = character(),
                                    stringsAsFactors = FALSE)

for (run in 1:num_runs) {
  random_otus <- real_comparison %>%
    pull(Dict_Network) %>%
    sample(sample_size)
  
    diversity_df_random <- calculate_diversity(otu_table, random_otus) %>% 
    setNames(gsub("^X", "", names(.))) 
  
    p_value_df_random <- data.frame(Keystone = character(), Diversity_Metric = character(), 
                                  P_Value = numeric(), with_keystone = numeric(), 
                                  without_keystone = numeric(), Origin = character(),
                                  stringsAsFactors = FALSE)
  
  for (diversity_metric in diversity_metrics) {
    p_value_random <- calculate_p_values(diversity_df_random, random_otus, diversity_metric) %>% 
      mutate(Origin = paste0('Random_', run))
    
    p_value_df_random <- rbind(p_value_df_random, p_value_random)
  }
  
  p_value_df_random_all <- rbind(p_value_df_random_all, p_value_df_random)
}

write.csv(p_value_df_random_all, 'Results/Alpha_Diversity_df_p_val_random_all.csv')

significant_counts <- p_value_df_random_all %>%
  group_by(Origin, Diversity_Metric, Significant) %>%
  summarise(Count = n())
print(significant_counts)

significant_counts_real <- p_value_df %>%
  group_by(Diversity_Metric, Significant) %>%
  summarise(Count = n()) %>% 
  mutate(Origin = 'Keystone')

comparison_df <- significant_counts %>%
  left_join(significant_counts_real, by = c("Diversity_Metric", "Significant"), suffix = c("_random", "_real")) %>%
  mutate(Difference = Count_real - Count_random)

colnames(comparison_df)

summary_stats <- comparison_df %>%
  filter(Significant == 'yes') %>% 
  group_by(Diversity_Metric) %>%
  summarise(Mean_Difference = mean(Difference, na.rm = TRUE),
            SD_Difference =  sd(Difference, na.rm = TRUE),
            Min_Difference = min(Difference, na.rm = TRUE),
            Max_Difference = max(Difference, na.rm = TRUE))
summary_stats
write.csv(p_value_df_random_all, 'Results/Alpha_Diversity_random_subsets_real_comparison.csv')

## Compare Beta-Diversity-Analysis to random samples ####################################################################
calculate_beta_p_values <- function(otu_list, diversity_df, otu_table, taxonomy) {
  p_value_df <- data.frame(keystone = character(), pval = numeric(), stringsAsFactors = FALSE)
  
  for (keystone_species in otu_list) {
    
    distmat_dist <- vegdist(otu_table, method = "bray", binary = FALSE)
    distmat_dudi <- dudi.pco(distmat_dist, scannf = FALSE)
    
    dist_df_bray <- as.data.frame(as.matrix(distmat_dist))
    dist_df_bray$Description <- rownames(dist_df_bray)
    
    meta_dist_df_bray <- inner_join(diversity_df, dist_df_bray, by = 'Description')
    meta_dist_df_bray$keystone_temp <- factor(meta_dist_df_bray[[keystone_species]], levels = c("yes", "no"))
    
    dist_df_bray <- dist_df_bray %>% select(-Description)
    
    permutations <- 10000
    permanova_test_bray <- adonis2(dist_df_bray ~ keystone_temp, data = meta_dist_df_bray, permutations = permutations, method = "bray")
    pval <- permanova_test_bray$'Pr(>F)'[1]
    p_value_df <- rbind(p_value_df, data.frame(keystone = keystone_species, pval = pval, stringsAsFactors = FALSE))
  }
  
  p_value_df <- p_value_df %>%
    mutate(Significant = ifelse(pval < 0.05, "yes", "no")) %>% 
    mutate(keystone = str_replace(keystone, "^X", ""))
  
  return(p_value_df)
}


calculate_beta_p_values_random <- function(otu_list, diversity_df, otu_table, taxonomy) {
  p_value_df <- data.frame(keystone = character(), pval = numeric(), stringsAsFactors = FALSE)
  
  for (keystone_species in otu_list) {
    
    distmat_dist <- vegdist(otu_table, method = "bray", binary = FALSE)
    distmat_dudi <- dudi.pco(distmat_dist, scannf = FALSE)
    
    dist_df_bray <- as.data.frame(as.matrix(distmat_dist))
    dist_df_bray$Description <- rownames(dist_df_bray)
    
    sanitized_otu_table_in <- otu_table
    colnames(sanitized_otu_table_in) <- sanitize_colnames(colnames(otu_table))
    otus_list_sanitized <- sanitize_colnames(otu_list)
    
    for (otu in otus_list_sanitized) {
      if (otu %in% colnames(sanitized_otu_table_in)) {
        diversity_df[[otu]] <- ifelse(sanitized_otu_table_in[[otu]] > 0, "yes", "no")
      } else {
        diversity_df[[otu]] <- "no"
      }}
    
    meta_dist_df_bray <- inner_join(diversity_df, dist_df_bray, by = 'Description')
    meta_dist_df_bray$keystone_temp <- factor(meta_dist_df_bray[[keystone_species]], levels = c("yes", "no"))
    
    dist_df_bray <- dist_df_bray %>% select(-Description)
    
    permutations <- 10000
    permanova_test_bray <- adonis2(dist_df_bray ~ keystone_temp, data = meta_dist_df_bray, permutations = permutations, method = "bray")
    pval <- permanova_test_bray$'Pr(>F)'[1]
    p_value_df <- rbind(p_value_df, data.frame(keystone = keystone_species, pval = pval, stringsAsFactors = FALSE))
  }
  
  p_value_df <- p_value_df %>%
    mutate(Significant = ifelse(pval < 0.05, "yes", "no")) %>% 
    mutate(keystone = str_replace(keystone, "^X", ""))
  
  return(p_value_df)
}

# Set seed for reproducibility
set.seed(19980607)

# Sanitize column names
names(diversity_df) <- sanitize_colnames(names(diversity_df))
keystone_list <- real_comparison %>%
  filter(Keystone_candidate == "yes") %>%
  pull(Dict_Network) %>%
  as.list()
keystone_list_sanitized <- sanitize_colnames(keystone_list)
names(otu_table) <- sanitize_colnames(names(otu_table))

p_value_df_beta_real <- calculate_beta_p_values(otu_list = keystone_list_sanitized, diversity_df, otu_table, taxonomy)

num_random_samples <- length(keystone_list)
p_value_df_beta_random <- data.frame(keystone = character(), pval = numeric(), Significant = character(), Feature.ID = character(), Taxonomy = character(), Origin = character(), stringsAsFactors = FALSE)

for (i in 1:10) {
  random_otus <- sample(names(otu_table), num_random_samples)
  random_otus <- sanitize_colnames(random_otus)
  
  p_value_df_random <- calculate_beta_p_values_random(random_otus, diversity_df, otu_table, taxonomy) %>%
    mutate(Origin = paste0("Random_", i))
  
  p_value_df_beta_random <- rbind(p_value_df_beta_random, p_value_df_random)
}

p_value_df_beta_real <- p_value_df_beta_real %>% 
  mutate_all(~ coalesce(as.character(.), "NA")) %>% 
  group_by(Significant) %>%
  summarise(Count = n())

p_value_df_beta_random <- p_value_df_beta_random %>% 
  mutate_all(~ coalesce(as.character(.), "NA")) %>% 
  group_by(Significant, Origin) %>%
  summarise(Count = n())

comparison_df <- p_value_df_beta_random %>%
  left_join(p_value_df_beta_real, by = c("Significant"), suffix = c("_random", "_real"))

comparison_df <- comparison_df %>%
  mutate(Difference = Count_real - Count_random)

summary_stats <- comparison_df %>%
  filter(Significant == 'yes') %>% 
  group_by(Significant) %>%
  summarise(Mean_Difference = mean(Difference, na.rm = TRUE),
            SD_Difference =  sd(Difference, na.rm = TRUE),
            Min_Difference = min(Difference, na.rm = TRUE),
            Max_Difference = max(Difference, na.rm = TRUE))
summary_stats
write.csv(p_value_df_random_all, 'Results/Beta_Diversity_random_subsets_real_comparison.csv')

colnames(real_comparison)
## Testing the keystoneness rank from other papers
real_comparison <- real_comparison %>% 
  mutate(keysone_indicator = Degree + Closeness - Betweenness)
  
