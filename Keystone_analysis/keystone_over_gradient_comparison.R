library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(grid)
  
setwd("D:/master_keystone_project/workflows/")

keystone_collection <- read.csv(paste0('keystone_analysis_run_1_data/Results/Real_Comparison.csv')) %>% 
  select(c(Dict_Network, Keystone_candidate))

colnames(keystone_collection) <- c('Dict_Network', paste0('Keystone_candidate_run_1'))

for(i in 9:17){
  real_comparison_df <- read.csv(paste0('keystone_analysis_run_', i, '_data/Results/Real_Comparison.csv'))

  real_comparison_df <- real_comparison_df %>%
    select(c(Dict_Network, Keystone_candidate)) 
    
  colnames(real_comparison_df) <- c('Dict_Network', paste0('Keystone_candidate_run_', i))
  
  keystone_collection <- merge(keystone_collection, real_comparison_df, by='Dict_Network', all = TRUE)
}

OTU_Numbers <- colSums(!is.na(keystone_collection))
Keystone_Numbers <- colSums(keystone_collection == 'yes', na.rm = TRUE)

run_parameters <- data.frame(
  run = c(1, 9:17),
  filter_param = c(10, seq(2, 8, 2), seq(12, 20, 2))
) %>% mutate(across(everything(), as.character))

summary <- as.data.frame(rbind(OTU_Numbers, Keystone_Numbers)) %>% 
  select(-Dict_Network) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(., var='Keystone') %>% 
  mutate(run = (str_split_i(Keystone, '_', 4))) %>% 
  inner_join(run_parameters, by='run') %>% 
  mutate(OTU_Numbers = as.numeric(OTU_Numbers),
         Keystone_Numbers = as.numeric(Keystone_Numbers),
         filter_param = as.numeric(filter_param))

sparsity_results <- data.frame(filter_param = integer(), sparsity = numeric(), stringsAsFactors = FALSE)

for (i in seq(2,20,2)){
  table <- read.csv(paste0('input_flashweave_complete/wl_otu_fp', i, 'p.tsv'), sep = '\t')
  
  total_elements <- prod(dim(table))
  zero_na_elements <- sum(table == 0, na.rm = TRUE) + sum(is.na(table))
  sparsity <- zero_na_elements / total_elements
  
  sparsity_results <- rbind(sparsity_results, data.frame(filter_param = i, sparsity = sparsity))
}

summary <- inner_join(summary, sparsity_results, by='filter_param')
  
keystone_num_OTU_num <- ggplot(summary, aes(x=OTU_Numbers, y=Keystone_Numbers)) +
  geom_point(size = 2, colour = "slateblue") +
  labs(x = "Number of OTUs", y = "Nr. of Keystone Genera") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        axis.title.y = element_text(size = 9))

keystone_num_filterparam <- ggplot(summary, aes(x=filter_param, y=Keystone_Numbers)) +
  geom_point(size = 2, colour = "slateblue") +
  labs(x = "Filter paramenter [%]", y = "") +
  scale_x_reverse() +
  theme(axis.title.y = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        axis.title.y = element_text(size = 9))

keystone_num_sparsity <- ggplot(summary, aes(x=sparsity, y=Keystone_Numbers)) +
  geom_point(size = 2, colour = "slateblue") +
  labs(x = "Sparsity of the OTU table", y = "") +
  theme(axis.title.y = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        axis.title.y = element_text(size = 9))


title <- textGrob("Number of Keystones in whole leaf networks", gp = gpar(fontsize = 12, fontface = "bold"))

pdf(file = 'final_plots/KeystoneNumbers_real_wl.pdf', width = 6.7, height = 2)
grid.arrange(title, 
             arrangeGrob(keystone_num_filterparam, keystone_num_OTU_num, keystone_num_sparsity, ncol = 3), 
             ncol = 1, heights = c(0.1, 0.9))
dev.off()


keystone_collection_short <- keystone_collection %>% 
  column_to_rownames(., var = 'Dict_Network') %>% 
  t() %>% 
  as.data.frame() %>% 
  select_if(~ any(. == "yes")) %>% 
  rownames_to_column(., var = 'Keystone')

taxonomy <- read.table("taxonomy.tsv", sep="\t", h=T) %>% 
  select(c(Feature.ID, Order, Family, Genus))
  
summary_with_keystone <- inner_join(summary, keystone_collection_short, by = 'Keystone') %>% 
  replace_na(as.list(setNames(rep("NA", ncol(.)), names(.))))

summary_with_keystone_long_fp <- summary_with_keystone %>%
  select(-c(OTU_Numbers, Keystone_Numbers, run, Keystone, sparsity)) %>% 
  pivot_longer(cols = !filter_param, names_to = "Keystone_ID", values_to = "value") %>% 
  left_join(., taxonomy, by = c('Keystone_ID' = 'Feature.ID')) %>% 
  mutate(name = paste0('f__', Family, ',g__', Genus),
         value = as.factor(value))

wl_full <- ggplot(summary_with_keystone_long_fp, aes(x = filter_param, y = name, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c('yes' = 'chartreuse2', 'no' = 'coral1', 'NA' = 'gray60' )) +
  scale_x_continuous(breaks = seq(2, 20, by = 2)) +
  labs(x = "Filter Parameter", y = " ", fill = "Keystone?") +
  theme_minimal() +
  theme(legend.position="bottom")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <-g_legend(wl_full)
wl_full

pdf(file = 'final_plots/KeystoneGenera_fullGradient_wl.pdf', width = 6.7, height = 6.7)
grid.arrange(wl_full + theme(legend.position="none"), mylegend,
             nrow = 2, heights = c(0.95, 0.05))
dev.off()

summary_with_keystone_long_fp_6t020 <- summary_with_keystone %>%
  select(-c(OTU_Numbers, Keystone_Numbers, run, Keystone)) %>%
  filter(filter_param > 8) %>% 
  mutate(filter_param = as.numeric(filter_param)) %>%
  select(filter_param, everything()) %>%
  select(filter_param, where(~ any(. == "yes"))) %>%
  pivot_longer(cols = !filter_param, names_to = "Keystone_ID", values_to = "value") %>% 
  left_join(., taxonomy, by = c('Keystone_ID' = 'Feature.ID')) %>% 
  mutate(name = paste0('f__', Family, ',g__', Genus),
         value = as.factor(value))

wl_filtered <- ggplot(summary_with_keystone_long_fp_6t020, aes(x = filter_param, y = name, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c('yes' = 'chartreuse2', 'no' = 'coral1', 'NA' = 'gray60' )) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  labs(x = "Filter Parameter", y = " ", fill = "Keystone?") +
  theme_minimal() +
  theme(legend.position="bottom")

mylegend <-g_legend(wl_filtered)

pdf(file = 'final_plots/KeystoneGenera_filteredGradient_wl.pdf', width = 6.7, height = 3)
grid.arrange(wl_filtered + theme(legend.position="none"), mylegend,
             nrow = 2, heights = c(0.90, 0.1))
dev.off()

##########################################################################
# SAME FOR ENDOPHYTES ####################################################
##########################################################################
setwd("D:/master_keystone_project/workflows/")

keystone_collection <- read.csv('keystone_analysis_run_5_data/Results/Real_Comparison.csv') %>% 
  select(c(Dict_Network, Keystone_candidate))

colnames(keystone_collection) <- c('Dict_Network', paste0('Keystone_candidate_run_5'))

for(i in 18:26){
  real_comparison_df <- read.csv(paste0('keystone_analysis_run_', i, '_data/Results/Real_Comparison.csv'))
  
  real_comparison_df <- real_comparison_df %>%
    select(c(Dict_Network, Keystone_candidate)) 
  
  colnames(real_comparison_df) <- c('Dict_Network', paste0('Keystone_candidate_run_', i))
  
  keystone_collection <- merge(keystone_collection, real_comparison_df, by='Dict_Network', all = TRUE)
}

OTU_Numbers <- colSums(!is.na(keystone_collection))
Keystone_Numbers <- colSums(keystone_collection == 'yes', na.rm = TRUE)

run_parameters <- data.frame(
  run = c(5, 18:26),
  filter_param = c(10, seq(2, 8, 2), seq(12, 20, 2))
) %>% mutate(across(everything(), as.character))

summary <- as.data.frame(rbind(OTU_Numbers, Keystone_Numbers)) %>% 
  select(-Dict_Network) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(., var='Keystone') %>% 
  mutate(run = (str_split_i(Keystone, '_', 4))) %>% 
  inner_join(run_parameters, by='run') %>% 
  mutate(OTU_Numbers = as.numeric(OTU_Numbers),
         Keystone_Numbers = as.numeric(Keystone_Numbers),
         filter_param = as.numeric(filter_param))

sparsity_results <- data.frame(filter_param = integer(), sparsity = numeric(), stringsAsFactors = FALSE)

for (i in seq(2,20,2)){
  table <- read.csv(paste0('input_flashweave_complete/endo_otu_fp', i, 'p.tsv'), sep = '\t')
  
  total_elements <- prod(dim(table))
  zero_na_elements <- sum(table == 0, na.rm = TRUE) + sum(is.na(table))
  sparsity <- zero_na_elements / total_elements
  
  sparsity_results <- rbind(sparsity_results, data.frame(filter_param = i, sparsity = sparsity))
}

summary <- inner_join(summary, sparsity_results, by='filter_param')

keystone_num_OTU_num <- ggplot(summary, aes(x=OTU_Numbers, y=Keystone_Numbers)) +
  geom_point(size = 2, colour = "slateblue") +
  labs(x = "Number of OTUs", y = "Nr. of Keystone Genera") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        axis.title.y = element_text(size = 9))

keystone_num_filterparam <- ggplot(summary, aes(x=filter_param, y=Keystone_Numbers)) +
  geom_point(size = 2, colour = "slateblue") +
  labs(x = "Filter paramenter [%]", y = "") +
  scale_x_reverse() +
  theme(axis.title.y = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        axis.title.y = element_text(size = 9))

keystone_num_sparsity <- ggplot(summary, aes(x=sparsity, y=Keystone_Numbers)) +
  geom_point(size = 2, colour = "slateblue") +
  labs(x = "Sparsity of the OTU table", y = "") +
  theme(axis.title.y = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        axis.title.y = element_text(size = 9))


title <- textGrob("Number of Keystones in endophytic networks", gp = gpar(fontsize = 12, fontface = "bold"))

pdf(file = 'final_plots/KeystoneNumbers_real_endo.pdf', width = 6.7, height = 2)
grid.arrange(title, 
             arrangeGrob(keystone_num_filterparam, keystone_num_OTU_num, keystone_num_sparsity, ncol = 3), 
             ncol = 1, heights = c(0.1, 0.9))
dev.off()

keystone_collection_short <- keystone_collection %>% 
  column_to_rownames(., var = 'Dict_Network') %>% 
  t() %>% 
  as.data.frame() %>% 
  select_if(~ any(. == "yes")) %>% 
  rownames_to_column(., var = 'Keystone')

taxonomy <- read.table("taxonomy.tsv", sep="\t", h=T) %>% 
  select(c(Feature.ID, Order, Family, Genus))

summary_with_keystone <- inner_join(summary, keystone_collection_short, by = 'Keystone') %>% 
  replace_na(as.list(setNames(rep("NA", ncol(.)), names(.))))

summary_with_keystone_long_fp <- summary_with_keystone %>%
  select(-c(OTU_Numbers, Keystone_Numbers, run, Keystone, sparsity)) %>% 
  pivot_longer(cols = !filter_param, names_to = "Keystone_ID", values_to = "value") %>% 
  left_join(., taxonomy, by = c('Keystone_ID' = 'Feature.ID')) %>% 
  mutate(name = paste0('f__', Family, ',g__', Genus),
         value = as.factor(value))

wl_full <- ggplot(summary_with_keystone_long_fp, aes(x = filter_param, y = name, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c('yes' = 'chartreuse2', 'no' = 'coral1', 'NA' = 'gray60' )) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  labs(x = "Filter Parameter", y = " ", fill = "Keystone?") +
  theme_minimal() +
  theme(legend.position="bottom")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <-g_legend(wl_full)

pdf(file = 'final_plots/KeystoneGenera_fullGradient_endo.pdf', width = 6.7, height = 6.7)
grid.arrange(wl_full + theme(legend.position="none"), mylegend,
             nrow = 2, heights = c(0.95, 0.05))
dev.off()

summary_with_keystone_long_fp_6t020 <- summary_with_keystone %>%
  select(-c(OTU_Numbers, Keystone_Numbers, run, Keystone)) %>%
  filter(filter_param > 8) %>% 
  mutate(filter_param = as.numeric(filter_param)) %>%
  select(filter_param, everything()) %>%
  select(filter_param, where(~ any(. == "yes"))) %>%
  pivot_longer(cols = !filter_param, names_to = "Keystone_ID", values_to = "value") %>% 
  left_join(., taxonomy, by = c('Keystone_ID' = 'Feature.ID')) %>% 
  mutate(name = paste0('f__', Family, ',g__', Genus),
         value = as.factor(value))

wl_filtered <- ggplot(summary_with_keystone_long_fp_6t020, aes(x = filter_param, y = name, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c('yes' = 'chartreuse2', 'no' = 'coral1', 'NA' = 'gray60' )) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  labs(x = "Filter Parameter", y = " ", fill = "Keystone?") +
  theme_minimal() +
  theme(legend.position="bottom")

###
wl_filtered <- ggplot(summary_with_keystone_long_fp_6t020, aes(x = filter_param, y = name, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c('yes' = 'chartreuse2', 'no' = 'coral1', 'NA' = 'gray60' )) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  labs(x = "Filter Parameter", y = " ", fill = "Keystone?") +
  theme_minimal() +
  theme(legend.position="bottom")

# Create annotations for the additional axis titles
sparcity_text <- textGrob("Sparcity", x = unit(0.5, "npc"), y = unit(0.95, "npc"), gp = gpar(fontsize = 12))
otu_number_text <- textGrob("OTU Number", x = unit(0.5, "npc"), y = unit(0.90, "npc"), gp = gpar(fontsize = 12))

# Combine the plot and the annotations
grid.arrange(
  wl_filtered,
  arrangeGrob(sparcity_text, otu_number_text, ncol = 1),
  heights = c(1, 0.1, 0.1)
)

mylegend <-g_legend(wl_filtered)

pdf(file = 'final_plots/KeystoneGenera_filteredGradient_endo.pdf', width = 6.7, height = 3)
grid.arrange(wl_filtered + theme(legend.position="none"), mylegend,
             nrow = 2, heights = c(0.90, 0.1))
dev.off()





#################
sparsity_results <- data.frame(Compartment = character(), sparsity = numeric(), stringsAsFactors = FALSE)

for (item in c('Root', 'Soil', 'RP', 'RS')){
  table <- read.csv(paste0('input_flashweave_complete/Root_', item, '_otu_fp10p.tsv'), sep = '\t')
  
  total_elements <- prod(dim(table))
  zero_na_elements <- sum(table == 0, na.rm = TRUE) + sum(is.na(table))
  sparsity <- zero_na_elements / total_elements
  
  sparsity_results <- rbind(sparsity_results, data.frame(Compartment = item, sparsity = sparsity))
}
sparsity_results
