setwd("D:/master_project_networks/workflows/otu_table")

library(ggplot2)
library(ggrepel)
library(gtable)
library(readr)
library(gridExtra)
library(tidyr)
library("dplyr") 
library(stringr)

# Read in the otu tables and the taxonomy tables
net_analysis_wl <- read_csv("network_output_fastspar/wl_nw_fp10/wl_fp10_p005_corNOT0_nw_node_preprocessed.csv")
net_analysis_endo <- read_csv("network_output_fastspar/endo_nw_fp10/endo_fp10_p005_corNOT0_nw_node_preprocessed.csv")
tax_wl <- read.table('raw_otu_and_tax/Approach6/taxonomy.tsv', header=TRUE, sep='\t')
tax_endo <- read.table('raw_otu_and_tax/Approach6/taxonomy.tsv', header=TRUE, sep='\t')

# Setting the quantile threshold
n <- 5
## Calculation for wl ####################
top5_BC_wl <- net_analysis_wl[net_analysis_wl$BetweennessCentrality > 
                                quantile(net_analysis_wl$BetweennessCentrality, prob=1-n/100),]$name
top5_BC_wl_min <- min(net_analysis_wl[net_analysis_wl$BetweennessCentrality > 
                                quantile(net_analysis_wl$BetweennessCentrality,prob=1-n/100),]$BetweennessCentrality)

top5_CC_wl <- net_analysis_wl[net_analysis_wl$ClosenessCentrality > 
                                quantile(net_analysis_wl$ClosenessCentrality,prob=1-n/100),]$name
top5_CC_wl_min<-min(net_analysis_wl[net_analysis_wl$ClosenessCentrality > 
                               quantile(net_analysis_wl$ClosenessCentrality,prob=1-n/100),]$ClosenessCentrality)

# top5_deg_wl <- net_analysis_wl[net_analysis_wl$Degree > 
#                                 quantile(net_analysis_wl$degree,prob=1-n/100),]$name
# top5_deg_wl_min <- min(net_analysis_wl[net_analysis_wl$degree > 
#                                       quantile(net_analysis_wl$degree,prob=1-n/100),]$Degree)

p5_Intersection_BC_CC_wl <- top5_BC_wl[which((top5_BC_wl %in% top5_CC_wl) == TRUE)]
net_analysis_wl$hub_BC_CC <- ifelse(net_analysis_wl$name %in% p5_Intersection_BC_CC_wl,'yes','no')

## Calculation for endo ######################
top5_BC_endo <- net_analysis_endo[net_analysis_endo$BetweennessCentrality > 
                                quantile(net_analysis_endo$BetweennessCentrality, prob=1-n/100),]$name
top5_BC_endo_min <- min(net_analysis_endo[net_analysis_endo$BetweennessCentrality > 
                                        quantile(net_analysis_endo$BetweennessCentrality,prob=1-n/100),]$BetweennessCentrality)

top5_CC_endo <- net_analysis_endo[net_analysis_endo$ClosenessCentrality > 
                                quantile(net_analysis_endo$ClosenessCentrality,prob=1-n/100),]$name
top5_CC_endo_min<-min(net_analysis_endo[net_analysis_endo$ClosenessCentrality > 
                                      quantile(net_analysis_endo$ClosenessCentrality,prob=1-n/100),]$ClosenessCentrality)

top5_deg_endo <- net_analysis_endo[net_analysis_endo$Degree > 
                                 quantile(net_analysis_endo$Degree,prob=1-n/100),]$name
top5_deg_endo_min <- min(net_analysis_endo[net_analysis_endo$Degree > 
                                         quantile(net_analysis_endo$Degree,prob=1-n/100),]$Degree)

p5_Intersection_endo <- top5_BC_endo[which((top5_BC_endo %in% top5_CC_endo) == TRUE)]
net_analysis_endo$hub_BC_CC <- ifelse(net_analysis_endo$name %in% p5_Intersection_endo,'yes','no')

## names for the plots # WL ###############
tax_wl <- tax_wl %>% select(c('Feature.ID', 'Genus')) %>% rename('name' = 'Feature.ID')

net_analysis_wl$name <- gsub("-", "_", net_analysis_wl$name)

net_analysis_wl <- inner_join(net_analysis_wl, tax_wl, by='name')
net_analysis_wl <- mutate(net_analysis_wl, name_graph = paste(name, Genus, sep = "-"))
net_analysis_wl$Kingdom <- sub("^[^_]+_", "", net_analysis_wl$name)
net_analysis_wl <- as.data.frame(net_analysis_wl)


## names for plots # ENDO #################
tax_endo <- tax_endo %>% select(c('Feature.ID', 'Genus')) %>% rename('name' = 'Feature.ID')

net_analysis_endo <- inner_join(net_analysis_endo, tax_endo, by='name')
net_analysis_endo <- mutate(net_analysis_endo, name_graph = paste(name, Genus, sep = "-"))
net_analysis_endo$Kingdom <- sub("^[^_]+_", "", net_analysis_endo$name)
net_analysis_endo <- as.data.frame(net_analysis_endo)


## plot the hubs
wl_BC_CC <- ggplot(net_analysis_wl, aes(x = ClosenessCentrality, y= BetweennessCentrality , label=Genus, shape = Kingdom)) +
  geom_point(aes(colour = Kingdom) ,stat= "identity",position="identity",alpha= 1 ) +
  scale_colour_manual(values = c('#8c564b', '#7f7f7f', '#ff7f0e', '#e377c2' )) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +  # Adjust shapes accordingly
  geom_vline(xintercept = top5_CC_wl_min, linetype="dashed", linewidth = 0.3) + 
  geom_hline(yintercept = top5_BC_wl_min, linetype= "dashed", linewidth = 0.3) + 
  labs(title = "Hubs - whole leaf samples") +
  theme_minimal() + 
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom") +
  geom_text_repel(data=subset(net_analysis_wl, hub_BC_CC =="yes")
                  ,size = 3,
                  box.padding   = 1,
                  point.padding = 0.2,
                  force         = 50,
                  segment.size  = 0.2,
                  segment.color = "darkgreen") 

wl_BC_CC

net_analysis_wl$print_name <- net_analysis_wl$hub_BC_CC
net_analysis_wl$print_name <- ifelse(grepl("Albugo", net_analysis_wl$Genus), "yes", net_analysis_wl$print_name)
net_analysis_wl$print_name <- ifelse(grepl("Dioszegia", net_analysis_wl$Genus), "yes", net_analysis_wl$print_name)

wl_BC_CC <- ggplot(net_analysis_wl, aes(x = ClosenessCentrality, y= BetweennessCentrality , label=Genus, shape = Kingdom)) +
  geom_point(aes(colour = Kingdom) ,stat= "identity",position="identity",alpha= 1 ) +
  scale_colour_manual(values = c('#8c564b', '#7f7f7f', '#ff7f0e', '#e377c2' )) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +  # Adjust shapes accordingly
  geom_vline(xintercept = top5_CC_wl_min, linetype="dashed", linewidth = 0.3) + 
  geom_hline(yintercept = top5_BC_wl_min, linetype= "dashed", linewidth = 0.3) + 
  labs(title = "Hubs - whole leaf samples") +
  theme_minimal() + 
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom") +
  geom_text_repel(data=subset(net_analysis_wl, print_name =="yes")
                  ,size = 3,
                  box.padding   = 1,
                  max.overlaps = 20,
                  point.padding = 0.2,
                  force         = 50,
                  segment.size  = 0.2,
                  segment.color = "darkgreen") 

wl_BC_CC

#---------------------------------------------------------------------------------------------------------------------------#
endo_BC_CC <- ggplot(net_analysis_endo, aes(x = ClosenessCentrality, y= BetweennessCentrality , label=Genus, shape = Kingdom)) +
  geom_point(aes(colour = Kingdom) ,stat= "identity",position="identity",alpha= 1 ) +
  scale_colour_manual(values = c('#8c564b', '#7f7f7f', '#ff7f0e', '#e377c2' )) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +  # Adjust shapes accordingly
  geom_vline(xintercept = top5_CC_endo_min, linetype="dashed", linewidth = 0.3) + 
  geom_hline(yintercept = top5_BC_endo_min, linetype= "dashed", linewidth = 0.3) + 
  labs(title = "Hubs - endophyte samples") +
  theme_minimal() + 
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom") +
  theme(plot.caption = element_text(hjust = 0.5, margin = margin(b = 10, unit = "pt"))) +
  geom_text_repel(data=subset(net_analysis_endo, hub_BC_CC =="yes")
                  ,size = 3,
                  box.padding   = 1,
                  point.padding = 0.2,
                  force         = 50,
                  segment.size  = 0.2,
                  segment.color = "darkgreen")

endo_BC_CC

net_analysis_endo$print_name <- net_analysis_endo$hub_BC_CC
net_analysis_endo$print_name <- ifelse(grepl("Albugo", net_analysis_endo$Genus), "yes", net_analysis_endo$print_name)
net_analysis_endo$print_name <- ifelse(grepl("Dioszegia", net_analysis_endo$Genus), "yes", net_analysis_endo$print_name)

endo_BC_CC <- ggplot(net_analysis_endo, aes(x = ClosenessCentrality, y= BetweennessCentrality , label=Genus, shape = Kingdom)) +
  geom_point(aes(colour = Kingdom) ,stat= "identity",position="identity",alpha= 1 ) +
  scale_colour_manual(values = c('#8c564b', '#7f7f7f', '#ff7f0e', '#e377c2' )) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +  # Adjust shapes accordingly
  geom_vline(xintercept = top5_CC_endo_min, linetype="dashed", linewidth = 0.3) + 
  geom_hline(yintercept = top5_BC_endo_min, linetype= "dashed", linewidth = 0.3) + 
  labs(title = "Hubs - endophyte samples") +
  theme_minimal() + 
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom") +
  theme(plot.caption = element_text(hjust = 0.5, margin = margin(b = 10, unit = "pt"))) +
  geom_text_repel(data=subset(net_analysis_endo, print_name =="yes")
                  ,size = 3,
                  box.padding   = 0.2,
                  point.padding = 0.2,
                  force         = 100,
                  segment.size  = 0.4,
                  segment.color = "darkgreen")

endo_BC_CC

#---Save-the-plots--------------------------------------------------------------------------------------------------#
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <-g_legend(endo_BC_CC)

pdf(file = 'D:/master_project_networks/workflows/otu_table/div_plots_final/hub_analysis.pdf',
    width = 6.7,
    height = 4)

grid.arrange(arrangeGrob(endo_BC_CC + theme(legend.position="none"), wl_BC_CC + theme(legend.position="none"),
                         ncol = 2),
             mylegend, nrow=2,heights=c(10, 1))

dev.off()
