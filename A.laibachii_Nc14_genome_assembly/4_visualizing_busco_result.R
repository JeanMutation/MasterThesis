library(ggplot2)
library(reshape2)
library(tidyverse)
library(scales)

# Import the BUSCO table

busco_df <- read.csv("quality_assesment_assembly/busco.csv", header = TRUE)

# Organize and rearrange the imported table

busco_df$Strain <- as.factor(busco_df$Strain)
busco_df.melted <- melt(busco_df, id.vars = "Strain")
busco_df.melted$variable <-relevel(busco_df.melted$variable, "Missing")

# Create a stacked bar plot for the BUSCO outputs

busco_plot <- ggplot(busco_df.melted, aes(x=Strain, fill=fct_rev(variable), y=value)) +
  geom_bar(position= "stack", width = 0.7, stat="identity") +
  labs(x = "Strain", y = "BUSCO", fill = "Type") +
  scale_y_continuous(labels=comma) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=12))

# Save the plot as “busco.pdf”

pdf("quality_assesment_assembly/busco.pdf",width=8,height=5,paper='special')

print(busco_plot)

dev.off()
