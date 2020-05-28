# This script draws boxplots showing enrichment for rare variants and NMD variants
# amoung single tissue eOutliers, sOutliers, and aseOutliers
# 
# Jonah Einson
# jeinson@nygenome.org

library(tidyverse)
library(ggpubr)
library(cowplot)

# Stuff just for this project
class_colors <- c(
  ASE=rgb(121, 92, 129, maxColorValue = 255), 
  Splicing=rgb(23, 50, 75, maxColorValue = 255), 
  TotalExpression=rgb(193, 205, 222, maxColorValue = 255)
)

gtex_v8_figure_theme <- function() {
  return(
    theme(
      plot.title = element_text(face="plain",size=8), 
      text = element_text(size=8),
      axis.text=element_text(size=7), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"), 
      legend.text = element_text(size=7), 
      legend.title = element_text(size=8)
    )
  )
}

#### Panel A: Enrichment of rare variants in single tissue outlirs 
results <- read_tsv("processed_input_data/figureS23/single_tissue_variant_enrichment_v8_Vg.tsv")
outlier_threshold <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9)
results$sig <- factor(results$sig, levels = outlier_threshold)

names(class_colors) <- c("ASE", "Splicing", "TotalExpression")
panel_a <- 
  ggplot(results, aes(sig, ratio, fill = outlier_type)) +
  geom_boxplot(size = .3, outlier.size = .3)  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = class_colors, labels = c("aseOutlier", "sOutlier", "eOutlier")) + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("FDR-corrected p-value outlier threshold") + ylab("Relative Risk") +
  theme(legend.position = c(.2, .9)) +
  gtex_v8_figure_theme() +
  theme(legend.title = element_blank())


#### Panel B: Enrichment of rare variants with NMD annotations in single tissue outliers ####
dat <- read_tsv("processed_input_data/figureS23/single_tissue_NMD_variant_enrichment_v8_Vg.tsv", comment = "#")

dat <- 
  dat %>%
  group_by(anno) %>%
  mutate(med_ratio = median(ratio)) %>%
  ungroup %>%
  arrange(med_ratio) %>%
  mutate(anno = factor(anno, levels = rev(unique(anno))), 
         outlier_type = factor(outlier_type, levels = c( "TotalExpression", "ASE", "Splicing"))) 

panel_b <- 
  ggplot(dat, aes(x = anno, y = ratio, fill = outlier_type)) +
  geom_boxplot(size = .3, outlier.size = .3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(), 
        legend.position = "none") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = class_colors, labels = c("aseOutlier", "sOutlier", "eOutlier")) +
  ylab("Relative Risk") +
  yscale("log2", .format = F) +
  gtex_v8_figure_theme()



#### Combine with cowplot ####
png("generated_figures/figS23.png", width = 6, height = 4, units = "in", res = 300)
plot_grid(panel_a, panel_b, ncol = 2, align = "h")
dev.off()
