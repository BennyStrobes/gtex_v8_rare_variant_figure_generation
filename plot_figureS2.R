# ANEVA-DOT summary stats figure
# Jonah Einson
# jeinson@nygenome.org

library(tidyverse)
library(magrittr)
library(ggpubr)
gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}
gtex_colors <- read_tsv("processed_input_data/gtex_colors.txt")
attach(gtex_colors)
names(tissue_site_detail) <- tissue_abbrv

############# Panel A. Dots per indv per tissue ##########
dots_per_indv_per_tiss <- read_rds("processed_input_data/figureS2/dots_per_indv_per_tiss_v8_Vg.rds")
genes_w_Vg_per_tiss <- read_rds("processed_input_data/figureS2/genes_w_Vg_per_tiss.rds")
median_dots_per_indv_per_tiss <- map_dbl(dots_per_indv_per_tiss, ~ median(.x$count))

# Save some things that will be useful 
tiss_order <- levels(genes_w_Vg_per_tiss$TissueID)
text_positions <- genes_w_Vg_per_tiss$genes_w_Vg_per_tiss + 600

Vg_coverage_tbl <- genes_w_Vg_per_tiss
Vg_coverage_tbl$`Median genes tested \nper indvididual` <- 
  median_dots_per_indv_per_tiss[tiss_order]
Vg_coverage_tbl$`Median genes tested \nper indvididual`[is.na(Vg_coverage_tbl$`Median genes tested \nper indvididual`)] <- 0

coverage_range <- 
  map_dfr(dots_per_indv_per_tiss, ~ range(.x$count)) %>% 
  t %>% as.data.frame %>%
  rownames_to_column %>%
  set_colnames(c("TissueID", "lower", "upper"))

indvs_covered <- 
  map_dfr(dots_per_indv_per_tiss, nrow) %>%
  t %>% as.data.frame %>%
  rownames_to_column %>%
  set_colnames(c("TissueID", "total"))

Vg_coverage_tbl %<>% 
  left_join(coverage_range) %>%
  left_join(indvs_covered) %>%
  mutate(TissueID = factor(TissueID, levels = tiss_order))

# get total samples available per tissue
total <- map_dbl(dots_per_indv_per_tiss, nrow)
total <- total[tiss_order]

# merge
Vg_coverage_tbl %<>% 
  mutate(`All Vg score available` = genes_w_Vg_per_tiss - `Median genes tested \nper indvididual`) %>%
  .[,-2] %>%
  gather(key, genes_w_Vg_per_tiss, -TissueID, -lower, - upper, -total)

############ Plot the results ###############
axis_labels <- tissue_site_detail[as.character(Vg_coverage_tbl$TissueID)]

Vg_coverage_tbl = Vg_coverage_tbl[complete.cases(Vg_coverage_tbl),]

Vg_coverage_plt <- 
  ggplot(Vg_coverage_tbl, aes(TissueID, genes_w_Vg_per_tiss, fill = key)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .3) +
  annotate("text", 
           y = text_positions, 
           x = 1:length(text_positions), 
           label = total, 
           angle = 270, 
           size = 3
           ) +
  ylab("Genes") + 
  xlab("Tissue") +
  
  # Add the arrow
  annotate("segment", 
           lineend = "butt", 
           linejoin = "mitre", 
           x = 25, y = 10000, xend = 19, yend = 8800, 
           arrow = arrow(length = unit(0.07, "inches"))) +
  annotate("text", label = "Total samples available", x = 25.5, y = 10100, size = 2.5, hjust = 0) +
  
  # Theme and aesthetic stuff
  gtex_v8_figure_theme() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1), 
    legend.position = c(.8, .9),
    legend.title = element_blank(), 
    plot.margin = margin(10,10,10,10)
  ) +
  scale_fill_brewer(
    palette = "Accent", 
    labels = c("Number of genes with a \nVg score available", "Median genes tested \nper individual")) +
  
  # Fix the axis labels, per Nicole's request
  scale_x_discrete(label = axis_labels)


########## B. Outlier Times ###########
outlier_times <- read_tsv("processed_input_data/figureS2/outlier_times_v8_Vg.tsv")

library("FSA")
frequent.outlier <- 
  outlier_times %$%
  binCI(times.gene.significant, times.gene.tested, type = "exact")[,1] %>% is_weakly_greater_than(.01)

outlier <- ifelse(frequent.outlier, "yes", "no") %>% as.factor
outlier_times_plt <- 
  ggplot(outlier_times, aes(times.gene.tested, times.gene.significant)) +
  geom_point(aes(colour = outlier)) +
  scale_color_manual(name = "Global Outlier", values = c("black", "red")) +
  xlab("Times gene is tested") +
  ylab("Times gene is significant") +
  theme(legend.position = c(0.2, 0.9)) +
  gtex_v8_figure_theme()

#outlier_times_plt


########## Panel C. Total Median DOTs per individual plot #############
# Load data from Rstudio server to show coding genes and highly expressed coding genes coverage distribution
totals_facet <- read_rds("processed_input_data/figureS2/tests_per_indv_faceted_plt_v8_Vg.rds")
totals_facet$method <- relevel(totals_facet$method, "TE")
#levels(totals_facet$method) <- c("TE", "ASE", "AS")

class_colors <- c(
  rgb(121, 92, 129, maxColorValue = 255), 
  rgb(23, 50, 75, maxColorValue = 255), 
  rgb(193, 205, 222, maxColorValue = 255)
)
names(class_colors) <- c("ASE", "AS", "TE")

# Save for plotting
total_med_test_across_methods_facet <- 
  ggplot(totals_facet, aes(x = total)) +
  geom_histogram(aes(fill = method), bins = 70) +
  scale_fill_manual(values = class_colors, labels = c("eOutliers", "aseOutliers", "sOutliers")) +
  xlim(0, NA) +
  xlab("Number of tested genes per individual") + ylab("Individuals") + 
  facet_wrap("filter", nrow = 2, scales = "free") +
  gtex_v8_figure_theme() +
  theme(legend.position = c(.42, .9), 
        legend.title = element_blank(),legend.direction = 'horizontal', 
        strip.background = element_rect(fill="white"))


#### Combine everything together with Cowplot ####
library(cowplot)
top_row <- plot_grid(Vg_coverage_plt, labels = 'A')
bottom_row <- plot_grid(outlier_times_plt, total_med_test_across_methods_facet, labels = c('B', 'C'))

png("generated_figures/fig_S2.png", width = 7, height = 7, units = "in", res = 300)
plot_grid(top_row, bottom_row, nrow = 2, rel_heights = c(1.6, 1))
dev.off()

