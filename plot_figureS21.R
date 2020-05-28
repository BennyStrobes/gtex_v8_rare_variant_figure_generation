# This script produces the heatmaps used to show sharing of expression 
# outlier status across tissues
# 
# Jonah Einson
# jeinson@nygenome.org

library(ggplot2)
library(ggdendro)
library(tidyverse)
library(magrittr)
library(cowplot)
library(RColorBrewer)

# Attach color information
gtex_key <- read_tsv("processed_input_data/gtex_colors.txt")

attach(gtex_key)
names(tissue_abbrv) <- tissue_id
color_hex <- str_c("#", color_hex)
names(color_hex) <- tissue_abbrv


# Load the results files
# _OA are the tables that have Only Available cross tissue gene comparisons
input_path = "processed_input_data/figureS21/heatmap_tables/"
ASE_sharing_OA <- read.csv(paste0(input_path, "ASE_sharing_heatmap_only_available_v8_Vg_no_go_indv_gene.tsv"), row.names = 1, check.names = F, sep = "\t")
TE_sharing_OA <- read.csv(paste0(input_path, "TE_sharing_heatmap_only_available.tsv"), row.names = 1, check.names = F, sep = "\t")
AS_sharing_OA <- read.csv(paste0(input_path, "AS_sharing_heatmap_only_available.tsv"), row.names = 1, check.names = F, sep = "\t")
# _NA files include NAs as a non-shared outlier status
ASE_sharing_NA <- read.csv(paste0(input_path, "ASE_sharing_heatmap_include_NA_v8_Vg_no_go_indv_gene.tsv"), row.names = 1, check.names = F, sep = "\t")
TE_sharing_NA <- read.csv(paste0(input_path, "TE_sharing_heatmap_include_NA.tsv"), row.names = 1, check.names = F, sep = "\t")
AS_sharing_NA <- read.csv(paste0(input_path, "AS_sharing_heatmap_include_NA.tsv"), row.names = 1, check.names = F, sep = "\t")

# Remove rows and columns with missing data
remove_missing <- function(table){
  to.keep <-  names(which(rowSums(is.na(table)) != nrow(table)))
  table[to.keep, to.keep]
}

ASE_sharing_OA %<>% remove_missing
TE_sharing_OA %<>% remove_missing
AS_sharing_OA %<>% remove_missing

ASE_sharing_NA %<>% remove_missing
TE_sharing_NA %<>% remove_missing
AS_sharing_NA %<>% remove_missing

# theme information 
gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

hmcol <- colorRampPalette(c("white","#d62222"))(100)

# Make the clustering order based on the ASE_sharing_OA
dat <- ASE_sharing_OA
dat[dat == 1] <- NA
dat[dat > .2] <- .2
distances_ASE_OA <- dist(dat, method = "canberra")
clustering_ASE_OA <- hclust(distances_ASE_OA)
row_order_ASE_OA <- colnames(dat)[clustering_ASE_OA$order]



auto_heatmap <- function(dat, exclude_legend = F, extract_legend = F){
  dat[dat == 1] <- NA
  dat[dat > .2] <- .2
  distances <- distances_ASE_OA
  clustering <- clustering_ASE_OA
  
  # Order the rows and columsn according to the hierarchical clustering
  row_order <- row_order_ASE_OA
  
  dat_df <- dat %>%
    rownames_to_column("Test") %>%
    gather("Discovery", "val", -"Test")
  
  dat_df$Test <- factor(dat_df$Test, levels = rev(row_order), ordered = T)
  dat_df$Discovery <- factor(dat_df$Discovery, levels = row_order, ordered = T)
  axis_colors <- color_hex[row_order]
  
  heatmap_plt <- 
    ggplot(dat_df, aes(y = Test, x = Discovery, fill = val)) +
    geom_tile(color = "grey") + 
    scale_fill_gradient("Sharing \nfraction", 
                        low = "white", high = "#d62222", 
                        na.value = "black", limits = c(0, .201)) + 
    scale_x_discrete(labels = rep("•", ncol(dat))) +
    scale_y_discrete(labels = rep("•", nrow(dat))) +
    gtex_v8_figure_theme() +
    theme(axis.text.x = element_text(color = axis_colors, size = 15, vjust = 1), 
          axis.text.y = element_text(color = rev(axis_colors), size = 15), 
          plot.margin = unit(c(.75 ,.25, .25 ,.25), "cm"))
  
  if(exclude_legend){
    heatmap_plt <- heatmap_plt + theme(legend.position = "none")
  }

  if(extract_legend){
    return(cowplot::get_legend(heatmap_plt))
  }
  # Combine together
  return(heatmap_plt)
}

############## Plot the Main Heatmaps #######################
# The only available heatmaps
top_margin <- unit(c(.5, .1, .1, .1), units = "cm")
OA_heatmaps <- 
  plot_grid(
    auto_heatmap(TE_sharing_OA, exclude_legend = T) + theme(plot.margin = top_margin) + xlab(""),
    auto_heatmap(ASE_sharing_OA, exclude_legend = T) + theme(plot.margin = top_margin) + xlab("") + ylab(""),
    auto_heatmap(AS_sharing_OA, exclude_legend = T) + theme(plot.margin = top_margin) + xlab("") + ylab(""),
    nrow = 1, labels = c("eOutliers", "aseOutliers", "sOutliers"), 
    label_x = .5, label_size = 8, hjust = c(0, 0, .25), vjust = 1.5
)


# The include NA heamaps
bottom_margin <- unit(c(.1, .1, .1, .1), units = "cm")
NA_heatmaps <- 
  plot_grid(
  auto_heatmap(TE_sharing_NA, exclude_legend = T) + theme(plot.margin = bottom_margin)+ xlab(""),
  auto_heatmap(ASE_sharing_NA, exclude_legend = T) + theme(plot.margin = bottom_margin) + ylab(""),
  auto_heatmap(AS_sharing_NA, exclude_legend = T) + theme(plot.margin = bottom_margin) + ylab("") + xlab(""),
  nrow = 1
)

# Top half of the plot #
top_half <- plot_grid(OA_heatmaps, NA_heatmaps, nrow = 2, align = "v", labels = "A")
legend <- auto_heatmap(ASE_sharing_OA, extract_legend = T)

# Add a quick label plot
plot_guide <- 
  ggplot() +
  annotate("text", x = 0, y = 1, angle = 90, label = "Yes") +
  annotate("text", x = 0, y = -1, angle = 90, label = "No") +
  annotate("text", x = -1, y = 0, angle = 90, size = 4, label = "Limit to co-available genes") +
  ylim(-2, 2) + xlim(-2, 1) +
  theme_void()

# plot_guide

top_half <- plot_grid(plot_guide, top_half, legend, rel_widths = c(.05, 1, .1), nrow = 1)



############### Panel B: NA Decrease Plot ######################

# How much does the percentage decrease by counting NA genes as non-shared?
mean_perc <- function(x){
  tmp <- x
  tmp <- tmp[tmp < 1]
  median(tmp, na.rm = T)
}

median_boot_ci <- function(x, R = 999){ # A bootstrap function that actually works!
  x <- x[x < 1]
  x <- x[!is.na(x)]
  
  low <- floor(R * .05)
  high <- ceiling(R * .95)
  
  out <- rep(0, R)
  for(i in 1:length(out)){
    out[i] <- median(sample(x, length(x), replace = T))
  }
  out <- sort(out)
  out[c(low, high)]
  
}

NA_decrease_tbl <- tibble(
  perc = c(
    mean_perc(TE_sharing_OA),
    mean_perc(TE_sharing_NA),
    
    mean_perc(ASE_sharing_OA),
    mean_perc(ASE_sharing_NA),
    
    mean_perc(AS_sharing_OA),
    mean_perc(AS_sharing_NA)
  ), 
  method = factor(c("eOutliers", "eOutliers", "aseOutliers", "aseOutliers", "sOutliers", "sOutliers"), levels = c("eOutliers", "aseOutliers", "sOutliers")),
  incl = factor(c("No", "Yes", "No", "Yes", "No", "Yes"), levels = c("No", "Yes"), ordered = T)
)

# This information from this part will go in the main text of the paper!
NA_decrease_tbl %>% filter(incl == "No")

intervals <- do.call(
  "rbind", 
  list(
    median_boot_ci(TE_sharing_OA), 
    median_boot_ci(TE_sharing_NA), 
    
    median_boot_ci(ASE_sharing_OA), 
    median_boot_ci(ASE_sharing_NA), 
    
    median_boot_ci(AS_sharing_OA), 
    median_boot_ci(AS_sharing_NA))
)
intervals <- as.data.frame(intervals)
colnames(intervals) <- c("lower", "upper")

NA_decrease_tbl %<>% cbind(intervals)

NA_decrease_plt <- 
  ggplot(NA_decrease_tbl, aes(method, perc, fill = incl)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  # scale_x_discrete(breaks = c("ASE", "AS", "TE"), 
  #                    labels = c("ASE", "Splicing", "Expression")
  #                    )+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .2, 
                position = position_dodge(.9)) +
  theme(plot.margin = unit(c(0,0,0,.5), units = "cm")) +
  scale_fill_manual(values = c("#d8b365", "#5ab4ac")) +
  ylab("Median Percent Sharing \n Across All Tissues") + 
  gtex_v8_figure_theme()

#NA_decrease_plt

######################## Panel C: Condensed Heatmap ###################
# Make an aggregated sharing percentage plot using BH corected ANEVA-DOT p values
# Make a barplot with median across columns on top and median across rows on bottom

ASE_sharing_BH <- read.csv(paste0(input_path, "ASE_sharing_heatmap_BH_v8_Vg.tsv"), sep = "\t", row.names = 1)

tmp <- ASE_sharing_BH
tmp[tmp == 1] <- NA

to_remove = which(colSums(is.na(tmp)) == nrow(tmp))
tmp = tmp[-to_remove, -to_remove]

vec_to_df <- function(vector){
  df_colnames <- names(vector)
  tibble(names = df_colnames, values = vector)
}

# discovery tissue
# How much, on average, do outliers in X replicate in another tissue?
discovery_median <- vec_to_df(apply(tmp, 2, median, na.rm = T))
discovery_range  <- apply(tmp, 2, range, na.rm = T)

discovery_median$lower <- discovery_range[1,]
discovery_median$upper <- discovery_range[2,]

discovery_median %<>% arrange(desc(values))
discovery_median$names <- factor(discovery_median$names, levels = discovery_median$names)
fill_cols <- color_hex[levels(discovery_median$names)]

discovery_plt <- 
  ggplot(discovery_median, aes(x = names, y = values)) +
  geom_bar(stat = "identity", fill = fill_cols) +
  scale_fill_manual(color_hex[levels(discovery_median$names)]) +
  geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.1) +
  ylab("Median % replication \n of  outliers in \n test tissue") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(.5, .5, 0. ,.5), units = "cm")) +
  gtex_v8_figure_theme()

# Test tissue
# (the row medians for the ASE heatmap from panel A)
test_median <- vec_to_df(apply(tmp, 1, median, na.rm = T))
test_range  <- apply(tmp, 1, range, na.rm = T)

test_median$lower <- test_range[1,]
test_median$upper <- test_range[2,]

test_median$names <- factor(test_median$names, levels = discovery_median$names)
fill_cols <- color_hex[levels(test_median$names)]

names(tissue_site_detail) <- tissue_abbrv

test_plt <- 
  ggplot(test_median, aes(x = names, y = values)) +
  geom_bar(stat = "identity", fill = fill_cols) +
  scale_fill_manual(color_hex[levels(test_median$names)]) +
  geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.1) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 270, hjust = 0),
        plot.margin = unit(c(0,.5,.5,.5), units = "cm")) +
  scale_x_discrete(labels = tissue_site_detail[levels(test_median$names)]) +
  #scale_x_discrete(labels = levels(test_median$names)) +
  ylab("Median % replication \n of outlier status \n from discovery tissues") +
  scale_y_reverse() +
  gtex_v8_figure_theme()

#test_plt

#### Generate the Color Key (Panel D) ####
gtex_key_plt <- 
  gtex_key %>%
  mutate(tiss_order = 1:nrow(.)) %>%
  mutate(tissue_site_detail = factor(tissue_site_detail, levels = tissue_site_detail))

gtex_key_plt$color_hex <- color_hex
color_key <-
  ggplot(gtex_key_plt, aes(0, tiss_order)) +
  geom_point(size = 2.5, color = color_hex) +
  geom_text(label = tissue_site_detail, aes(x = .02, hjust = 0), size = 2.5) +
  scale_y_reverse() +
  xlim(-.025, .35) +
  gtex_v8_figure_theme() +
  theme_void()
  

#color_key

#### Bottom Half ####
condensed_heatmap <- plot_grid(discovery_plt, test_plt, nrow = 2, rel_heights = c(1.3,2), align = "v", axis = "l")
bottom_left <- plot_grid(NA_decrease_plt, condensed_heatmap, nrow = 2, align = "none", rel_heights = c(1, 3), labels = c("B", "C"))
bottom_half <- plot_grid(bottom_left, color_key, nrow = 1, rel_widths = c(1, .4), labels = c(NA, "D"))


#### Assemble the whole thing!
whole_fig <- plot_grid(top_half, bottom_half, nrow = 2, rel_heights = c(1, 1.35))

png("generated_figures/Fig_S21.png", width = 8.5, height = 11, units = "in", res = 300)
whole_fig
dev.off()
