library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)

get_color_vector <- function(tissue_colors, tissue_names) {
	colors <- c()
	for (iter in 1:length(tissue_names)) {
		tissue_name <- as.character(tissue_names[iter])
		if (tissue_name == "Brain_Spinal_cord_cervical_c-1") {
			tissue_name = "Brain_Spinal_cord_cervical_c1"
		}
		if (tissue_name == "Cells_EBV-transformed_lymphocytes") {
			tissue_name = "Cells_EBVtransformed_lymphocytes"
		}

		indices = tissue_name==tissue_colors$tissue_id

		color <- tissue_colors$tissue_color_hex[indices]
		
		colors <- c(colors, paste0('#',color))
	}
	return(colors)
}


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

tbt_variant_outlier_enrichment_errorbar_plot <- function(enrichments, color_vector) {
	num_tissues <- dim(enrichments)[1]

	# Initialize vectors
	tissue_names <- c()
	odds_ratios <- c()
	lower_bounds <- c()
	upper_bounds <- c()

	# Loop through tissues
	for (tissue_number in 1:num_tissues) {
		# Extract tissue name for this line
		tissue_name <- as.character(enrichments[tissue_number, 1])
		# Compute odds ratios for this line
		a <- enrichments[tissue_number, 2]
		b <- enrichments[tissue_number, 3]
		c <- enrichments[tissue_number, 4]
		d <- enrichments[tissue_number, 5]
		orat <- (a/b)/(c/d)
		# Compute error bars for this orat
		log_orat <- log(orat)
		log_bounds <- 1.96*sqrt((1.0/a) - (1.0/b) + (1.0/c) - (1.0/d))
		#upper_bound <- exp(log_orat + log_bounds)
		#lower_bound <- exp(log_orat - log_bounds) 
		upper_bound <- orat*exp(log_bounds)
		lower_bound <- orat*exp(-log_bounds)

		# Add information to vectors
		tissue_names <- c(tissue_names, tissue_name)
		odds_ratios <- c(odds_ratios, orat)
		lower_bounds <- c(lower_bounds, lower_bound)
		upper_bounds <- c(upper_bounds, upper_bound)
	}

	tissue_names <- gsub("_", " ", tissue_names)
	# Add information to data frame
	df <- data.frame(tissue_names=factor(tissue_names), odds_ratios=odds_ratios, lower_bounds=lower_bounds, upper_bounds=upper_bounds)
	error_bar_plot <- ggplot() + geom_errorbar(data=df, mapping=aes(x=tissue_names,ymin=lower_bounds, ymax=upper_bounds), width=0, colour=color_vector) +
					geom_point(data=df, mapping=aes(x=tissue_names, y=odds_ratios), colour=color_vector) +
					labs(x = "Tissue", y = "Relative risk") +
					geom_hline(yintercept=1) + 
					gtex_v8_figure_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))

	return(error_bar_plot)
}

odds_ratio_concensus_boxplot_across_tissues <- function(df) {
	plotter <- ggplot(df, aes(x=type, y=odds_ratio, fill=type)) + geom_boxplot() + 
		labs(x="", y="Junction usage", fill="") + 
		scale_fill_manual(values=c("dodgerblue", "violetred1")) +
		theme(legend.position="none") + 
		geom_hline(yintercept=0) + 
		gtex_v8_figure_theme()
	return(plotter)
}

odds_ratio_ppt_boxplot_across_tissues <- function(df) {
	plotter <- ggplot(df, aes(x=type, y=odds_ratio, fill=type)) + geom_boxplot() + 
		theme(text = element_text(size=12),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=12), legend.title = element_text(size=11)) +
		labs(x="", y="Junction usage", fill="") + 
		scale_fill_manual(values=c("dodgerblue", "violetred1")) +
		theme(legend.position="none") + 
		geom_hline(yintercept=0) +
		gtex_v8_figure_theme()
	return(plotter)
}


##############
# Input data
tissue_names_file = "processed_input_data/figureS24/figS24_tissue_names_input_data.txt"
tissue_colors_file = "processed_input_data/figureS24/figS24_tissue_colors_input_data.txt"
tbt_variant_outlier_overlap_file = "processed_input_data/figureS24/figS24_tissue_by_tissue_variant_sOutlier_overlap_input_data.txt"
tbt_concensus_jxn_usage_file = "processed_input_data/figureS24/figS24_tissue_by_tissue_concensus_jxn_usage_input_data.txt"
tbt_ppt_jxn_usage_file = "processed_input_data/figureS24/figS24_tissue_by_tissue_ppt_jxn_usage_input_data.txt"


########################
# Extract tissue color information
########################
# Extract vector tissue names
tissue_names <- as.character(unlist(read.table(tissue_names_file,header=FALSE), use.names=FALSE))
tissue_colors = read.table(tissue_colors_file, header = T, stringsAsFactors = F, sep = "\t")
# Get vector of hex colors in correct order
color_vector <- get_color_vector(tissue_colors, tissue_names)

##############
# Load in data
tbt_variant_outlier_df <- read.table(tbt_variant_outlier_overlap_file, header=TRUE, sep="\t")
tbt_concensus_jxn_usage_df <- read.table(tbt_concensus_jxn_usage_file, header=TRUE, sep="\t")
tbt_ppt_jxn_usage_df <- read.table(tbt_ppt_jxn_usage_file, header=TRUE, sep="\t")



##############
# Make panel A
figS24a <- tbt_variant_outlier_enrichment_errorbar_plot(tbt_variant_outlier_df, color_vector)

##############
# Make panel B
figS24b <- odds_ratio_concensus_boxplot_across_tissues(tbt_concensus_jxn_usage_df)

##############
# Make panel C
figS24c <- odds_ratio_ppt_boxplot_across_tissues(tbt_ppt_jxn_usage_df)

##############
# Combine panels with cowplot
figS24 <- plot_grid(figS24a, plot_grid(figS24b, figS24c, labels=c("B","C"), ncol=2, rel_widths=c(.5,.5)), ncol=1, rel_heights=c(.6,.44), labels=c("A",""))


# Save to output file
output_file <- "generated_figures/figureS24.pdf"
ggsave(figS24, file=output_file, width=7.2, height=6.0, units="in")
