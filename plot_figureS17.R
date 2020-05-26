library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

make_variant_enrichment_within_splicing_outlier_error_bar_plot <- function(variant_enrichment_df) {
	# compute odds ratios
	orat = (variant_enrichment_df$num_outlier_clusters_with_rv/variant_enrichment_df$num_outlier_clusters)/(variant_enrichment_df$num_inlier_clusters_with_rv/variant_enrichment_df$num_inlier_clusters)

	# Compute error bounds on odds ratios
	log_bounds <- 1.96*sqrt((1.0/variant_enrichment_df$num_outlier_clusters_with_rv) - (1.0/variant_enrichment_df$num_outlier_clusters) + (1.0/variant_enrichment_df$num_inlier_clusters_with_rv) - (1.0/variant_enrichment_df$num_inlier_clusters))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)

	# Make data frame for plotting
	distance_ordering <- c("[0,2]", "[3,4]", "[5,6]", "[7,8]", "[9,10]", "[11,100]")
	df <- data.frame(pvalues=factor(variant_enrichment_df$pvalues), distance=factor(variant_enrichment_df$distance, levels=distance_ordering), odds_ratios=orat, lower_bounds=lower_bound, upper_bounds=upper_bound)
	# Plot
	dodge <- position_dodge(width=0.5)
	error_bar_plot <- ggplot() + geom_errorbar(data=df, mapping=aes(x=distance,ymin=lower_bounds, ymax=upper_bounds, colour=pvalues), position=dodge,width=0.0) +
					geom_point(data=df, mapping=aes(x=distance, y=odds_ratios, colour=pvalues), position=dodge) +
					labs(x = "Base pair window around splice site", y = "Relative risk", colour="p-value") +
					geom_hline(yintercept=1) + 
					gtex_v8_figure_theme() + theme(legend.position = c(0.8, 0.8))

	return(error_bar_plot)
}


make_junction_usage_concensus_read_count_boxplot_across_tissues <- function(junction_usage_df) {
	plotter <- ggplot(junction_usage_df, aes(x=type, y=odds_ratio, fill=type)) + geom_boxplot() + 
		labs(x="", y="Junction usage", fill="") + 
		scale_fill_manual(values=c("dodgerblue", "violetred1")) +
		theme(legend.position="none") + 
		geom_hline(yintercept=0) + 
		gtex_v8_figure_theme()
	return(plotter)

}

##############
# Input data for figure S17A
variant_enrichment_input_file = "processed_input_data/figureS17/figS17a_input_data.txt"
# Input data for figure S17B
junction_usage_input_file = "processed_input_data/figureS17/figS17b_input_data.txt"



##############
# Load in variant enrichment data
variant_enrichment_df <- read.table(variant_enrichment_input_file, header=TRUE, sep="\t")
# Load in junction usage data
junction_usage_df <- read.table(junction_usage_input_file, header=TRUE, sep="\t")



##############
# plot S17a
figS17a <- make_variant_enrichment_within_splicing_outlier_error_bar_plot(variant_enrichment_df)
# plot S17b
figS17b <- make_junction_usage_concensus_read_count_boxplot_across_tissues(junction_usage_df)


# Merge two panels together with cowplot
figS17 <- plot_grid(figS17a, figS17b, labels=c("A","B"), ncol=2, rel_widths=c(.6,.4))



# Save to output file
output_file <- "generated_figures/figureS17.pdf"
ggsave(figS17, file=output_file, width=7.2, height=3, units="in")
