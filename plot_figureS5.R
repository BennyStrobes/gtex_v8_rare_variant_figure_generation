library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

# Scatter plot showing corrected and uncorrected splicing outlier -log10(pvalues) at gene level colored by number of clusters mapped to the gene
corrected_vs_uncorrected_gene_level_pvalue_scatter_plot <- function(df, tissue_type) {
	df$corrected_pvalue <- -log10(df$corrected_pvalue + .000001)
	df$uncorrected_pvalue <- -log10(df$uncorrected_pvalue + .000001)

	df <- df[sample(nrow(df), 100000),]

	scatter <- ggplot(df, aes(x=uncorrected_pvalue, y=corrected_pvalue, colour=number_of_clusters)) + geom_point(size=.1) + 
				labs(x = "-log10(uncorrected pvalue)", y = "-log10(corrected pvalue)", colour="Number clusters") +
				gtex_v8_figure_theme()
	return(scatter)
}

# histogram showing corrected and uncorrected splicing outlier pvalues at gene level
corrected_vs_uncorrected_gene_level_pvalue_histogram <- function(df, tissue_type) {
	pvalues <- c()
	version <- c()

	pvalues <- c(pvalues, df$corrected_pvalue)
	version <- c(version, rep("corrected", length(df$corrected_pvalue)))
	pvalues <- c(pvalues, df$uncorrected_pvalue)
	version <- c(version, rep("uncorrected", length(df$uncorrected_pvalue)))

	df2 <- data.frame(pvalue=pvalues, version=factor(version))

	histo <- ggplot(df2, aes(x=pvalue, fill=version)) + geom_histogram(alpha=.35,position="identity",breaks = seq(0,1,.01)) + 
				labs(x = "pvalue", colour="", fill="") + 
				gtex_v8_figure_theme()
	return(histo)
}


##############
# Input data
splicing_outlier_pvalue_input_file = "processed_input_data/figureS5/figS5_input_data.txt"

##############
# Load in data
splicing_outlier_df <- read.table(splicing_outlier_pvalue_input_file, header=TRUE)


# Analysis restricted to the following tissue type
tissue_type="Muscle_Skeletal"


# Scatter plot showing corrected and uncorrected splicing outlier -log10(pvalues) at gene level colored by number of clusters mapped to the gene
cluster_correction_scatter <- corrected_vs_uncorrected_gene_level_pvalue_scatter_plot(splicing_outlier_df, tissue_type)

# histogram showing corrected and uncorrected splicing outlier pvalues at gene level
clustter_correction_histogram <- corrected_vs_uncorrected_gene_level_pvalue_histogram(splicing_outlier_df, tissue_type)

# Merge figure panels together with COWPLOT
figS5 <- plot_grid(cluster_correction_scatter, clustter_correction_histogram, ncol=1, labels=c("A","B"))

# Save to output file
output_file <- "generated_figures/figureS5.pdf"
ggsave(figS5, file=output_file, width=7.2, height=5, units="in")
