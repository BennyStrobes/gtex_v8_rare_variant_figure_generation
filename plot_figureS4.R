args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


# Make bar plot showing number unique genes (with jxns mapped to it) in each tissue
barplot_showing_number_of_genes_per_tissue <- function(df) {
	# Order tissues alphabetically 
	df$tissue = factor(df$tissue, levels=rev(as.character(df$tissue)))
	# Make barplot
	bar_plot <- ggplot(data=df, aes(x=tissue,y=number_of_genes)) + geom_bar(stat="identity",fill="darkcyan") + coord_flip() +
		labs(x = "", y = "Number of genes") +
		gtex_v8_figure_theme()
	return(bar_plot)
}

# Make bar plot showing number unique LeafCutter Clusters (with jxns mapped to it) in each tissue
barplot_showing_number_of_clusters_per_tissue <- function(df) {
	# Order tissues alphabetically 
	df$tissue = factor(df$tissue, levels=rev(as.character(df$tissue)))
	# Make barplot
	bar_plot <- ggplot(data=df, aes(x=tissue,y=number_of_clusters)) + geom_bar(stat="identity",fill="mediumpurple") + coord_flip() +
				labs(x = "", y = "Number of clusters") +
				gtex_v8_figure_theme()
	return(bar_plot)
}


barplot_showing_number_of_jxns_per_tissue <- function(df) {
	# Order tissues alphabetically 
	df$tissue = factor(df$tissue, levels=rev(as.character(df$tissue)))
	# Make bar plot
	bar_plot <- ggplot(data=df, aes(x=tissue,y=number_of_junctions)) + geom_bar(stat="identity",fill="steelblue3") + coord_flip() +
				labs(x = "", y = "Number of junctions") +
				gtex_v8_figure_theme()
	return(bar_plot)

}


##############
# Input data
jxn_count_data_input_file = "processed_input_data/figureS4/figS4_input_data.txt"


##############
# Load in input data
jxn_count_df <- read.table(jxn_count_data_input_file, header=TRUE, sep="\t")


##############
# Make bar plot showing number unique genes (with junctions mapped to it) in each tissue
genes_per_tissue_bar_plot <- barplot_showing_number_of_genes_per_tissue(jxn_count_df)

##############
# Make bar plot showing number unique LeafCutter Clusters (with junctions mapped to it) in each tissue
clusters_per_tissue_bar_plot <- barplot_showing_number_of_clusters_per_tissue(jxn_count_df)

##############
# Make bar plot showing number unique exon-exon junctions in each tissue
jxns_per_tissue_bar_plot <- barplot_showing_number_of_jxns_per_tissue(jxn_count_df)

##############
# Use cowplot to combine panels
figS4 <- plot_grid(jxns_per_tissue_bar_plot, clusters_per_tissue_bar_plot + theme(axis.text.y=element_blank()), genes_per_tissue_bar_plot + theme(axis.text.y=element_blank()), labels=c("A","B","C"), ncol=3, rel_widths=c(2.0,1.0,1.0))

##############
# Save to output
ggsave(figS4, file=paste0("generated_figures/figureS4.pdf"),width=7.2, height=5, units="in")
