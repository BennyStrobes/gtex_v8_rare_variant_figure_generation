library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


make_figS6_panel <- function(x_axis_pvalues, y_axis_pvalues, color_vector, x_axis_label, y_axis_label) {
	# randomly sample indices for viz
	random_indices <- sample (c(1:length(x_axis_pvalues)), size=length(x_axis_pvalues)/100.0, replace=F)
	x_axis_pvalues = x_axis_pvalues[random_indices]
	y_axis_pvalues = y_axis_pvalues[random_indices]
	color_vector = color_vector[random_indices]

	df <- data.frame(stadard_pvalue=-log10(x_axis_pvalues+1e-6), alt_pvalue=-log10(y_axis_pvalues+1e-6), fraction=color_vector)
	p <- ggplot(df, aes(x=stadard_pvalue, y=alt_pvalue,colour=fraction)) + 
  			geom_point(size=.0001) +
  			gtex_v8_figure_theme() +
  			theme(legend.position="bottom") + 
  			labs(x = paste0("-log10(p-value) [ ", x_axis_label, " ]"), y=paste0("-log10(p-value) [ ", y_axis_label, " ]"),colour="Fraction of reads\nfrom one junction") 
  	return(p)
}


##############
# Input data
figS6_input_file = "processed_input_data/figureS6/figS6_input_data.txt"

##############
# Load in data
figS6_df <- read.table(figS6_input_file, header=TRUE)


figS6a <- make_figS6_panel(figS6_df$standard_prior_20K_reads, figS6_df$standard_prior_10K_reads, figS6_df$fraction, "20000 reads", "10000 reads")
figS6b <-make_figS6_panel(figS6_df$standard_prior_20K_reads, figS6_df$standard_prior_100K_reads, figS6_df$fraction, "20000 reads", "100000 reads")
figS6c <-make_figS6_panel(figS6_df$standard_prior_20K_reads, figS6_df$no_prior_20K_reads, figS6_df$fraction, "Standard prior", "No prior")

# Make combined plot
figS6 <- plot_grid(figS6a, figS6b, figS6c, ncol=2, labels=c("A","B", "C"))

# Save to output file
output_file <- "generated_figures/figureS6.pdf"
ggsave(figS6, file=output_file, width=7.2, height=7, units="in")
