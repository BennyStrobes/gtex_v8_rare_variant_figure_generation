library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

purine_pyrimidine_variant_changes_in_specific_region <- function(outlier_distances, inlier_distances, ss_type, ppt_start, ppt_end) {
	inlier_distances_region <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance > ppt_start & inlier_distances$distance <= ppt_end & as.character(inlier_distances$annotated_splice_site)=="annotated",]
	outlier_distances_region <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance > ppt_start & outlier_distances$distance <= ppt_end & as.character(outlier_distances$annotated_splice_site)=="annotated",]

	num_outliers <- dim(outlier_distances_region)[1]
	num_inliers <- dim(inlier_distances_region)[1]

	odds_ratios <- c()
	lower_bounds <- c()
	upper_bounds <- c()
	variant_class <- c()

	# Purine -> pyrimidine
	pu_py_outlier <- sum((outlier_distances_region$major_allele == "A" | outlier_distances_region$major_allele == "G") & (outlier_distances_region$variant_allele == "C" | outlier_distances_region$variant_allele == "T"))
	pu_py_inlier <- sum((inlier_distances_region$major_allele == "A" | inlier_distances_region$major_allele == "G") & (inlier_distances_region$variant_allele == "C" | inlier_distances_region$variant_allele == "T"))
	orat <- (pu_py_outlier/num_outliers)/(pu_py_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/pu_py_outlier) - (1.0/num_outliers) + (1.0/pu_py_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Purine to Pyrimidine")

	# Purine -> purine
	pu_pu_outlier <- sum((outlier_distances_region$major_allele == "A" | outlier_distances_region$major_allele == "G") & (outlier_distances_region$variant_allele == "A" | outlier_distances_region$variant_allele == "G"))
	pu_pu_inlier <- sum((inlier_distances_region$major_allele == "A" | inlier_distances_region$major_allele == "G") & (inlier_distances_region$variant_allele == "A" | inlier_distances_region$variant_allele == "G"))
	orat <- (pu_pu_outlier/num_outliers)/(pu_pu_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/pu_pu_outlier) - (1.0/num_outliers) + (1.0/pu_pu_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Purine to Purine")

	# pyrimidine -> purine
	py_pu_outlier <- sum((outlier_distances_region$major_allele == "C" | outlier_distances_region$major_allele == "T") & (outlier_distances_region$variant_allele == "A" | outlier_distances_region$variant_allele == "G"))
	py_pu_inlier <- sum((inlier_distances_region$major_allele == "C" | inlier_distances_region$major_allele == "T") & (inlier_distances_region$variant_allele == "A" | inlier_distances_region$variant_allele == "G"))
	orat <- (py_pu_outlier/num_outliers)/(py_pu_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/py_pu_outlier) - (1.0/num_outliers) + (1.0/py_pu_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Pyrimidine to Purine")


	# pyrimidine -> pyrimidine
	py_py_outlier <- sum((outlier_distances_region$major_allele == "C" | outlier_distances_region$major_allele == "T") & (outlier_distances_region$variant_allele == "C" | outlier_distances_region$variant_allele == "T")) 
	py_py_inlier <- sum((inlier_distances_region$major_allele == "C" | inlier_distances_region$major_allele == "T") & (inlier_distances_region$variant_allele == "C" | inlier_distances_region$variant_allele == "T")) 
	orat <- (py_py_outlier/num_outliers)/(py_py_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/py_py_outlier) - (1.0/num_outliers) + (1.0/py_py_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Pyrimidine to Pyrimidine")

	df <- data.frame(odds_ratio=odds_ratios, dist=1:length(variant_class), lower_bound=lower_bounds, upper_bound=upper_bounds, variant_class=factor(variant_class))



	error_bar_plot_annotated <-  ggplot() + geom_errorbar(data=df, mapping=aes(x=dist,ymin=lower_bound, ymax=upper_bound),color="darkorchid", width=0) +
					geom_point(data=df, mapping=aes(x=dist, y=odds_ratio), color="darkorchid") +
					labs(x = "", y = "Relative risk", title="Annotated splice site") +
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					gtex_v8_figure_theme() + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
					scale_x_continuous(breaks=1:length(df$dist), labels=df$variant_class)



	inlier_distances_region <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance > ppt_start & inlier_distances$distance <= ppt_end & as.character(inlier_distances$annotated_splice_site)=="novel",]
	outlier_distances_region <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance > ppt_start & outlier_distances$distance <= ppt_end & as.character(outlier_distances$annotated_splice_site)=="novel",]

	num_outliers <- dim(outlier_distances_region)[1]
	num_inliers <- dim(inlier_distances_region)[1]

	odds_ratios <- c()
	lower_bounds <- c()
	upper_bounds <- c()
	variant_class <- c()

	# Purine -> pyrimidine
	pu_py_outlier <- sum((outlier_distances_region$major_allele == "A" | outlier_distances_region$major_allele == "G") & (outlier_distances_region$variant_allele == "C" | outlier_distances_region$variant_allele == "T"))
	pu_py_inlier <- sum((inlier_distances_region$major_allele == "A" | inlier_distances_region$major_allele == "G") & (inlier_distances_region$variant_allele == "C" | inlier_distances_region$variant_allele == "T"))
	orat <- (pu_py_outlier/num_outliers)/(pu_py_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/pu_py_outlier) - (1.0/num_outliers) + (1.0/pu_py_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Purine to Pyrimidine")

	# Purine -> purine
	pu_pu_outlier <- sum((outlier_distances_region$major_allele == "A" | outlier_distances_region$major_allele == "G") & (outlier_distances_region$variant_allele == "A" | outlier_distances_region$variant_allele == "G"))
	pu_pu_inlier <- sum((inlier_distances_region$major_allele == "A" | inlier_distances_region$major_allele == "G") & (inlier_distances_region$variant_allele == "A" | inlier_distances_region$variant_allele == "G"))
	orat <- (pu_pu_outlier/num_outliers)/(pu_pu_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/pu_pu_outlier) - (1.0/num_outliers) + (1.0/pu_pu_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Purine to Purine")

	# pyrimidine -> purine
	py_pu_outlier <- sum((outlier_distances_region$major_allele == "C" | outlier_distances_region$major_allele == "T") & (outlier_distances_region$variant_allele == "A" | outlier_distances_region$variant_allele == "G"))
	py_pu_inlier <- sum((inlier_distances_region$major_allele == "C" | inlier_distances_region$major_allele == "T") & (inlier_distances_region$variant_allele == "A" | inlier_distances_region$variant_allele == "G"))
	orat <- (py_pu_outlier/num_outliers)/(py_pu_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/py_pu_outlier) - (1.0/num_outliers) + (1.0/py_pu_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Pyrimidine to Purine")


	# pyrimidine -> pyrimidine
	py_py_outlier <- sum((outlier_distances_region$major_allele == "C" | outlier_distances_region$major_allele == "T") & (outlier_distances_region$variant_allele == "C" | outlier_distances_region$variant_allele == "T")) 
	py_py_inlier <- sum((inlier_distances_region$major_allele == "C" | inlier_distances_region$major_allele == "T") & (inlier_distances_region$variant_allele == "C" | inlier_distances_region$variant_allele == "T")) 
	orat <- (py_py_outlier/num_outliers)/(py_py_inlier/num_inliers)
	log_bounds <- 1.96*sqrt((1.0/py_py_outlier) - (1.0/num_outliers) + (1.0/py_py_inlier) - (1.0/num_inliers))
	upper_bound <- orat*exp(log_bounds)
	lower_bound <- orat*exp(-log_bounds)
	odds_ratios <- c(odds_ratios, orat)
	lower_bounds <- c(lower_bounds, lower_bound)
	upper_bounds <- c(upper_bounds, upper_bound)
	variant_class <- c(variant_class, "Pyrimidine to Pyrimidine")

	df <- data.frame(odds_ratio=odds_ratios, dist=1:length(variant_class), lower_bound=lower_bounds, upper_bound=upper_bounds, variant_class=factor(variant_class))



	error_bar_plot_novel <-  ggplot() + geom_errorbar(data=df, mapping=aes(x=dist,ymin=lower_bound, ymax=upper_bound),color="darkorchid", width=0.0) +
					geom_point(data=df, mapping=aes(x=dist, y=odds_ratio), color="darkorchid") +
					labs(x = "", y = "Relative risk", title="Novel splice site") +
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					gtex_v8_figure_theme() + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
					scale_x_continuous(breaks=1:length(df$dist), labels=df$variant_class)

	error_bar_plot <- plot_grid(error_bar_plot_annotated, error_bar_plot_novel, ncol=2, labels=c("A","B"))

	return(error_bar_plot)

}

##############
# Input data
outlier_variant_input_file = "processed_input_data/figureS20/figS20_outliers_input_data.txt"
inlier_variant_input_file = "processed_input_data/figureS20/figS20_inliers_input_data.txt"

##############
# Load in data
outlier_variant_df <- read.table(outlier_variant_input_file, header=TRUE)
inlier_variant_df <- read.table(inlier_variant_input_file, header=TRUE)


##############
# make plot
ppt_start <- -35
ppt_end <- -5
ss_type <- "acceptor"

figS20 <- purine_pyrimidine_variant_changes_in_specific_region(outlier_variant_df, inlier_variant_df, ss_type, ppt_start, ppt_end)





###################
# Save to output file
output_file <- "generated_figures/figureS20.pdf"
ggsave(figS20, file=output_file, width=7.2, height=3.2, units="in")
