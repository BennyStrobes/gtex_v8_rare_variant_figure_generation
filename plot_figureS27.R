library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)
theme_set(theme_cowplot())

gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

generate_figureS27 <- function(df) {
	# Specify ordering
	df$prediction_type = factor(df$prediction_type, levels=c("Watershed", "CADD"))
	# Make ASE PLOT
	outlier_type <- "ase"
	outlier_name <- "ASE"
  	plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision", group="", linetype="", colour="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                scale_linetype_manual(values=c("solid", "solid")) +
                gtex_v8_figure_theme() + 
                draw_label(outlier_name,x=.5,y=.95,size=8)
    # Make SPlicing plot
	outlier_type <- "splicing"
	outlier_name <- "Splicing"
  	plotter_splice <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision",group="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                scale_linetype_manual(values=c("solid", "solid")) +
                gtex_v8_figure_theme() +
                 draw_label(outlier_name,x=.5,y=.95,size=8)
    # Make expression plot
	outlier_type <- "total_expression"
	outlier_name <- "Expression"
  	plotter_te <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision", group="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                scale_linetype_manual(values=c("solid", "solid")) +
                gtex_v8_figure_theme() + 
                draw_label(outlier_name,x=.5,y=.95,size=8)

    legend <- get_legend(plotter_ase + theme(legend.position="bottom"))
    combined_plots <- plot_grid(plotter_te+ theme(legend.position="none"), plotter_ase + theme(legend.position="none"), plotter_splice+ theme(legend.position="none"), rel_widths=c(1,1,1), nrow=1)

	combined <- ggdraw() + draw_plot(combined_plots,0,.07,1,.9) + draw_plot(legend,.44,-0.42,1,1)

	return(combined)

}

##############
# Input data
figS27_input_file = "processed_input_data/figureS27/figS27_input_data.txt"



##############
# Load in input data
figS27_df <- read.table(figS27_input_file, header=TRUE, sep="\t")


################
# generate figure S27
figS27 <- generate_figureS27(figS27_df)

# Save to output file
output_file <- "generated_figures/figureS27.pdf"
ggsave(figS27, file=output_file, width=10, height=3.0, units="in")
