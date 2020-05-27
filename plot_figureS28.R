library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)
theme_set(theme_cowplot())


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

generate_figureS28_panel <- function(df, title, use_legend=TRUE, legend_height=-.38) {
    df$prediction_type = factor(df$prediction_type, levels=c("Watershed", "RIVER", "GAM"))
    outlier_type <- "ase"
    outlier_name <- "ASE"
    plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision", group="", linetype="", colour="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "steelblue3", "firebrick4")) +
                scale_linetype_manual(values=c("solid", "dotted", "solid")) +
                gtex_v8_figure_theme() + 
                draw_label(outlier_name,x=.5,y=.95,size=8)

    outlier_type <- "splicing"
    outlier_name <- "Splicing"
    plotter_splice <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision",group="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "steelblue3", "firebrick4")) +
                scale_linetype_manual(values=c("solid", "dotted", "solid")) +
                gtex_v8_figure_theme() +
                 draw_label(outlier_name,x=.5,y=.95,size=8)

    outlier_type <- "total_expression"
    outlier_name <- "Expression"
    plotter_te <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision", group="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "steelblue3", "firebrick4")) +
                scale_linetype_manual(values=c("solid", "dotted", "solid")) +
                gtex_v8_figure_theme() + 
                draw_label(outlier_name,x=.5,y=.95,size=8)

    legend <- get_legend(plotter_ase + theme(legend.position="bottom"))
    

    if (use_legend==TRUE) {
        combined_plots <- plot_grid(plotter_te+ theme(legend.position="none"), plotter_ase + theme(legend.position="none"), plotter_splice+ theme(legend.position="none"), rel_widths=c(1,1,1), nrow=1)
        combined <- ggdraw() + draw_text(title,size=8, x=.53,y=.92) + draw_plot(combined_plots,0,0.03,1,.83) + draw_plot(legend,.38,legend_height,1,1)
    } else {
        combined_plots <- plot_grid(plotter_te + labs(x="") + theme(legend.position="none"), plotter_ase + labs(x="") + theme(legend.position="none"), plotter_splice + labs(x="") + theme(legend.position="none"), rel_widths=c(1,1,1), nrow=1)
        combined <- ggdraw() + draw_text(title,size=8, x=.53,y=.92) + draw_plot(combined_plots,0,-.13,1,1.03)
    }

    return(combined)
}


##############
# Input files
figS28a_input_file = "processed_input_data/figureS28/figS28a_input_data.txt"
figS28b_input_file = "processed_input_data/figureS28/figS28b_input_data.txt"
figS28c_input_file = "processed_input_data/figureS28/figS28c_input_data.txt"
figS28d_input_file = "processed_input_data/figureS28/figS28d_input_data.txt"
figS28e_input_file = "processed_input_data/figureS28/figS28e_input_data.txt"
figS28f_input_file = "processed_input_data/figureS28/figS28f_input_data.txt"



##############
# Load in input data
figS28a_df <- read.table(figS28a_input_file, header=TRUE, sep="\t")
figS28b_df <- read.table(figS28b_input_file, header=TRUE, sep="\t")
figS28c_df <- read.table(figS28c_input_file, header=TRUE, sep="\t")
figS28d_df <- read.table(figS28d_input_file, header=TRUE, sep="\t")
figS28e_df <- read.table(figS28e_input_file, header=TRUE, sep="\t")
figS28f_df <- read.table(figS28f_input_file, header=TRUE, sep="\t")


###############
# Titles of each panel
title_a <- "Train models on genes where each of 3 outlier signals has at least one outlier individual (median p-value < 0.05)"
title_b <- "Train models on genes where each of 3 outlier signals has at least one outlier individual (median p-value < 0.1)"
title_c <- "Train models on genes where at least 1 outlier signal has at least one outlier individual (median p-value < 0.01)"
title_d <- "Evaluate models on genes where each of 3 outlier signals has at least one outlier individual (median p-value < 0.05)"
title_e <- "Evaluate models on genes where each of 3 outlier signals has at least one outlier individual (median p-value < 0.1)"
title_f <- "Evaluate models on genes where at least 1 outlier signal has at least one outlier individual (median p-value < 0.01)"



################
# generate figure S28 panels
figS28a <- generate_figureS28_panel(figS28a_df, title_a, use_legend=FALSE)
figS28b <- generate_figureS28_panel(figS28b_df, title_b, use_legend=FALSE)
figS28c <- generate_figureS28_panel(figS28c_df, title_c, use_legend=FALSE)
figS28d <- generate_figureS28_panel(figS28d_df, title_d, use_legend=FALSE)
figS28e <- generate_figureS28_panel(figS28e_df, title_e, use_legend=FALSE)
figS28f <- generate_figureS28_panel(figS28f_df, title_f, legend_height=-.45)

################
# Merge panels together with cowplot
figS28 <- plot_grid(figS28a, figS28b, figS28c, figS28d, figS28e, figS28f, ncol=1, rel_heights=c(1,1,1,1,1,1.28), labels=c("A","B","C", "D", "E", "F"))


# Save to output file
output_file <- "generated_figures/figureS28.pdf"
ggsave(figS28, file=output_file, width=7.2, height=9.0, units="in")
