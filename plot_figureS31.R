library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
theme_set(theme_cowplot())


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

generate_figure_S31a <- function(df) {
    # Set ordering of outliers
    df$outlier_class=factor(df$outlier_class, levels=c("Expression", "ASE", "Splicing"))
    # Make scatter-plot
    plotter <- ggplot(df, aes(x=exact_betas, y=approximate_betas, colour=outlier_class)) + geom_point() +
            geom_abline() + 
            theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)) +
            labs(x="Beta (exact)", y="Beta (approximate)",colour="") + 
            scale_color_manual(values=c("#BFCDE0", "#7F5A83", "#0D324D")) +
            gtex_v8_figure_theme()
    return(plotter)
}

generate_figure_S31b <- function(df) {
    # Set ordering of inference methdos
    df$inference=factor(df$inference, levels=c("Exact","Approximate"))
    # Set ordering of models
    df$prediction_type=factor(df$prediction_type, levels=c("Watershed","RIVER"))

    #####################
    # PR curve for ASE
    outlier_type <- "ase"
    plotter_ase <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type, linetype=inference)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", linetype="", title="ASE") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                gtex_v8_figure_theme()
                #theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))
    #####################
    # PR curve for Splicing
    outlier_type <- "splicing"
    plotter_splice <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type, linetype=inference)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", linetype="", title="Splicing") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                gtex_v8_figure_theme()
                #theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))
    #####################
    # PR curve for Expression
    outlier_type <- "total_expression"
    plotter_te <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, colour=prediction_type,linetype=inference)) + geom_line() + 
                labs(x="Recall", y="Precision", colour="", linetype="", title="Expression") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "black")) +
                gtex_v8_figure_theme()
                #theme(text = element_text(size=14),axis.text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=14))

    legend <- get_legend(plotter_ase)
    combined_plots <- plot_grid(plotter_te + theme(legend.position="none"), plotter_ase + theme(legend.position="none"),plotter_splice+ theme(legend.position="none"), nrow=1)


    return(plot_grid(combined_plots,legend,ncol=1,rel_heights=c(1,.06)))
}


##############
# Input files
figS31a_input_file = "processed_input_data/figureS31/figS31a_input_data.txt"
figS31b_input_file = "processed_input_data/figureS31/figS31b_input_data.txt"


##############
# Load in data
figS31a_df <- read.table(figS31a_input_file, header=TRUE, sep="\t")
figS31b_df <- read.table(figS31b_input_file, header=TRUE, sep="\t")

################
# generate figure S31 panels
figS31a <- generate_figure_S31a(figS31a_df)
figS31b <- generate_figure_S31b(figS31b_df)


################
# Merge panels together with cowplot
figS31 <- plot_grid(figS31a, figS31b, ncol=1,labels=c("A","B"))

################
# Save to output file
output_file <- "generated_figures/figureS31.pdf"
ggsave(figS31, file=output_file, width=7.2, height=5.0, units="in")
