library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
library(reshape)
theme_set(theme_cowplot())


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

generate_figureS33_panel <- function(df, outlier_type) {
    cols <- c( "c1" = "steelblue3", "c2" = "firebrick4" )
    plotter <- ggplot(df) +
               geom_segment(aes(x=tissues_position, xend=tissues_position, y=auc_watershed, yend=auc_river), color="grey") +
               geom_point( aes(x=tissues_position, y=auc_watershed, color="c1"), size=1.5) +
               geom_point( aes(x=tissues_position, y=auc_river, color="c2"), size=1.5) +
               scale_color_manual(name="", breaks=c("c1","c2"), values=cols, labels=c("tissue Watershed", "tissue RIVER")) +
               ggtitle(outlier_type) +
                xlab("") +
               ylab("AUC (PR)") + 
               theme(legend.position="bottom") +
               theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
               gtex_v8_figure_theme() +
               scale_x_continuous(breaks=1:length(df$tissue),labels=gsub("_"," ",df$tissue,fixed=TRUE)) 
    return(plotter)
}


##############
# Input files (one for splicing, ase, and expression)
figS33_expression_input_file = "processed_input_data/figureS33/figS33_expression_input_data.txt"
figS33_ase_input_file = "processed_input_data/figureS33/figS33_ase_input_data.txt"
figS33_splicing_input_file = "processed_input_data/figureS33/figS33_splicing_input_data.txt"


##############
# Load in data
expression_df <- read.table(figS33_expression_input_file, header=TRUE, sep="\t")
ase_df <- read.table(figS33_ase_input_file, header=TRUE, sep="\t")
splicing_df <- read.table(figS33_splicing_input_file, header=TRUE, sep="\t")




################
# generate figure S32 panels
figS33_expression <- generate_figureS33_panel(expression_df, "Expression")
figS33_ase <- generate_figureS33_panel(ase_df, "ASE")
figS33_splicing <- generate_figureS33_panel(splicing_df, "Splicing")



################
# Merge panels together with cowplot

legend <- get_legend(figS33_ase)
figS33_tmp <- plot_grid(figS33_expression + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), figS33_ase + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), figS33_splicing +theme(legend.position="none",plot.margin=grid::unit(c(0,0,0,0), "mm")), rel_heights=c(1,1,2), ncol=1)
figS33 <- ggdraw() + draw_plot(figS33_tmp,0,0,1,1) + draw_plot(legend,.7,.47,1,1)


################
# Save to output file
output_file <- "generated_figures/figureS33.pdf"
ggsave(figS33,file=output_file,width=7.2, height=7.2, units="in")


