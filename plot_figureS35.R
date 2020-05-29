library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
library(reshape)
theme_set(theme_cowplot())


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

generate_figureS35_panel <- function(df, outlier_type) {
  #df <- data.frame(auc_watershed = auc_vi, auc_river=auc_exact, tissue=ordered_tissue_names, tissues_position=1:length(ordered_tissue_names))
    cols <- c( "c1" = "steelblue3", "c2" = rgb(0.2,0.7,0.1,0.5) )
    plotter <- ggplot(df) +
           geom_segment(aes(x=tissues_position, xend=tissues_position, y=auc_watershed, yend=auc_river), color="grey") +
           geom_point( aes(x=tissues_position, y=auc_watershed, color="c1"), size=1.5) +
               geom_point( aes(x=tissues_position, y=auc_river, color="c2"), size=1.5) +
               scale_color_manual(name="", breaks=c("c1","c2"), values=cols, labels=c("tissue Watershed", "River")) +
               ggtitle(outlier_type) +
                xlab("") +
               ylab("AUC (PR)") + 
               theme(legend.position="bottom") +
               theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
               gtex_v8_figure_theme() +
               #scale_y_continuous(expand = c(0, 0), limits = c(0, 1.06)) + 
               scale_x_continuous(breaks=1:length(df$tissue),labels=gsub("_"," ",df$tissue,fixed=TRUE)) 
    return(plotter)
}


##############
# Input files (one for splicing, ase, and expression)
figS35_expression_input_file = "processed_input_data/figureS35/figS35_expression_input_data.txt"
figS35_ase_input_file = "processed_input_data/figureS35/figS35_ase_input_data.txt"
figS35_splicing_input_file = "processed_input_data/figureS35/figS35_splicing_input_data.txt"


##############
# Load in data
expression_df <- read.table(figS35_expression_input_file, header=TRUE, sep="\t")
ase_df <- read.table(figS35_ase_input_file, header=TRUE, sep="\t")
splicing_df <- read.table(figS35_splicing_input_file, header=TRUE, sep="\t")




################
# generate figure S32 panels
figS35_expression <- generate_figureS35_panel(expression_df, "Expression")
figS35_ase <- generate_figureS35_panel(ase_df, "ASE")
figS35_splicing <- generate_figureS35_panel(splicing_df, "Splicing")



################
# Merge panels together with cowplot

legend <- get_legend(figS35_ase)
figS35_tmp <- plot_grid(figS35_expression + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), figS35_ase + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), figS35_splicing +theme(legend.position="none",plot.margin=grid::unit(c(0,0,0,0), "mm")), rel_heights=c(1,1,2), ncol=1)
figS35 <- ggdraw() + draw_plot(figS35_tmp,0,0,1,1) + draw_plot(legend,.7,.47,1,1)
################
# Save to output file
output_file <- "generated_figures/figureS35.pdf"
ggsave(figS35,file=output_file,width=7.2, height=7.2, units="in")


