library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
library(reshape)
theme_set(theme_cowplot())


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

generate_figureS36_panel <- function(df, outlier_type) {
  #df <- data.frame(tissue=ordered_tissue_names, tissues_position=1:length(ordered_tissue_names), delta_auc=delta_auc_vec, lower_bound=lower_bound_vec, upper_bound=upper_bound_vec)
  error_bar_plot <- ggplot() + geom_point(data=df, mapping=aes(x=tissues_position, y=delta_auc), color="#0D324D") +
          geom_errorbar(data=df, mapping=aes(x=tissues_position,ymin=lower_bound, ymax=upper_bound),color="#0D324D",width=0.0) +
          gtex_v8_figure_theme() +
          geom_hline(yintercept = 0.0, size=.00001,linetype="dashed") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=0.5)) +
          scale_x_continuous(breaks=1:length(df$tissue), labels=df$tissue) + 
          xlab("") + 
          ylab("tissue Watershed AUC (PR) - RIVER AUC (PR)")
  
  return(error_bar_plot + draw_text(outlier_type, x=11, y=.68,size=8))
}


##############
# Input files (one for splicing, ase, and expression)
figS36_expression_input_file = "processed_input_data/figureS36/figS36_expression_input_data.txt"
figS36_ase_input_file = "processed_input_data/figureS36/figS36_ase_input_data.txt"
figS36_splicing_input_file = "processed_input_data/figureS36/figS36_splicing_input_data.txt"


##############
# Load in data
expression_df <- read.table(figS36_expression_input_file, header=TRUE, sep="\t")
ase_df <- read.table(figS36_ase_input_file, header=TRUE, sep="\t")
splicing_df <- read.table(figS36_splicing_input_file, header=TRUE, sep="\t")




################
# generate figure S32 panels
figS36_expression <- generate_figureS36_panel(expression_df, "Expression")
figS36_ase <- generate_figureS36_panel(ase_df, "ASE")
figS36_splicing <- generate_figureS36_panel(splicing_df, "Splicing")



################
# Merge panels together with cowplot
figS36 <- plot_grid(figS36_expression + ylab("        ")  + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), figS36_ase + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), figS36_splicing+ ylab("        ") +theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), legend.position="none",plot.margin=grid::unit(c(0,0,0,0), "mm")), rel_heights=c(1.24,1.24,1.9), ncol=1)
################
# Save to output file
output_file <- "generated_figures/figureS36.pdf"
ggsave(figS36,file=output_file,width=7.2, height=6.5, units="in")


