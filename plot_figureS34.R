library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
library(reshape)
theme_set(theme_cowplot())


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

generate_figureS34_panel <- function(df, outlier_type) {
  error_bar_plot <- ggplot() + geom_point(data=df, mapping=aes(x=tissues_position, y=delta_auc), color="#0D324D") +
          geom_errorbar(data=df, mapping=aes(x=tissues_position,ymin=lower_bound, ymax=upper_bound),color="#0D324D",width=0.0) +
          gtex_v8_figure_theme() +
          geom_hline(yintercept = 0.0, size=.00001,linetype="dashed") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=0.5)) +
          scale_x_continuous(breaks=1:length(df$tissue), labels=df$tissue) + 
          xlab("") + 
          ylab("tissue Watershed AUC (PR) - tissue RIVER AUC (PR)")
  
  return(error_bar_plot + draw_text(outlier_type, x=11, y=.68,size=8))
}


##############
# Input files (one for splicing, ase, and expression)
figS34_expression_input_file = "processed_input_data/figureS34/figS34_expression_input_data.txt"
figS34_ase_input_file = "processed_input_data/figureS34/figS34_ase_input_data.txt"
figS34_splicing_input_file = "processed_input_data/figureS34/figS34_splicing_input_data.txt"


##############
# Load in data
expression_df <- read.table(figS34_expression_input_file, header=TRUE, sep="\t")
ase_df <- read.table(figS34_ase_input_file, header=TRUE, sep="\t")
splicing_df <- read.table(figS34_splicing_input_file, header=TRUE, sep="\t")




################
# generate figure S32 panels
figS34_expression <- generate_figureS34_panel(expression_df, "Expression")
figS34_ase <- generate_figureS34_panel(ase_df, "ASE")
figS34_splicing <- generate_figureS34_panel(splicing_df, "Splicing")



################
# Merge panels together with cowplot
figS34 <- plot_grid(figS34_expression + ylab("        ")  + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), figS34_ase + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="none",axis.text.x=element_blank()), figS34_splicing+ ylab("        ") +theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), legend.position="none",plot.margin=grid::unit(c(0,0,0,0), "mm")), rel_heights=c(1.24,1.24,1.9), ncol=1)

################
# Save to output file
output_file <- "generated_figures/figureS34.pdf"
ggsave(figS34,file=output_file,width=7.2, height=6.5, units="in")


