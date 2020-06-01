library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
library(reshape)
theme_set(theme_cowplot())


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

generate_figureS37_panel <- function(df, outlier_type) {
  high_watershed_label <- paste0("GTEx Watershed posterior > .8")
  med_watershed_label <- paste0(".5 > GTEx Watershed posterior <= .8")
  low_watershed_label <- paste0("GTEx Watershed posterior < .01")

  df2 <- df[as.character(df$type) != low_watershed_label, ]


  # Set order
  df$type = factor(df$type, levels=c(low_watershed_label, med_watershed_label, high_watershed_label))
  df$model=factor(df$model, levels=c("Watershed", "GAM"))
  df2$type = factor(df2$type, levels=c(med_watershed_label, high_watershed_label))
  df2$model=factor(df2$model, levels=c("Watershed", "GAM"))
    p <- ggplot() + 
      geom_boxplot(data=df, aes(x=type, y=pval,fill=model), outlier.size=.1, alpha=0.22) +
      geom_dotplot(data=df2, aes(x=type, y=pval, fill=model), binaxis='y', stackdir='center', position=position_dodge(width=.75)) +
      scale_fill_manual(values=c("steelblue3", "firebrick4")) + 
      scale_colour_manual(values=c("steelblue3", "firebrick4")) + 
      gtex_v8_figure_theme() + 
      labs(title=outlier_type,fill="", colour="",x="", y=paste0(outlier_type, " median -log10(p-value)\n in ASMAD cohort"))
    return(p) 

}


##############
# Input files (one for splicing, ase, and expression)
figS37_expression_input_file = "processed_input_data/figureS37/figS37_expression_input_data.txt"
figS37_ase_input_file = "processed_input_data/figureS37/figS37_ase_input_data.txt"
figS37_splicing_input_file = "processed_input_data/figureS37/figS37_splicing_input_data.txt"


##############
# Load in data
expression_df <- read.table(figS37_expression_input_file, header=TRUE, sep="\t")
ase_df <- read.table(figS37_ase_input_file, header=TRUE, sep="\t")
splicing_df <- read.table(figS37_splicing_input_file, header=TRUE, sep="\t")



################
# generate figure S37 panels
figS37_expression <- generate_figureS37_panel(expression_df, "Expression")
figS37_ase <- generate_figureS37_panel(ase_df, "ASE")
figS37_splicing <- generate_figureS37_panel(splicing_df, "ASE")



################
# Merge panels together with cowplot
legend <- get_legend(figS37_splicing + theme(legend.position="bottom"))
combined_plots <- plot_grid(figS37_expression + theme(legend.position="none"), figS37_ase + theme(legend.position="none"), figS37_splicing+ theme(legend.position="none"), ncol=1)
figS37 <- ggdraw() + draw_plot(combined_plots,0,.07,1,.9) + draw_plot(legend,.42,-0.43,1,1)


################
# Save to output file
output_file <- "generated_figures/figureS37.pdf"
ggsave(figS37,file=output_file,width=7.2, height=7.0, units="in")


