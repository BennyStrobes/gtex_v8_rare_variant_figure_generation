library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
theme_set(theme_cowplot())


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

make_absolute_risk_panel <- function(df, watershed_threshold) {
    # Set ordering of outlier types
    df$outlier_type = factor(df$outlier_type, levels=c("eOutlier", "aseOutlier", "sOutlier"))
    # plot figure
    p <- ggplot(data=df, aes(x=model_type, y=absolute_risk, fill=outlier_type)) +
        geom_bar(stat="identity", color="black", position=position_dodge())+
        ylim(0,.72) + 
        gtex_v8_figure_theme() + 
        labs(x="", y="Proportion of variants\nleading to outlier", fill="", title=paste0("Watershed posterior > ", watershed_threshold)) + 
        scale_fill_manual(values=c("#BFCDE0", "#7F5A83", "#0D324D")) 
    return(p)
}

##############
# Input files
figS30_input_file = "processed_input_data/figureS30/figS30_input_data.txt"


##############
# Load in data
figS30_df <- read.table(figS30_input_file, header=TRUE, sep="\t")

################
# generate figure S30 panels
threshold=0.5
figS30_a <- make_absolute_risk_panel(figS30_df[figS30_df$watershed_threshold==threshold, ], threshold)
threshold=0.7
figS30_b <- make_absolute_risk_panel(figS30_df[figS30_df$watershed_threshold==threshold, ], threshold)
threshold=0.9
figS30_c <- make_absolute_risk_panel(figS30_df[figS30_df$watershed_threshold==threshold, ], threshold)

legend <- get_legend(figS30_a + theme(legend.position="bottom"))

################
# Merge panels together with cowplot
figS30 <- plot_grid(figS30_a + theme(legend.position="none"), figS30_b + theme(legend.position="none"), figS30_c + theme(legend.position="none"), legend, labels=c("A","B","C"), rel_heights=c(1,1,1,.13), nrow=4)

################
# Save to output file
output_file <- "generated_figures/figureS30.pdf"
ggsave(figS30, file=output_file, width=7.2, height=6.0, units="in")
