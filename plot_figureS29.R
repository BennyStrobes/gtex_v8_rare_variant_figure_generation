library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)
library(reshape)
theme_set(theme_cowplot())


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

load_in_confusion_matrix_for_method <- function(input_file) {
    aa <- read.table(input_file, header=TRUE, sep="\t")
    colnames(aa) = rownames(aa)
    return(as.matrix(aa))
}

visualize_confusion_matrix <- function(confusion_matrix, titler) {
    confusion_matrix = confusion_matrix/rowSums(confusion_matrix)
    melted_corr <- melt(confusion_matrix)

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1, levels=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"))
    melted_corr$X2 <- factor(melted_corr$X2, levels=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"))
    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) #+ scale_fill_gradient(low="grey",high="plum2")

    heatmap <- heatmap + scale_fill_distiller(palette = "Blues", direction=1, limits=c(0,1))

    #heatmap <- heatmap + theme(text = element_text(size=12),axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11), axis.text.x = element_text(angle = 0, vjust=.5)) 
    heatmap <- heatmap + labs(x = "Observed Class", y = "Predicted Class",fill="",title=titler)
    heatmap <- heatmap + gtex_v8_figure_theme() + theme(axis.text.x = element_text(angle = 0, vjust=.5))
    heatmap <- heatmap + scale_x_discrete(breaks=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"),labels=c("0 0 0", "1 0 0", "0 1 0", "0 0 1", "1 0 1","0 1 1", "1 1 0", "1 1 1"))
    heatmap <- heatmap + scale_y_discrete(breaks=c("1 1 1", "2 1 1", "1 2 1", "1 1 2", "2 1 2", "1 2 2", "2 2 1", "2 2 2"),labels=c("0 0 0", "1 0 0", "0 1 0", "0 0 1", "1 0 1","0 1 1", "1 1 0", "1 1 1"))

    return(heatmap)
}


##############
# Input files
figS29_river_input_file = "processed_input_data/figureS29/figS29_river_input_data.txt"
figS29_watershed_exact_input_file = "processed_input_data/figureS29/figS29_watershed_exact_input_data.txt"
figS29_watershed_approximate_input_file = "processed_input_data/figureS29/figS29_watershed_approximate_input_data.txt"


##############
# Load in confusion matrices for each method
river_confusion_matrix <- load_in_confusion_matrix_for_method(figS29_river_input_file)
watershed_exact_confusion_matrix <- load_in_confusion_matrix_for_method(figS29_watershed_exact_input_file)
watershed_approximate_confusion_matrix <- load_in_confusion_matrix_for_method(figS29_watershed_approximate_input_file)


################
# generate figure S29 panels
river_confusion_plot <- visualize_confusion_matrix(river_confusion_matrix, "RIVER")
watershed_confusion_plot <- visualize_confusion_matrix(watershed_exact_confusion_matrix, "Watershed (Exact)")
watershed_approximate_confusion_plot <- visualize_confusion_matrix(watershed_approximate_confusion_matrix, "Watershed (Approximate)")



################
# Merge panels together with cowplot
figS29 <- plot_grid(river_confusion_plot, watershed_confusion_plot, watershed_approximate_confusion_plot, ncol=1)

################
# Save to output file
output_file <- "generated_figures/figureS29.pdf"
ggsave(figS29, file=output_file, width=7.2, height=7.0, units="in")
