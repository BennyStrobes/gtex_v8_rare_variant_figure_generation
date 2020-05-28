library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
library(reshape)
theme_set(theme_cowplot())


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

load_in_theta_matrix <- function(input_file) {
    aa <- read.table(input_file, header=TRUE, sep="\t")
    colnames(aa) = rownames(aa)
    return(as.matrix(aa))
}

generate_theta_pair_heatmap <- function(theta_pair_mat, outlier_type) {
    order <- hclust( dist(theta_pair_mat, method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(theta_pair_mat)
    # Axis labels are factors
    melted_mat$X1 <- factor(melted_mat$X1)
    melted_mat$X2 <- factor(melted_mat$X2)


    #  Use factors to represent covariate and pc name
    melted_mat$X1 <- factor(melted_mat$X1, levels = rownames(theta_pair_mat)[order])
    melted_mat$X2 <- factor(melted_mat$X2, levels = rownames(theta_pair_mat)[order])


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    heatmap <- heatmap + scale_fill_gradient2() #+ scale_fill_distiller(palette="RdBu")
    heatmap <- heatmap + theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5,size=6), axis.text.y = element_text(size=6))
    heatmap <- heatmap + gtex_v8_figure_theme() + theme(axis.text.x=element_blank())
    heatmap <- heatmap + labs(x = "Tissue", y = "Tissue",fill="Edge Weight", title=outlier_type)
    return(heatmap)

}

##############
# Input files (one for splicing, ase, and expression)
figS32_expression_input_file = "processed_input_data/figureS32/figS32_expression_input_data.txt"
figS32_ase_input_file = "processed_input_data/figureS32/figS32_ase_input_data.txt"
figS32_splicing_input_file = "processed_input_data/figureS32/figS32_splicing_input_data.txt"


##############
# Load in data
expression_theta_matrix <- load_in_theta_matrix(figS32_expression_input_file)
ase_theta_matrix <- load_in_theta_matrix(figS32_ase_input_file)
splicing_theta_matrix <- load_in_theta_matrix(figS32_splicing_input_file)



################
# generate figure S32 panels
figS32_expression <- generate_theta_pair_heatmap(expression_theta_matrix, "Expression")
figS32_ase <- generate_theta_pair_heatmap(ase_theta_matrix, "ASE")
figS32_splicing <- generate_theta_pair_heatmap(splicing_theta_matrix, "Splicing")


################
# Merge panels together with cowplot
figS32 <- plot_grid(figS32_expression, figS32_ase, figS32_splicing, ncol=1)

################
# Save to output file
output_file <- "generated_figures/figureS32.pdf"
ggsave(figS32, file=output_file, width=7.2, height=11.0, units="in")


