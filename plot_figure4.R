library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
library(reshape)
theme_set(theme_cowplot())


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

make_panel_4b <- function(theta_pair) {
    # Initialize mat to keep track of edge weights
    # first row corresponds to expression
    # secnod row corresponds to ASE
    # thrid row corresponds to splicing
	mat <- matrix(0, 3, 3)

    # Fill in mat
    # Splicing-ASE pair
	mat[3,2] <- theta_pair[1,2] 
	mat[2,3] <- theta_pair[1,2]
    # Expression-ASE pair
	mat[1,2] <- theta_pair[1,3] 
	mat[2,1] <- theta_pair[1,3]
    # Splicing-Expression pair
	mat[1,3] <- theta_pair[1,1] 
	mat[3,1] <- theta_pair[1,1]

    # Melt mat to easily make geom_tile
	melted_corr <- melt(mat)	
    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1)
    melted_corr$X2 <- factor(melted_corr$X2)

    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) #+ scale_fill_gradient(low="grey",high="plum2")
    heatmap <- heatmap + scale_fill_distiller(palette = "Blues", direction=1, na.value = "white")
    heatmap <- heatmap + theme(axis.text.x = element_text(angle = 0, vjust=.5),legend.position="bottom") 
    heatmap <- heatmap + gtex_v8_figure_theme()
    heatmap <- heatmap + labs(fill="Edge weight",x = "", y="")
    heatmap <- heatmap + scale_x_discrete(breaks=c("1", "2", "3"),labels=c("Expresssion", "ASE", "Splicing"))
    heatmap <- heatmap + scale_y_discrete(breaks=c("1", "2", "3"),labels=c("Expression", "ASE", "Splicing"))
    # Extract legend
    legend <- get_legend(heatmap)

    # Merge plot and legend in specified orientation using cowplot
    combined <- ggdraw() + draw_plot(heatmap+theme(legend.position="none"),0,.15,1,.85) + draw_plot(legend,.2,-.33,1,1)

	return(combined)
}

make_panel_4c <- function(df) {
	# Set order of outliers
	df$outlier_type=factor(df$outlier_type, levels=c("eOutlier", "aseOutlier","sOutlier"))
    # Set order of model types
	df$model_type=factor(df$model_type, levels=c("GAM", "Watershed"))
	
    # Make barplot
	p <- ggplot(data=df, aes(x=model_type, y=absolute_risk, fill=outlier_type)) +
		geom_bar(stat="identity", color="black", position=position_dodge())+
  		gtex_v8_figure_theme() + 
  		labs(x="", y="Proportion of variants\nleading to outlier", fill="") + 
  		scale_fill_manual(values=c("#BFCDE0", "#7F5A83", "#0D324D"))

    # Extract legend
  	legend <- get_legend(p)
    # Merge plot and legend in specified orientation using cowplot
  	combined <- ggdraw() + draw_plot(p + theme(legend.position="none"),0, -.1, 1,1.1) + draw_plot(legend, .27, .33, 1,1)

	return(combined)
}

make_panel_4d <- function(df) {
    # Set ordering of model types
	df$prediction_type=factor(df$prediction_type, levels=c("Watershed","RIVER", "GAM"))

    # Make Precision-recall curve for ASE
    # This plot contains 3 curves (one for Watershed, RIVER, and GAM)
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

    # Make Precision-recall curve for splicing
    # This plot contains 3 curves (one for Watershed, RIVER, and GAM)
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

    # Make Precision-recall curve for Expression
    # This plot contains 3 curves (one for Watershed, RIVER, and GAM)
	outlier_type <- "total_expression"
	outlier_name <- "Expression"
  	plotter_expression <- ggplot(data=df[as.character(df$outlier_type)==outlier_type,], aes(x=recall, y=precision, group=prediction_type)) + geom_line(aes(linetype=prediction_type, colour=prediction_type)) + 
                labs(x="Recall", y="Precision", group="") +
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
                theme(legend.position="bottom") +
                theme(panel.spacing = unit(2, "lines")) +
                scale_color_manual(values=c("steelblue3", "steelblue3", "firebrick4")) +
                scale_linetype_manual(values=c("solid", "dotted", "solid")) +
                gtex_v8_figure_theme() + 
                draw_label(outlier_name,x=.5,y=.95,size=8)

    # Get legend from ASE Plot
    legend <- get_legend(plotter_ase + theme(legend.position="bottom"))

    # Combine PR-curves from Expression, ASE, and splicing into 1 plot via cowplot
    combined_plots <- plot_grid(plotter_expression + theme(legend.position="none"), plotter_ase + theme(legend.position="none"), plotter_splice+ theme(legend.position="none"), rel_widths=c(1,1,1), nrow=1)

    # Merge combined plot with single legend via cowplot
	combined <- ggdraw() + draw_plot(combined_plots,0,.07,1,.9) + draw_plot(legend,.38,-.40,1,1)
	
	return(combined)
}

make_panel_4f <- function(df) {
	# Set order of Outlier types
	df$outlier_type=factor(df$outlier_type, levels=c("Expression", "ASE","Splicing"))
    # Set order of modles
	df$method=factor(df$method,levels=c("GAM", "RIVER", "Watershed"))

	# Make boxplot corresponding to AUC across tissues
	boxplot <- ggplot(df, aes(x=method, y=auc, fill=outlier_type)) + geom_boxplot() +
				gtex_v8_figure_theme() + 
				xlab("Tissue Model") +
               	ylab("AUC(PR)") + 
               	theme(legend.position="none") +
               	scale_fill_manual(values=c("#BFCDE0", "#7F5A83", "#0D324D"))  

	return(boxplot)
}

load_in_gtex_tissue_colors <- function(tissue_colors_input_file) {
	# Read in tissue colors and names
	tissue_colors = read.table(tissue_colors_input_file, header = T, stringsAsFactors = F, sep = "\t")
	# correction to minor mislabeling
	for (tiss_num in 1:length(tissue_colors$tissue_id)) {
		if (tissue_colors$tissue_id[tiss_num] == "Brain_Spinal_cord_cervical_c1") {
			tissue_colors$tissue_id[tiss_num] = "Brain_Spinal_cord_cervical_c.1"
		}
		if (tissue_colors$tissue_id[tiss_num] == "Cells_EBVtransformed_lymphocytes") {
			tissue_colors$tissue_id[tiss_num] = "Cells_EBV.transformed_lymphocytes"
		}
	}
	return(tissue_colors)
}

load_in_edge_weights <- function(input_file) {
    aa <- read.table(input_file, header=TRUE, sep="\t")
    colnames(aa) = rownames(aa)
    return(as.matrix(aa))
}


make.point.plot = function(tissuesdf, colors, vertical = TRUE){
    if (vertical) {
        p = ggplot(tissuesdf, aes(x = 1, y = order, label = tissue_id))
    } else {
        p = ggplot(tissuesdf, aes(x = order, y = 1))
    }
    p = p + geom_point(aes(colour = tissue_id), size = .4) +
        scale_colour_manual(values = colors) + guides(colour = FALSE) + 
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank())
    return(p)
}

make_panel_4e <- function(theta_pair_mat, tissue_colors) {
	# Get ordered names of tissues according theta_pair_mat (theta_pair_mat is the edge weight matrix)
    tissue_names <- colnames(theta_pair_mat)
    # Perform hierarchial clustering on tissues according to their edge weights
	order <- hclust( dist(theta_pair_mat, method = "euclidean"), method = "ward.D" )$order

	#######################
    # Reformat tissue colors
	tissue_colors$tissue_id = factor(tissue_colors$tissue_id, levels = as.character(tissue_names[order]))
	#print(tissue_colors$tissue_id)
	tissue_colors$order = as.numeric(tissue_colors$tissue_id)
	tissue_colors = tissue_colors[!is.na(tissue_colors$order),]
	colors = paste0("#",tissue_colors$tissue_color_hex)
	names(colors) = tissue_colors$tissue_id


    # Make dot plots (that will be the boarders of the heatmap) corresponding to tissue types
	colors.vertical = make.point.plot(tissue_colors, colors)
	colors.horizontal = make.point.plot(tissue_colors, colors, vertical = FALSE)

    # Add row and column names to theta_pair_mat
	rownames(theta_pair_mat) <- tissue_names
	colnames(theta_pair_mat) <- tissue_names

    # Melt theta_pair_mat to make data compatible with geom_tile
	melted_mat <- melt(theta_pair_mat)
    # Axis labels are factors
    melted_mat$X1 <- factor(melted_mat$X1)
    melted_mat$X2 <- factor(melted_mat$X2)

    #  Order tissue factors by their hierarchial clustering order
    melted_mat$X1 <- factor(melted_mat$X1, levels = rownames(theta_pair_mat)[order])
    melted_mat$X2 <- factor(melted_mat$X2, levels = rownames(theta_pair_mat)[order])


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    heatmap <- heatmap + scale_fill_gradient2() #+ scale_fill_distiller(palette="RdBu")
    heatmap <- heatmap + theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5,size=7), axis.text.y = element_text(size=7))
    heatmap <- heatmap + gtex_v8_figure_theme() + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
    heatmap <- heatmap + labs(x = "Tissue", y = "Tissue",fill="Edge Weight")
    
    # Merge heatmap with tissue-colors dot-plots via Cowplot
   	combined = ggdraw() +
        draw_plot(heatmap, .05,.05,.95,0.95) +
        draw_plot(colors.vertical, -.344, .084, .903, .946) + 
        draw_plot(colors.horizontal, 0.063,-.192, .76,.69)

    return(combined)



}

##############
# Input data
fig4b_input_file = "processed_input_data/figure4/figure4b_watershed_edge_weights_input_data.txt"
fig4c_input_file = "processed_input_data/figure4/figure4c_absolute_risk_input_data.txt"
fig4d_input_file = "processed_input_data/figure4/figure4d_pr_curve_input_data.txt"
fig4e_input_file = "processed_input_data/figure4/figure4e_tissue_watershed_edge_weights_input_data.txt"
tissue_colors_input_file = "processed_input_data/figure4/figure4e_gtex_tissue_colors_input_data.txt"
fig4f_input_file = "processed_input_data/figure4/figure4f_tissue_watershed_pr_auc_input_data.txt"

##############
# Load in data
# File contains mapping from GTEx tissue name to HEX color (standardized for gtex)
gtex_tissue_colors <- load_in_gtex_tissue_colors(tissue_colors_input_file)
# Input files for each panel
watershed_edge_weights <- read.table(fig4b_input_file, header=FALSE, sep="\t")
fig4c_df <- read.table(fig4c_input_file, header=TRUE, sep="\t")
fig4d_df <- read.table(fig4d_input_file, header=TRUE, sep="\t")
tissue_watershed_edge_weights <- load_in_edge_weights(fig4e_input_file)
fig4f_df <- read.table(fig4f_input_file, header=TRUE, sep="\t")


##############
# make panels
fig4b <- make_panel_4b(watershed_edge_weights)
fig4c <- make_panel_4c(fig4c_df)
fig4d <- make_panel_4d(fig4d_df)
fig4e <- make_panel_4e(tissue_watershed_edge_weights, gtex_tissue_colors)
fig4f <- make_panel_4f(fig4f_df)

###############
# Combine panels with cowplot
row_1 <- plot_grid(NULL, fig4b, fig4c, ncol=3, labels=c("A","B","C"), rel_widths=c(1,1.1,1))
row_2 <- plot_grid(fig4d, labels=c("D"))
row_3 <- plot_grid(fig4e, fig4f,ncol=2, rel_widths=c(1,.75), labels=c("E", "F"))
fig4 <- plot_grid(row_1, row_2, row_3, ncol=1, rel_heights=c(.55,.5,.7))

###################
# Save to output file
output_file <- "generated_figures/figure4.pdf"
ggsave(fig4, file=output_file, width=7.2, height=6.0, units="in")
