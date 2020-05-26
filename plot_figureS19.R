library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


broad_mutation_type_bar_plot <- function(inlier_distances, outlier_distances, ss_type, position, title,concensus_allele, y_axis) {
	position_outlier_distances <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance == position,]
	position_inlier_distances <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance == position,]
	# Initialize output arrays
	variant_types <- c()
	outlier_status <- c()
	counts <- c()

	###############
	# Variant changes to concensus allele
	###############
	variant_type <- "Concencus Created"
	a <- sum(position_outlier_distances$variant_allele == concensus_allele) 
	c <- sum(position_inlier_distances$variant_allele == concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "outlier")
	counts <- c(counts, a)
	# Add inliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "inlier")
	counts <- c(counts, c)

	###############
	# Variant changes from concensus allele
	###############
	variant_type <- "Concencus Destroyed"
	a <- sum(position_outlier_distances$major_allele == concensus_allele) 
	c <- sum(position_inlier_distances$major_allele == concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "outlier")
	counts <- c(counts, a)
	# Add inliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "inlier")
	counts <- c(counts, c)

	###############
	# Variant changes neither to or from concensus allele
	###############
	variant_type <- "Neither"
	a <- sum(position_outlier_distances$major_allele != concensus_allele & position_outlier_distances$variant_allele != concensus_allele) 
	c <- sum(position_inlier_distances$major_allele != concensus_allele & position_inlier_distances$variant_allele != concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "outlier")
	counts <- c(counts, a)
	# Add inliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "inlier")
	counts <- c(counts, c)



	
	# Organize into compact data frame
	df <- data.frame(variant_types=factor(variant_types,levels=c("Neither", "Concencus Destroyed", "Concencus Created")), outlier_status = factor(outlier_status), counts=counts)
	# Count number of outliers
	outlier_counts <- dim(position_outlier_distances)[1] 
	# Count number of inliers
	inlier_counts <- dim(position_inlier_distances)[1]

	plotter <- ggplot(df,aes(x=outlier_status, y=counts, fill=variant_types)) + 
    	geom_bar(stat="identity", position="fill") +
    	labs(x="",y="Proportion",fill="", title="") + 
    	scale_fill_manual(values=c("grey", "palevioletred3","skyblue3")) + 
    	gtex_v8_figure_theme()
    	#y,x
    x_pos_1 = .61
    x_pos_2 = .84
    x_pos_3 = .75
    if (y_axis == FALSE) {
    	plotter <- plotter + theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank())
        x_pos_1 = .32
        x_pos_2 = .76
        x_pos_3 = .55
    }
    plotter <- ggdraw(plotter + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme(legend.position="none")) +
     		   draw_label(title, x=x_pos_3, y=.08, size=8)
    return(plotter)
}
get_broad_mutation_type_bar_plot_legend <- function(inlier_distances, outlier_distances, ss_type, position, title,concensus_allele, y_axis) {
	position_outlier_distances <- outlier_distances[outlier_distances$splice_site_type==ss_type & outlier_distances$distance == position,]
	position_inlier_distances <- inlier_distances[inlier_distances$splice_site_type==ss_type & inlier_distances$distance == position,]
	# Initialize output arrays
	variant_types <- c()
	outlier_status <- c()
	counts <- c()

	###############
	# Variant changes to concensus allele
	###############
	variant_type <- "Consensus Created"
	a <- sum(position_outlier_distances$variant_allele == concensus_allele) 
	c <- sum(position_inlier_distances$variant_allele == concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "outlier")
	counts <- c(counts, a)
	# Add inliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "inlier")
	counts <- c(counts, c)

	###############
	# Variant changes from concensus allele
	###############
	variant_type <- "Consensus Destroyed"
	a <- sum(position_outlier_distances$major_allele == concensus_allele) 
	c <- sum(position_inlier_distances$major_allele == concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "outlier")
	counts <- c(counts, a)
	# Add inliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "inlier")
	counts <- c(counts, c)

	###############
	# Variant changes neither to or from concensus allele
	###############
	variant_type <- "Neither"
	a <- sum(position_outlier_distances$major_allele != concensus_allele & position_outlier_distances$variant_allele != concensus_allele) 
	c <- sum(position_inlier_distances$major_allele != concensus_allele & position_inlier_distances$variant_allele != concensus_allele) 
	# Add outliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "outlier")
	counts <- c(counts, a)
	# Add inliers
	variant_types <- c(variant_types, variant_type)
	outlier_status <- c(outlier_status, "inlier")
	counts <- c(counts, c)



	
	# Organize into compact data frame
	df <- data.frame(variant_types=factor(variant_types,levels=c("Consensus Created", "Consensus Destroyed", "Neither")), outlier_status = factor(outlier_status), counts=counts)
	# Count number of outliers
	outlier_counts <- dim(position_outlier_distances)[1] 
	# Count number of inliers
	inlier_counts <- dim(position_inlier_distances)[1]
	plotter <- ggplot(df,aes(x=outlier_status, y=counts, fill=variant_types)) + 
    	geom_bar(stat="identity", position="fill") +
    	labs(x="",y="proportion",fill="", title="") + 
    	scale_fill_manual(values=c("skyblue3", "palevioletred3","grey")) + 
    	gtex_v8_figure_theme()
    	#y,x

    return(get_legend(plotter + theme(legend.text = element_text(size=7),legend.key.size = unit(0.1, "cm"), legend.position="bottom")))
}


make_mutation_type_bar_plot_across_consensus_sites <- function(outlier_distances, inlier_distances) {
	# Further seperate data frame based on whether nearby splice site was annotated
	inlier_distances_novel <- inlier_distances[as.character(inlier_distances$annotated_splice_site) == "novel", ]
	outlier_distances_novel <- outlier_distances[as.character(outlier_distances$annotated_splice_site) == "novel", ]

	inlier_distances_annotated <- inlier_distances[as.character(inlier_distances$annotated_splice_site) == "annotated", ]
	outlier_distances_annotated <- outlier_distances[as.character(outlier_distances$annotated_splice_site) == "annotated", ]

	print(summary(outlier_distances_novel))

	################
	# Annotated
	################
	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	concensus <- "T"
	d_6_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	concensus <- "G"
	d_5_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	concensus <- "A"
	d_4_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	concensus <- "A"
	d_3_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -2
	ss_name <- "D+2"
	concensus <- "T"
	d_2_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -1
	ss_name <- "D+1"
	concensus <- "G"
	d_1_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	concensus <- "G"
	d_minus_1_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, TRUE)

	ss_type <- "acceptor"
	position <- -2
	ss_name <- "A-2"
	concensus <- "A"
	a_minus_2_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- -1
	ss_name <- "A-1"
	concensus <- "G"
	a_minus_1_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	concensus <- "G"
	a_1_plot_ano <- broad_mutation_type_bar_plot(inlier_distances_annotated, outlier_distances_annotated, ss_type, position, ss_name,concensus, FALSE)

	################
	# NOVEL
	################
	ss_type <- "donor"
	position <- -6
	ss_name <- "D+6"
	concensus <- "T"
	d_6_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name, concensus,  FALSE)

	ss_type <- "donor"
	position <- -5
	ss_name <- "D+5"
	concensus <- "G"
	d_5_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -4
	ss_name <- "D+4"
	concensus <- "A"
	d_4_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -3
	ss_name <- "D+3"
	concensus <- "A"
	d_3_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -2
	ss_name <- "D+2"
	concensus <- "T"
	d_2_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- -1
	ss_name <- "D+1"
	concensus <- "G"
	d_1_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "donor"
	position <- 0
	ss_name <- "D-1"
	concensus <- "G"
	d_minus_1_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, TRUE)

	ss_type <- "acceptor"
	position <- -2
	ss_name <- "A-2"
	concensus <- "A"
	a_minus_2_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- -1
	ss_name <- "A-1"
	concensus <- "G"
	a_minus_1_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)

	ss_type <- "acceptor"
	position <- 0
	ss_name <- "A+1"
	concensus <- "G"
	a_1_plot <- broad_mutation_type_bar_plot(inlier_distances_novel, outlier_distances_novel, ss_type, position, ss_name,concensus, FALSE)




	mutation_bar_plot_legend <- get_broad_mutation_type_bar_plot_legend(inlier_distances_novel, outlier_distances_novel, "acceptor", 0, "A+1", "G", TRUE)

	gg_novel <- plot_grid(d_minus_1_plot, d_1_plot, d_2_plot, d_3_plot, d_4_plot, d_5_plot, d_6_plot, a_minus_2_plot, a_minus_1_plot, a_1_plot, nrow=1, rel_widths=c(1.7,1, 1, 1, 1, 1, 1, 1, 1, 1))
	gg_ano <- plot_grid(d_minus_1_plot_ano, d_1_plot_ano, d_2_plot_ano, d_3_plot_ano, d_4_plot_ano, d_5_plot_ano, d_6_plot_ano, a_minus_2_plot_ano, a_minus_1_plot_ano, a_1_plot_ano, nrow=1, rel_widths=c(1.7,1, 1, 1, 1, 1, 1, 1, 1, 1))

	combined_gg <- plot_grid(gg_ano, gg_novel, ncol=1)+
					draw_label("Annotated Splice Site", x=.5,y=.98,size=9) + 
					draw_label("Novel Splice Site", x=.5,y=.49,size=9) 

	return(plot_grid(combined_gg, mutation_bar_plot_legend, nrow=2, rel_heights=c(1,.05)))
}


##############
# Input data
outlier_variant_input_file = "processed_input_data/figureS19/figS19_outliers_input_data.txt"
inlier_variant_input_file = "processed_input_data/figureS19/figS19_inliers_input_data.txt"

##############
# Load in data
outlier_variant_df <- read.table(outlier_variant_input_file, header=TRUE)
inlier_variant_df <- read.table(inlier_variant_input_file, header=TRUE)


##############
# make plot
figS19 <- make_mutation_type_bar_plot_across_consensus_sites(outlier_variant_df, inlier_variant_df)





###################
# Save to output file
output_file <- "generated_figures/figureS19.pdf"
ggsave(figS19, file=output_file, width=7.2, height=3, units="in")
