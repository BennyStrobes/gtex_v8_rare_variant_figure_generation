library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)
library(ggseqlogo)


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}



mutation_type_pwm_plot_across_novel_concensus_sites <- function(outlier_distances) {
	alleles <- c("A", "C", "T", "G")
	# Seperate mutations string into reference and alternate allels
	mutations <- as.character(outlier_distances$mutations)
	refs <- c()
	alts <- c()
	for (iter in 1:length(mutations)) {
		mutation <- mutations[iter]
		ref <- strsplit(mutation, "->")[[1]][1]
		alt <- strsplit(mutation, "->")[[1]][2]
		refs <- c(refs, ref)
		alts <- c(alts, alt)
	}

	#####################################################
	# Outlier variant PWM plot for Donor splice sites
	#####################################################
	ss_names <- c("D-1","D+1", "D+2", "D+3", "D+4", "D+5", "D+6")
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		variant_alleles <- alts[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio > 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(variant_alleles==allele)
		}
	}
	
	outlier_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier variants at donor sites") + 
	        labs(x="",y="Frequency", title="") + 
	        gtex_v8_figure_theme() +
	       	theme(plot.title = element_text(hjust = 0.5)) +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names)


	#####################################################
	# Background variant PWM plot for Donor splice sites
	#####################################################
	ss_names <- c("D-1","D+1", "D+2", "D+3", "D+4", "D+5", "D+6")
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		outlier_background_alleles <- refs[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio > 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(outlier_background_alleles==allele)
		}
	}

	background_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Background at donor sites") + 
	        labs(x="",y="Frequency", title="") + 
	        theme(plot.title = element_text(hjust = 0.5)) +
	        gtex_v8_figure_theme() +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names)
	
	combined_donor_plotter <- plot_grid(background_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), outlier_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),nrow=2)



	#####################################################
	# Outlier variant PWM plot for acceptor splice sites
	#####################################################
	ss_names <- c("A-2","A-1", "A+1")
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		variant_alleles <- alts[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio > 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(variant_alleles==allele)
		}
	}
	
	outlier_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier variants at acceptor sites") + 
	        labs(x="",y="Frequency", title="") + 
	        gtex_v8_figure_theme() +
	       	theme(plot.title = element_text(hjust = 0.5)) +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names) +
	        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(),axis.ticks.y=element_blank())


	#####################################################
	# Background variant PWM plot for acceptor splice sites
	#####################################################
	ss_names <- c("A-2","A-1", "A+1")
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles


	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		outlier_background_alleles <- refs[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio > 1.0]

		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(outlier_background_alleles==allele)
		}
	}

	background_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        labs(x="",y="Frequency", title="") + 
	        #labs(x="",y="Frequency", title="Background at acceptor sites") + 
	        theme(plot.title = element_text(hjust = 0.5)) +
	        gtex_v8_figure_theme() +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names) +
	        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(),axis.ticks.y=element_blank())

	
	combined_acceptor_plotter <- plot_grid(background_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), outlier_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),nrow=2)

	# MERGE ALL 4 plots together
	combined <- plot_grid(combined_donor_plotter, combined_acceptor_plotter, ncol=2, rel_widths=c(7,3)) +
				draw_text("Background for splice sites with high junction usage", x=.55, y=.97, size=8)+
				draw_text("Outlier variants for splice sites with high junction usage", x=.55, y=.47, size=8)

	return(combined)
}


##############
# Input data for figure S17A
variant_type_input_file = "processed_input_data/figureS18/figS18_input_data.txt"




##############
# Load in variant enrichment data
variant_type_df <- read.table(variant_type_input_file, header=TRUE, sep="\t")




##############
# plot S17a
figS18 <- mutation_type_pwm_plot_across_novel_concensus_sites(variant_type_df)



# Save to output file
output_file <- "generated_figures/figureS18.pdf"
ggsave(figS18, file=output_file, width=7.2, height=3, units="in")
