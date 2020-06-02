library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(cowplot)
theme_set(theme_cowplot())
library(ggseqlogo)


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=8, hjust=0.5), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

make_panel_2a <- function(dist.data) {
	# Use these colors eOutliers, aseOutliers, and sOutliers, respectively
	dcols = c('#BFCDE0', '#7F5A83', '#0D324D')

	# Specify ordering of outliers
	dist.data$Method = factor(dist.data$Method, levels=c('eOutlier', 'aseOutlier', 'sOutlier'))
	# Specify ordering of variant types
	dist.data$Type = factor(dist.data$Type, levels=c('SNVs', 'indels', 'SVs'))

	# Filter data frame to SNVs
	dist.data.filtered1 = dist.data[dist.data$Riskratio > 1 & dist.data$Type=='SNVs',]
	# Make distance plot for SNVs
	fig_2A1 = ggplot(dist.data.filtered1, aes(x=Window, y=Riskratio,Group=Method)) +
		geom_point(size=1.7,aes(color=Method), position=position_dodge(width=0.5)) +
		geom_errorbar(aes(ymin = Lower, ymax = Upper),size=.2, width=0,position=position_dodge(width=0.5)) +
		theme_bw() + ylab('Relative risk') + xlab('') +
		geom_line(aes(group=Method),size=0.2) +
		geom_hline(yintercept=1,color='grey',size=.2) +
		scale_color_manual(values=dcols) + guides(color=F) +
		theme(axis.text.x=element_text(angle=45,hjust=1),strip.text.x=element_text(size=8)) +
		gtex_v8_figure_theme() +
		theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
			panel.border = element_blank()) +
		facet_wrap(Type~.,scales='free', ncol=3) + theme(strip.background = element_blank())

	# Filter data frame to Indels and SVs
	dist.data.filtered2 = dist.data[dist.data$Riskratio > 1 & dist.data$Type!='SNVs',]
	# Make distance plot for indels and SVs
	fig_2A2 = ggplot(dist.data.filtered2, aes(x=Window, y=Riskratio,Group=Method)) +
		geom_point(size=1.7,aes(color=Method), position=position_dodge(width=0.5)) +
		geom_errorbar(aes(ymin = Lower, ymax = Upper),size=.2, width=0,position=position_dodge(width=0.5)) +
		theme_bw() + ylab('') + geom_line(aes(group=Method),size=0.2) +
		#xlab('Distance upstream from gene (kb)') + ylab('') +
		xlab('') + ylab('') +
		scale_y_log10() + geom_hline(yintercept=1,size=.2,color='grey') +
		scale_color_manual(values=dcols) +
		theme(axis.text.x=element_text(angle=45,hjust=1),
        	axis.title.x=element_text(hjust=0),
        	strip.text.x=element_text(size=8),
        	legend.position=c(0.83,0.89),
        	panel.border = element_blank(),
        	legend.key.size=unit(0.05,'in')) +
        gtex_v8_figure_theme() +
        facet_wrap(Type~.,scales='free', ncol=3) + 
        theme(plot.margin = unit(c(0, .5, 0, 0), "cm")) +
        theme(legend.title=element_blank(), strip.background = element_blank())

    # Merge 3 distance plots (one for SNVs, indels, and SVs) together using cowplot
    fig_2A <- plot_grid(fig_2A1, fig_2A2, ncol=2, rel_widths=c(.96,2)) + draw_text('Distance upstream from gene (kb)', x=.53, y=.04,size=8)

   return(fig_2A)
}

make_panel_2b <- function(plot.tss.data) {
	# Specify colors for each motif
	vcols = c(brewer.pal(11,'Spectral')[1:7],'grey',brewer.pal(11,'Spectral')[8:11])
	names(vcols) = unique(plot.tss.data$promoter_motif)
	# Set ordering of eOutlier Bins
	plot.tss.data$ExpBin = factor(plot.tss.data$ExpBin, levels=c('Under', 'Control', 'Over'))
	# Set ordering of motifs
	plot.tss.data$promoter_motif = factor(plot.tss.data$promoter_motif, 
                                      levels=rev(c('Other motif', 'Cmyc', 'CTCF', 'E2F4', 'E2F6', 'ELF1', 'Gabp', 'Nrf1', 'Nrsf', 'PU1', 'SP1', 'Srf', 'Tr4', 'USF1', 'Yy1')))
	# Make stacked barplot
	fig_2B = ggplot(plot.tss.data, aes(x=ExpBin,y=NumBin)) + 
  		geom_bar(aes(fill=promoter_motif),stat='identity',color='black') + theme_bw() +
  		scale_fill_manual(values=vcols) + xlab('eOutlier bin') + 
  		ylab('Proportion with RV in promoter') +
  		theme(panel.border = element_blank(),
   	    	legend.key.height = unit(0.05, "in")) +
  		gtex_v8_figure_theme() + theme(legend.title=element_blank(),
                                 legend.position=c(0.5,0.6))
  	return(fig_2B)
}

# Get Odds ratio corresponding to every position around splice sites
extract_positional_odds_ratio_data_for_all_positions_in_window <- function(outlier_distances, inlier_distances, intron_start_pos, exon_start_pos) {
	dist_to_ss <- c()
	odds_ratios <- c()
	upper_bounds <- c()
	lower_bounds <- c()
	# Loop through each position
	for (distance_iter in -intron_start_pos:(exon_start_pos-1)) {
		# Compute elements of contingency table at this position
		rv_outlier <- sum(outlier_distances$distance == distance_iter)
		rv_inlier <- sum(inlier_distances$distance == distance_iter)
		no_rv_outlier <- sum(outlier_distances$distance != distance_iter)
		no_rv_inlier <- sum(inlier_distances$distance != distance_iter)

		# Add pseudo-count of 1 to each element of contingency table
		a <- rv_outlier  + 1
		b <- rv_outlier + no_rv_outlier + 2
		c <- rv_inlier +1
		d <- rv_inlier + no_rv_inlier + 2
		# Compute odds ratio
		orat <- (a/b)/(c/d)

		# Get 95% CIs
		log_bounds <- 1.96*sqrt((1.0/a) - (1.0/b) + (1.0/c) - (1.0/d))
		upper_bound <- orat*exp(log_bounds)
		lower_bound <- orat*exp(-log_bounds)

		# Add data to array
		odds_ratios <- c(odds_ratios, orat)
		dist_to_ss <- c(dist_to_ss, distance_iter)
		upper_bounds <- c(upper_bounds, upper_bound)
		lower_bounds <- c(lower_bounds, lower_bound)
	}
	# Store in compact data frame
	df <- data.frame(odds_ratio=odds_ratios, dist_to_ss=dist_to_ss, lower_bounds=lower_bounds, upper_bounds=upper_bounds)

	return(df)

}

make_panel_2d <- function(inlier_distances, outlier_distances) {
	# Create seperate a seperate data frame for the following 4 classes:
	### 1. sOutliers with variant nearest to donor splice site
	outlier_donor_distances <- outlier_distances[as.character(outlier_distances$splice_site_type)=="donor",]
	### 2. sOutliers with variant nearest to acceptor splice site
	outlier_acceptor_distances <- outlier_distances[as.character(outlier_distances$splice_site_type)=="acceptor",]
	### 3. sInliers with variant nearest to donor splice site
	inlier_donor_distances <- inlier_distances[as.character(inlier_distances$splice_site_type)=="donor",]
	### 4. sInliers with variant nearest to acceptor splice site
	inlier_acceptor_distances <- inlier_distances[as.character(inlier_distances$splice_site_type)=="acceptor",]

	# Get Odds ratio corresponding to every position in Donor concensus splice sites
	# 7 is the intron start position (D+7) and 3 is the exon start position (D-3)
	donor_df <- extract_positional_odds_ratio_data_for_all_positions_in_window(outlier_donor_distances, inlier_donor_distances, 7, 3)
	# Get Odds ratio corresponding to every position in Acceptor concensus splice sites
	# 3 is the intron start position (A-3) and 2 is the exon start position (A+2)
	acceptor_df <- extract_positional_odds_ratio_data_for_all_positions_in_window(outlier_acceptor_distances, inlier_acceptor_distances, 3,2)

	#Combine outliers and non-outliers into a compact data frame
	odds_ratio <- c(donor_df$odds_ratio, acceptor_df$odds_ratio)
	dist_to_ss <- c(donor_df$dist_to_ss, acceptor_df$dist_to_ss)
	ss_type <- c(rep("donor",length(donor_df$dist_to_ss)), rep("acceptor", length(acceptor_df$dist_to_ss)))
	lower_bounds <- c(donor_df$lower_bounds, acceptor_df$lower_bounds)
	upper_bounds <- c(donor_df$upper_bounds, acceptor_df$upper_bounds)
	df <- data.frame(odds_ratio=odds_ratio, dist_to_ss=dist_to_ss, lower_bound=lower_bounds, upper_bound=upper_bounds, ss_type=factor(ss_type, levels=c("donor","acceptor")))

	# Make seperate data frame for Donor and Acceptor splice sites
	df_acceptor <- df[as.character(df$ss_type) == "acceptor",]
	df_donor <- df[as.character(df$ss_type) == "donor", ]

	# Make error bar plot for positions around acceptor splices site
	acceptor_plot <-  ggplot() + geom_point(data=df_acceptor, mapping=aes(x=dist_to_ss, y=odds_ratio), color="#0D324D") +
					geom_errorbar(data=df_acceptor, mapping=aes(x=dist_to_ss,ymin=lower_bound, ymax=upper_bound),color="#0D324D",width=0.0) +
					labs(x = "", y = "") +
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					gtex_v8_figure_theme() +
					theme(axis.text.x=element_text(angle=45,hjust=1)) +
					scale_x_continuous(breaks=-3:1, labels=c("A-3","A-2","A-1","A+1", "A+2")) +
					theme(axis.title.y=element_blank())

	# Make error bar plot for positions around donor splice site
	donor_plot <-  ggplot() + geom_point(data=df_donor, mapping=aes(x=-dist_to_ss, y=odds_ratio), color="#0D324D") +
					geom_errorbar(data=df_donor, mapping=aes(x=-dist_to_ss,ymin=lower_bound, ymax=upper_bound),color="#0D324D",width=0.0) +
					labs(x = "", y = "") +
					geom_hline(yintercept = 1, size=.00001,linetype="dashed") +
					gtex_v8_figure_theme() + 
					theme(axis.text.x=element_text(angle=45,hjust=1)) +
					scale_x_continuous(breaks=-2:7, labels=c("D-3", "D-2", "D-1", "D+1", "D+2", "D+3","D+4", "D+5", "D+6", "D+7"))


	# Merge donor and acceptor error plots into 1 via Cowplot
	error_bar_plot <- plot_grid(donor_plot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),acceptor_plot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ncol=2, rel_widths=c(1.9,1))

	return(error_bar_plot)
}

make_panel_2e <- function(outlier_distances) {
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

	alleles <- c("A", "C", "T", "G")

	#####################################################
	# Make Outlier variant PWM plot for variants near donor site
	#####################################################
	# Positions corresponding to this plot
	ss_names <- c("D-1","D+1", "D+2", "D+3", "D+4", "D+5", "D+6")

	# Initialize PWM
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles

	# Loop through positions and fill in PWM
	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		variant_alleles <- alts[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio <= 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(variant_alleles==allele)
		}
	}
	
	# Make PWM plot with help of ggseqlogo
	outlier_donor_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier variants at donor sites") + 
	        labs(x="",y="Frequency",title="") + 
	        gtex_v8_figure_theme() +
	       	theme(plot.title = element_text(hjust = 0.5)) +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names)


	#####################################################
	# Make Background variant PWM plot for variants near donor site
	#####################################################
	# Positions corresponding to this plot
	ss_names <- c("D-1","D+1", "D+2", "D+3", "D+4", "D+5", "D+6")
	
	# Initialize PWM
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles

	# Loop through positions and fill in PWM
	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		outlier_background_alleles <- refs[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio <= 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(outlier_background_alleles==allele)
		}
	}

	# Make PWM plot with help of ggseqlogo
	background_donor_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier background at donor sites") + 
	        labs(x="",y="Frequency",title="") + 
	        theme(plot.title = element_text(hjust = 0.5)) +
	        gtex_v8_figure_theme() +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names)
	
	# Merge together background-donor-pwm plot and outlier-donor-pwm plot with cowplot
	combined_donor_plotter <- plot_grid(background_donor_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), outlier_donor_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),nrow=2)


	#####################################################
	# Make Outlier variant PWM plot for variants near acceptor site
	#####################################################
	# Positions corresponding to this plot
	ss_names <- c("A-2","A-1", "A+1")

	# Initialize PWM
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles

	# Loop through positions and fill in PWM
	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		variant_alleles <- alts[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio <= 1.0]
		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(variant_alleles==allele)
		}
	}

	# Make PWM plot with help of ggseqlogo
	outlier_acceptor_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier variants at acceptor sites") + 
	        labs(x="",y="Frequency", title="") + 
	        gtex_v8_figure_theme() +
	       	theme(plot.title = element_text(hjust = 0.5)) +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names) +
	        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(),axis.ticks.y=element_blank())


	#####################################################
	# Make Background variant PWM plot for variants near acceptor site
	#####################################################
	# Positions corresponding to this plot
	ss_names <- c("A-2","A-1", "A+1")

	# Initialize PWM
	pwm <- matrix(0, length(alleles), length(ss_names))
	rownames(pwm) = alleles

	# Loop through positions and fill in PWM
	for (ss_iter in 1:length(ss_names)) {
		ss_name <- ss_names[ss_iter]
		outlier_background_alleles <- refs[as.character(outlier_distances$distance_to_ss) == ss_name & outlier_distances$odds_ratio <= 1.0]

		for (allele_iter in 1:length(alleles)) {
			allele <- alleles[allele_iter]
			pwm[allele_iter, ss_iter] <- sum(outlier_background_alleles==allele)
		}
	}

	# Make PWM plot with help of ggseqlogo
	background_acceptor_plotter <- ggseqlogo( pwm, method = 'prob' ) +
	        #labs(x="",y="Frequency", title="Outlier background at acceptor sites") + 
	        labs(x="",y="Frequency", title="") + 
	        theme(plot.title = element_text(hjust = 0.5)) +
	        gtex_v8_figure_theme() +
	        scale_x_continuous(breaks=1:length(ss_names), labels=ss_names) +
	        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank(),axis.ticks.y=element_blank())

	
	# Merge together background-acceptor-pwm plot and outlier-acceptor-pwm plot with cowplot
	combined_acceptor_plotter <- plot_grid(background_acceptor_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), outlier_acceptor_plotter+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),nrow=2)

	# Merge together Donor and Acceptor PWM plots with cowlpot
	combined <- plot_grid(combined_donor_plotter, combined_acceptor_plotter, ncol=2, rel_widths=c(7,3)) +
				draw_text("Background for splice sites with low junction usage", x=.55, y=.97, size=8)+
				draw_text("sOutlier variants for splice sites with low junction usage", x=.55, y=.47, size=8)

	return(combined)

}

make_panel_2f <- function(fig2f_df) {
	# Simply make odds-ratio boxplots for each variant type with input data
	plotter <- ggplot(fig2f_df, aes(x=type, y=odds_ratio, fill=type)) + geom_boxplot() + 
		theme(text = element_text(size=12),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=12), legend.title = element_text(size=11)) +
		labs(x="", y="Junction usage", fill="") + 
		scale_fill_manual(values=c("dodgerblue", "violetred1")) +
		theme(legend.position="none") + 
		geom_hline(yintercept=0) +
		gtex_v8_figure_theme() +
		theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.border = element_blank())
	return(plotter)
}

##############
# Input data
fig2a_input_file = "processed_input_data/figure2/fig2a_input_data.txt"
fig2b_input_file = "processed_input_data/figure2/fig2b_input_data.txt"
fig2d_outliers_input_file = "processed_input_data/figure2/fig2d_outliers_input_data.txt"
fig2d_inliers_input_file = "processed_input_data/figure2/fig2d_inliers_input_data.txt"
fig2e_input_file = "processed_input_data/figure2/fig2e_input_data.txt"
fig2f_input_file = "processed_input_data/figure2/fig2f_input_data.txt"


##############
# Load in data
fig2a_df <- read.table(fig2a_input_file, header=TRUE, sep="\t")
fig2b_df <- read.table(fig2b_input_file, header=TRUE, sep="\t")
fig2d_outlier_df <- read.table(fig2d_outliers_input_file, header=TRUE, sep="\t")
fig2d_inlier_df <- read.table(fig2d_inliers_input_file, header=TRUE, sep="\t")
fig2e_df <- read.table(fig2e_input_file, header=TRUE, sep="\t")
fig2f_df <- read.table(fig2f_input_file, header=TRUE, sep="\t")


##############
# make panels
panel_2a <- make_panel_2a(fig2a_df)
panel_2b <- make_panel_2b(fig2b_df)
# Note that panel 2c is just a cartoon/graphic (so not included here). It was added later in Illustrator.
panel_2d <- make_panel_2d(fig2d_inlier_df, fig2d_outlier_df)
panel_2e <- make_panel_2e(fig2e_df)
panel_2f <- make_panel_2f(fig2f_df)


###############
# Combine panels with cowplot
first_row <- plot_grid(panel_2a, panel_2b, labels = c('A','B'), ncol=2)
second_row <- plot_grid(NULL, labels = c('C'), ncol=1)
third_row <- plot_grid(NULL,NULL, labels = c('D', ''),rel_heights=c(.04,1), ncol=1)
fourth_row <- plot_grid(NULL, NULL, panel_2e, NULL, labels = c("E","F", "", ""), ncol=2, rel_heights=c(.03,1), rel_widths=c(1,.56))
fig_2_tmp = plot_grid(first_row,second_row,third_row, fourth_row, nrow = 4, align='v', rel_heights=c(.98,.45,.63,.78))
fig2 <- ggdraw() + 
				draw_plot(fig_2_tmp,0,0,1,1) +
				draw_plot(panel_2f + labs(y=""),.632,-.019,.357,.255) +
				draw_plot(panel_2d,-.0145,.242,1.00,.22) +
				draw_label("Relative risk", 0.007, .39, size=8,angle=90) +
				draw_label("Junction usage", 0.655, .145, size=8,angle=90) 



###################
# Save to output file
output_file <- "generated_figures/figure2.pdf"
ggsave(fig2, file=output_file, width=7.2, height=7.0, units="in")
