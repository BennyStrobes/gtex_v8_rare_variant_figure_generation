## check variants in TOPMed

library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figureS13/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'

### panel A - observed and expected outlier pairs within a window ###
pair_counts = fread(paste0(data_dir, 'figS13A_input_data.txt'))
pair_counts$Type = factor(pair_counts$Type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
pair_counts$variable = factor(pair_counts$variable, levels=c('Observed', 'Expected'))
sfig_A = ggplot(pair_counts, aes(x=Window, y=value, Group=variable)) +
  geom_bar(stat='identity',aes(fill=variable),position='dodge') + 
  gtex_v8_figure_theme() + facet_wrap(.~Type) + xlab('Window size: log10(bp)') +
  ylab('Outlier pairs within window') +
  theme(legend.title=element_blank(),
        strip.background = element_blank(),
        legend.position=c(0.1,0.8),
        legend.key.height = unit(0.05, "cm"),
        panel.spacing = unit(2, "lines"))

### panel B - enrichment of outlier pairs across window sizes ###
pair_enrich = fread(paste0(data_dir, 'figS13B_input_data.txt'))
pair_enrich$Type = factor(pair_enrich$Type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_B = ggplot(pair_enrich, aes(x=Window, y=LogRatio)) +
  geom_bar(stat='identity') + geom_hline(yintercept=0, color='lightgrey') +
  gtex_v8_figure_theme() + facet_wrap(.~Type) + xlab('Window size: log10(bp)') +
  ylab('Expected / Observed') +
  theme(strip.background = element_blank(),
        panel.spacing = unit(2, "lines"))

### panel C - sOutlier enrichments w/o filtering ###
sp_enrich = fread(paste0(data_dir, 'figS13C_input_data.txt'))
sfig_C = ggplot(sp_enrich, aes(x=Window, y=LogRatio)) +
  geom_bar(stat='identity') + geom_hline(yintercept=0, color='lightgrey') +
  ggtitle('sOutlier enrichments w/o filtering') +
  gtex_v8_figure_theme() + xlab('Window size: log10(bp)') +
  ylab('Expected / Observed') +
  theme(strip.background = element_blank())

### panel D - enrichment of rare variants nearby eOutlier pairs ###
eoutlier_enrich = fread(paste0(data_dir, 'figS13D_input_data.txt'))
eoutlier_enrich$Window = factor(eoutlier_enrich$Window, levels=c('100kbp', '500kbp', '1000kbp', '5000kbp'))
sfig_D = ggplot(eoutlier_enrich, aes(x=Cat, y=Riskratio)) + 
  geom_point(size=3) + xlab('') + ylab('Relative risk') +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0) +
  gtex_v8_figure_theme() + scale_y_log10() + 
  ggtitle('eOutlier pairs') +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey")

### combine
third_row <- plot_grid(sfig_C, sfig_D, labels=c('C', 'D'), ncol=2, align='hv', axis='tblr')
sfig <- plot_grid(sfig_A, sfig_B, third_row, labels=c('A', 'B', NULL, NULL), nrow=3)

ggsave(sfig, file=paste0(out_dir, 'figS13.pdf'), width=7.2, height=9, units="in")



