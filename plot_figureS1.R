library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(viridis)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figureS1/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'

dcols = c('#7F5A83', '#BFCDE0', '#0D324D')
names(dcols) = c('aseOutliers', 'eOutliers', 'sOutliers')

### panel A - number of outliers per person split by self-reported population group ###
outlier_summary = fread(paste0(data_dir, 'figS1A_input_data.txt'))
outlier_summary$POP = factor(outlier_summary$POP, levels=c('African-American (n=78)', 'European-American (n=668)', 'Asian-American (n=10)', 'Other/Unknown (n=8)'))
outlier_summary$Type = factor(outlier_summary$Type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_1A = ggplot(outlier_summary, aes(x=POP,y=N,Group=Type)) + geom_boxplot(aes(fill=Type)) +
  xlab('') + ylab('# outliers/ind') + scale_fill_manual(values=dcols) +
  gtex_v8_figure_theme() + 
  theme(axis.text.x=element_text(hjust=1, angle=30),
        legend.position=c(0.9,0.8), 
        legend.title=element_blank(),
        legend.key.height = unit(0.05, "cm")) 

### panel B - number of eOutliers at a median Z-score threshold of 3 by over/under ###
medz_outliers = fread(paste0(data_dir, 'figS1B_input_data.txt'))
sfig_1B = ggplot(medz_outliers, aes(MedZ)) + geom_histogram(bins=40) + gtex_v8_figure_theme() +
  xlab('Median Z-score') + ylab('Number of eOutliers') +
  annotate("text", x=7, y=500, label='n=2215',cex=3) +
  annotate("text", x=-7, y=500, label='n=1509',cex=3)

### panel C - effect of expression data correction on rare variant enrichments ###
risk_estims = fread(paste0(data_dir, 'figS1C_input_data.txt'))
risk_estims = risk_estims %>% arrange(by=Riskratio)
risk_estims$Category = factor(risk_estims$Category, levels=unique(risk_estims$Category))
sfig_1C = ggplot(risk_estims, aes(x=Category, y=Riskratio)) + 
  geom_point(size=2, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('') +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.x=element_text(size=8)) +
  theme(legend.key.height = unit(0.05, "cm")) +
  theme(legend.title=element_blank(), legend.position=c(0.8,0.8)) + 
  facet_wrap(.~Type, scales='free') + theme(strip.background = element_blank())


## SFig 1D - enrichments across thresholds
threshold_data = fread(paste0(data_dir, 'figS1D_input_data.txt'))
threshold_data$Threshold = factor(threshold_data$Threshold, levels=c(0.01,0.001,1e-04,1e-05,1e-06,1e-07,1e-08))
threshold_data$Maf = factor(threshold_data$Maf, levels=c('novel', 'rare', 'low frequency'))
threshold_data$Method = factor(threshold_data$Method, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_1D = ggplot(threshold_data, aes(x=Threshold, y=Riskratio, Group=Method)) + 
  geom_point(size=2, aes(color=Threshold,shape=Method),position=position_dodge(width=0.5)) +
  scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('Outlier threshold') +
  scale_color_viridis(discrete=T,end=0.8) + guides(color=F) +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.x=element_text(size=8),
        legend.position=c(0.9,0.7), 
        legend.title=element_blank(),
        legend.key.height = unit(0.05, "cm")) +
  facet_wrap(.~Maf) + theme(strip.background = element_blank())

### combine
first_row <- plot_grid(sfig_1A, sfig_1B, labels = c('A', 'B'), ncol=2, align='hv', axis='tlbr')
second_row <- plot_grid(sfig_1C, sfig_1D, labels = c('C', 'D'), ncol=2, align='hv', axis='tlbr')
sfig <- plot_grid(first_row, second_row, nrow=2, rel_heights = c(1,1.4))

ggsave(sfig, file=paste0(out_dir, 'figS1.pdf'), width=7.2, height=6, units="in")





