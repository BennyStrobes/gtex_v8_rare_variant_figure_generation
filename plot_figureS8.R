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

data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figureS8/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'

dcols = c('#7F5A83', '#BFCDE0', '#0D324D')
names(dcols) = c('aseOutliers', 'eOutliers', 'sOutliers')

### panel A - distribution of proportion of tested tissues supporting multi-tissue outlier calls ###
outlier_summary = fread(paste0(data_dir, 'figS8AB_input_data.txt'))
outlier_summary$Category = factor(outlier_summary$Category, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))

sf1 = ggplot(outlier_summary, aes(Prop)) + geom_density(aes(fill=Category),alpha=0.6) + 
  theme_bw() + scale_fill_manual(values=dcols) + 
  xlab('Proportion tissues supporting outlier call') +
  gtex_v8_figure_theme() + guides(fill=F) +
  theme(legend.title=element_blank(),
        legend.position=c(0.2,0.9)) 

### panel B - distribution of number of tissues supporting multi-tissue outlier calls ###
sf2 = ggplot(outlier_summary, aes(NT)) + geom_density(aes(fill=Category),alpha=0.6) + 
  theme_bw() + scale_fill_manual(values=dcols) + 
  xlab('Number tissues supporting outlier call') +
  gtex_v8_figure_theme() + 
  theme(legend.title=element_blank(),
        legend.position=c(0.8,0.9)) 

### panel C - rare variant enrichments depending on tissue support ###
tissue_enrich = fread(paste0(data_dir, 'figS8C_input_data.txt'))
tissue_enrich$Threshold = factor(tissue_enrich$Threshold, levels=c('>1','>2','25','50','75','90','100'))
tissue_enrich$Category = factor(tissue_enrich$Category, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sf3 = ggplot(tissue_enrich %>% filter(!is.na(Upper)), aes(x=Threshold, y=Riskratio, Group=Category)) +
  geom_point(size=3, aes(color=Category),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  scale_color_manual(values=dcols) + xlab('Percentage tissues supporting outlier call') +
  geom_hline(yintercept=1, color='grey', linetype='dashed') + ylab('Relative risk') +
  theme_bw() + gtex_v8_figure_theme() + guides(color=F) +
  facet_wrap(.~Feature,scales='free') + theme(strip.background = element_blank())

### combine
first_row <- plot_grid(sf1, sf2, labels=c('A','B'), ncol=2, align='hv', axis='tlbr')
second_row <- plot_grid(sf3, labels=c('C'), ncol=1)
sfig <- plot_grid(first_row, second_row, nrow=2, align='hv', axis='tlbr')

ggsave(sfig, file=paste0(out_dir, 'figS8.pdf'), width=7.2, height=8, units="in")
