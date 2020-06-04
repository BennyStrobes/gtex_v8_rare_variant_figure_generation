library(ggplot2)
library(data.table)
library(dplyr)
library(cowplot)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = 'processed_input_data/figureS25/'
out_dir = 'generated_figures/'

### panel A - plot enrichment of nearby rare SNVs in single tissue eOutliers across thresholds ###
single_enrich = fread(paste0(data_dir, 'figS25_input_data.txt'))
single_enrich$ZT = factor(single_enrich$ZT, levels=unique(single_enrich$ZT))
single_enrich$Type = factor(single_enrich$Type, levels=c('SNVs', 'indels', 'SVs'))

sfig_A = ggplot(single_enrich %>% filter(Type == 'SNVs'), aes(x=ZT, y=Riskratio)) + geom_boxplot() + 
  xlab('|Z-score| threshold') + ylab('Relative risk') + ggtitle('SNVs') +
  gtex_v8_figure_theme() + geom_hline(yintercept=1, color='darkgrey')

### panel B - plot enrichment of nearby rare indels in single tissue eOutliers across thresholds ###
sfig_B = ggplot(single_enrich %>% filter(Type == 'indels'), aes(x=ZT, y=Riskratio)) + geom_boxplot() + 
  xlab('|Z-score| threshold') + ylab('') + ggtitle('indels') +
  gtex_v8_figure_theme() + geom_hline(yintercept=1, color='darkgrey')

### panel C - plot enrichment of nearby rare SVs in single tissue eOutliers across thresholds ###
sfig_C = ggplot(single_enrich %>% filter(Type == 'SVs'), aes(x=ZT, y=Riskratio)) + geom_boxplot() + 
  xlab('|Z-score| threshold') + ylab('') + ggtitle('SVs') +
  gtex_v8_figure_theme() + geom_hline(yintercept=1, color='darkgrey')

### combine
sfig <- plot_grid(sfig_A, sfig_B, sfig_C, labels = c('A','B', 'C'), ncol=3)

ggsave(sfig, file=paste0(out_dir, 'figS25.pdf'), width=7.2, height=4,units="in")

