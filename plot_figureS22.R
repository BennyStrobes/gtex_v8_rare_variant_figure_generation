library(gplots)
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


data_dir = 'processed_input_data/figureS22/'
out_dir = 'generated_figures/'

dcols = c('#7F5A83', '#BFCDE0', '#0D324D')
names(dcols) = c('aseOutliers', 'eOutliers', 'sOutliers')

### proportion of outliers from each tissue that replicate in clinically accessible tissues ###
all_rep_rates = fread(paste0(data_dir, 'figS22_input_data.txt'))
all_rep_rates$variable = factor(all_rep_rates$variable, levels=c('Whole Blood', 'Fibroblasts', 'LCLs', '>1 Tissue'))
all_rep_rates$Type = factor(all_rep_rates$Type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))

sfig = ggplot(all_rep_rates, aes(x=long, y=value)) +
  geom_bar(aes(fill=Type),stat='identity',position='dodge') +
  theme_bw() + gtex_v8_figure_theme() + ylab('Replication proportion') +
  scale_fill_manual(values=dcols) +
  theme(legend.title=element_blank(),
        legend.position=c(0.8,0.85),
        legend.key.height=unit(0.05,'in')) +
  xlab('') + coord_flip() + facet_wrap(.~variable) + theme(strip.background = element_blank())

ggsave(sfig, file=paste0(out_dir, 'figS22.pdf'), width=7.2, height=8, units='in')
