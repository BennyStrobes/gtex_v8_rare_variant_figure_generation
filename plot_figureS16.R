library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


data_dir = 'processed_input_data/figureS16/'
out_dir = 'generated_figures/'

### panel A - Z-score distribution across tissues for individual with fusion transcript ###
zscore_data = fread(paste0(data_dir, 'figS16a_input_data.txt'))
vcols = c('#3772FF', '#E2EF70')
names(vcols) = c('EML6', 'SPTBN1')
sfig_A = ggplot(zscore_data, aes(Z)) + geom_density(aes(fill=Category),alpha=0.8) + theme_bw() +
  gtex_v8_figure_theme() + xlab('Tissue Z-score') +
  scale_fill_manual(values=vcols) +
  theme(legend.title=element_blank(),
        legend.position=c(0.85,0.9))

### combined with Sashimi plots from IGV (panels B+C) in Illustrator
ggsave(sfig_A, file=paste0(out_dir, 'figS16.pdf'), width=7.2, height=3, units='in')
