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


data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figureS41/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'

### summary stats across 4 cholesterol related traits in MVP ###
chol_stats = fread(paste0(data_dir, 'figS41_input_data.txt'))
sfig = ggplot(all_chol_data, aes(x=GAF,y=abs(BETA))) +
  geom_point(aes(color=VOI),size=0.4) + theme_bw() + 
  geom_point(data = subset(all_chol_data, VOI == 1),
             aes(color = VOI),size=2) +
  scale_color_manual(values=c('black', '#CC00B4')) +
  xlab('gnomAD minor allele frequency') + guides(color=F,size=F) +
  gtex_v8_figure_theme() + facet_wrap('Trait', ncol=2) +
  ylab('|Effect size|') + 
  theme(strip.background = element_blank(),
        strip.text.x=element_text(size=8))

ggsave(sfig, file=paste0(out_dir, 'figS41.pdf'), width=7.2, height=4, units="in")

