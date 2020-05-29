library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figureS12/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'

dcols = c('#7F5A83', '#BFCDE0', '#0D324D')
names(dcols) = c('aseOutliers', 'eOutliers', 'sOutliers')

### plot rare variant enrichments at varying windows downstream from gene ###
window_data = fread(paste0(data_dir, 'figS12_input_data.txt'))
window_data$Method = factor(window_data$Method, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
window_data$Type = factor(window_data$Type, levels=c('SNVs', 'indels', 'SVs'))
sf = ggplot(window_data, aes(x=Window, y=Riskratio,Group=Method)) + 
  geom_point(size=3,aes(color=Method), position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('Distance downstream from gene (kb)') + 
  geom_line(aes(group=Method),size=0.5) +
  geom_hline(yintercept=1,color='grey') +
  scale_color_manual(values=dcols) + 
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.x=element_text(size=8),
        legend.title=element_blank()) +
  facet_wrap(Type~.,scales='free', ncol=3) + theme(strip.background = element_blank())

ggsave(sf, file=paste0(out_dir, 'figS12.pdf'), width=7.2, height=4, units="in")

