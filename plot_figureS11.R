library(ggplot2)
library(ggthemes)
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


data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figureS11/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'

dcols = c('#7F5A83', '#0D324D', '#BFCDE0')
names(dcols) = c('aseOutliers', 'sOutliers', 'eOutliers')

### panel A - comparison of relative risk of nearby rare variants between outlier types ###
risk_compare = fread(paste0(data_dir, 'figS11A_input_data.txt'))
plot_cols = readRDS(paste0(data_dir, 'variant_type_colors.rds'))
risk_compare$Type = factor(risk_compare$Type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_A = ggplot(risk_compare %>% filter(is.finite(FC)), aes(x=Cat, y=FC)) + 
  geom_point(aes(color=Cat,shape=Type),size=2) + theme_bw() + xlab('') +
  scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  ylab('Enrichment / max enrichment') +
  scale_color_manual(values=plot_cols) + 
  guides(color=F) + geom_hline(yintercept=1,color='grey',linetype='dashed') +
  gtex_v8_figure_theme() + theme(axis.text.x=element_text(hjust=1,angle=45),
                                 legend.title=element_blank(),
                                 legend.position=c(0.9,0.2))


### panel B - proportion of rare variants leading to outlier signal across types ###
prop_variants = fread(paste0(data_dir, 'figS11B_input_data.txt'))
prop_variants$variable = factor(prop_variants$variable, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_B = ggplot(prop_variants,aes(x=Feature,y=Risk,Group=variable)) +
  geom_bar(aes(fill=variable),stat='identity',color='black',position='dodge') +
  scale_y_sqrt() + scale_fill_manual(values=dcols) + theme_bw() + 
  xlab('') + ylab('Proportion variants') + 
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position=c(0.2,0.8),
        legend.title=element_blank(),
        legend.key.height = unit(0.1, "cm"))

## combine
sfig = plot_grid(sfig_A, sfig_B, nrow = 2, labels=c('A', 'B'), align='v')

ggsave(sfig, file=paste0(out_dir, 'figS11.pdf'), width=7.2, height=7.2,units="in")

