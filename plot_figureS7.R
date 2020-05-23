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

data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figureS7/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'


### panel A - results of linear model of outlier statistics vs rare variant nearby ###
continuous_enrich = fread(paste0(data_dir, 'figS7A_input_data.txt'))
continuous_enrich$Category = factor(continuous_enrich$Category, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sf1 = ggplot(continuous_enrich, aes(x=Category,y=Beta)) + geom_point(size=2) +
  geom_errorbar(aes(ymin = Beta - 1.95*SE, ymax = Beta + 1.95*SE), width=0,position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0,color='grey',linetype='dashed') + xlab('') +
  theme_bw() + gtex_v8_figure_theme() +
  facet_wrap(.~Variant,scales='free') + theme(strip.background = element_blank())

### panel B - continuous rare variant enrichments split by type ###
all_coefs = fread(paste0(data_dir, 'figS7B_input_data.txt'))
plot_cols = readRDS(paste0(data_dir, 'variant_type_colors.rds'))
all_coefs$Variant = factor(all_coefs$Variant, levels=c('no_variant','other_noncoding','TE', 'coding','TSS', 'conserved_noncoding','INV','BND','DEL','CNV','DUP', 'splice_region_variant', 'splice_acceptor_variant','frameshift','splice_donor_variant', 'stop'))
all_coefs$Category = factor(all_coefs$Category, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sf2 = ggplot(all_coefs, aes(x=Variant,y=Beta,Group=Category)) +
  geom_point(size=2, aes(color=Variant,shape=Category),position=position_dodge(width=0.5)) +
  scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  geom_errorbar(aes(ymin = Beta - 1.95*SE, ymax = Beta + 1.95*SE), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Beta') + xlab('') + ylim(c(-0.5,1.5)) +
  ggtitle('') + guides(color=F) +
  scale_color_manual(values=plot_cols) +
  gtex_v8_figure_theme() +
  geom_hline(yintercept=0, linetype="dashed", colour="darkgrey") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title=element_blank(),
        legend.position=c(0.1,0.8),
        legend.key.height = unit(0.02, "cm"))

### combine
sfig <- plot_grid(sf1, sf2, nrow=2, labels=c('A', 'B'), axis='tlbr')

ggsave(sfig, file=paste0(out_dir, 'figS7.pdf'), width=7.2, height=7, units="in")


