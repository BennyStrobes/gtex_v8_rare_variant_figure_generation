library(data.table)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}



data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figureS10/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'

### panel A - GO enrichment for non-outlier genes ###
no_outlier_GO = fread(paste0(data_dir, 'figS10A_input_data.txt'))
ccols = c('grey', alpha('maroon4',0.7))
names(ccols) = c('no', 'yes')

te_GO = no_outlier_GO %>% filter(Type == 'eOutliers')
te_GO$GO_term = factor(te_GO$GO_term, levels=unique(te_GO$GO_term))
te_plot = ggplot(te_GO, aes(x=GO_term,y=-log10(FDR))) + 
  geom_bar(stat='identity',aes(fill=Significant),width=0.5) + 
  theme_bw() + coord_flip() + xlab('') + 
  scale_fill_manual(values=ccols) + guides(fill=F) +
  ggtitle('eOutlier (n=261)') + gtex_v8_figure_theme()

ase_GO = no_outlier_GO %>% filter(Type == 'aseOutliers')
ase_GO$GO_term = factor(ase_GO$GO_term, levels=unique(ase_GO$GO_term))
ase_plot = ggplot(ase_GO, aes(x=GO_term,y=-log10(FDR))) + 
  geom_bar(stat='identity',aes(fill=Significant),width=0.5) + 
  theme_bw() + coord_flip() + xlab('') + 
  scale_fill_manual(values=ccols) + guides(fill=F) +
  ggtitle('aseOutlier (n=11573)') + gtex_v8_figure_theme()

sp_GO = no_outlier_GO %>% filter(Type == 'sOutliers')
sp_GO$GO_term = factor(sp_GO$GO_term, levels=unique(sp_GO$GO_term))
sp_plot = ggplot(sp_GO, aes(x=GO_term,y=-log10(FDR))) + 
  geom_bar(stat='identity',aes(fill=Significant),width=0.5) + 
  theme_bw() + coord_flip() + xlab('') + 
  scale_fill_manual(values=ccols) + guides(fill=F) +
  ggtitle('sOutlier (n=9548)') + gtex_v8_figure_theme()

sfigA = plot_grid(te_plot, ase_plot, sp_plot, nrow = 1, axis='tb')

### panel B - GO enrichment for extreme outlier genes ###
ex_outlier_GO = fread(paste0(data_dir, 'figS10B_input_data.txt'))

te_ex_GO = ex_outlier_GO %>% filter(Type == 'eOutliers')
te_ex_GO$GO_term = factor(te_ex_GO$GO_term, levels=unique(te_ex_GO$GO_term))

te_ex_plot = ggplot(te_ex_GO, aes(x=GO_term,y=-log10(Pval))) + 
  geom_bar(stat='identity',aes(fill=Significant),width=0.5) + theme_bw() + coord_flip() + xlab('') + 
  scale_fill_manual(values=ccols) + guides(fill=F) +
  ggtitle('eOutlier (n=127)') + gtex_v8_figure_theme()

ase_ex_GO = ex_outlier_GO %>% filter(Type == 'aseOutliers')
ase_ex_GO$GO_term = factor(ase_ex_GO$GO_term, levels=unique(ase_ex_GO$GO_term))

ase_ex_plot = ggplot(ase_ex_GO, aes(x=GO_term,y=-log10(Pval))) + 
  geom_bar(stat='identity',aes(fill=Significant),width=0.5) + theme_bw() + coord_flip() + xlab('') + 
  scale_fill_manual(values=ccols) + guides(fill=F) +
  ggtitle('aseOutlier (n=261)') + gtex_v8_figure_theme()

sp_ex_GO = ex_outlier_GO %>% filter(Type == 'sOutliers')
sp_ex_GO$GO_term = factor(sp_ex_GO$GO_term, levels=unique(sp_ex_GO$GO_term))

ase_ex_plot = ggplot(sp_ex_GO, aes(x=GO_term,y=-log10(Pval))) + 
  geom_bar(stat='identity',aes(fill=Significant),width=0.5) + theme_bw() + coord_flip() + xlab('') + 
  scale_fill_manual(values=ccols) + guides(fill=F) +
  ggtitle('sOutlier (n=389)') + gtex_v8_figure_theme()

sfigB = plot_grid(te_ex_plot, ase_ex_plot, sp_ex_plot, nrow = 1, axis='tb')

### combine
sfig = plot_grid(sfigA, sfigB, labels=c('A', 'B'), nrow = 2, axis='tb', align='hv')
ggsave(sfig, file=paste0(out_dir, 'figS10.pdf'), width=7.2, height=4, units="in")




