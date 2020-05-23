### analysis for supplemental figure for CADD

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


data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figureS40/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'

### panel A - distribution of watershed posteriors ###
variant_data = fread(paste0(data_dir, 'figS40AB_input_data.txt'))
sfig_A = ggplot(variant_data, aes(MaxWS)) + geom_histogram() + 
  xlab('Watershed posterior') + ylab('Count') + gtex_v8_figure_theme()

### panel B - distribution of CADD scores ###
sfig_B = ggplot(variant_data, aes(RawScore)) + geom_histogram() + 
  xlab('CADD score') + ylab('Count') + gtex_v8_figure_theme()

### panel C - comparison of Watershed and CADD scores ###
compare_data = fread(paste0(data_dir, 'figS40C_input_data.txt'))
sfig_C = ggplot(compare_data, aes(x=RawScore,y=ws_posterior)) + geom_point() +
  xlab('CADD Score') + ylab('Watershed posterior') +
  geom_vline(xintercept=2.3,color='blue') + 
  geom_hline(yintercept=0.5,color='blue') +
  gtex_v8_figure_theme()

### panel D - enrichment of variant types in high CADD vs high Watershed variant lists ###
variant_enrich = fread(paste0(data_dir, 'figS40D_input_data.txt'))
sfig_D = ggplot(variant_enrich, aes(x=Name, y=FC)) + geom_bar(stat='identity') + 
  scale_y_log10() + theme_bw() + xlab('') +
  ylab('Proportion watershed / Proportion CADD') +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(hjust=1,angle=45))

### panel E - proportion of variants at given CADD thresholds in top 25% of variant effect sizes in colocalized traits ###
hit_cadd_random_table = fread(paste0(data_dir, 'figS40E_input_data.txt'))
hit_cadd_random_table$PT = factor(hit_cadd_random_table$PT, levels=c(0.01,0.05,0.25,0.5,0.75,0.9))
sfig_E = ggplot(hit_cadd_random_table %>% filter(Cat == 'Random', Model == 'CADD'), aes(x=PT,y=Prop)) + geom_boxplot() +
  geom_point(data = filter(hit_cadd_random_table, Cat == 'Actual', Model == 'CADD'), aes(x=PT, y=Prop),color = 'red',size=2) +
  theme_bw() + xlab('Posterior threshold') + ylab('Proportion of variants in top 25%') +
  gtex_v8_figure_theme()

### combine
first_row <- plot_grid(sfig_A, sfig_B, sfig_C, ncol=3,labels=c('A', 'B', 'C'), align='hv', axis='tlbr')
second_row <- plot_grid(sfig_D, sfig_E, ncol=2, labels=c('D', 'E'), align='hv', axis='tlbr')
sfig <- plot_grid(first_row, second_row, nrow=2)

ggsave(sfig, file=paste0(out_dir, 'figS40.pdf'), width=7.2, height=8, units="in")
