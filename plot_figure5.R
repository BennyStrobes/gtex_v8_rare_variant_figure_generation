library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggpubr)


gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figure5/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'

dcols = c('#7F5A83', '#0D324D', '#BFCDE0')
names(dcols) = c('aseOutliers', 'sOutliers', 'eOutliers')

### panel A ###
outlier_counts = fread(paste0(data_dir, 'fig5A_input_data.txt'))
outlier_counts$variable = factor(outlier_counts$variable, levels=c('Outliers', 'Outlier RV', 'Coloc Outlier RV', 'Watershed > 0.5', 'Watershed > 0.7', 'Watershed > 0.9'))
outlier_counts$Type = factor(outlier_counts$Type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))

fig_5A = ggplot(outlier_counts, aes(x=variable,y=value+1)) + 
  geom_boxplot(aes(fill=Type),outlier.size = 0.75) + theme_bw() + # removed outlier.alpha=0.5
  xlab('') + ylab('Number per individual') +
  scale_fill_manual(values=dcols) +
  scale_y_log10() +
  gtex_v8_figure_theme() +
  theme(legend.position=c(0.1,0.85),
        legend.title=element_blank(),
        legend.key.height=unit(0.05,'in'),
        axis.text.x=element_text(hjust=1,angle=35))

### panel B ###
ukbb_any_coloc = fread(paste0(data_dir, 'fig5B_input_data.txt'))
my_comparisons <- list(c("Coloc Outlier", "Outlier"), c("Coloc Outlier", "Non-Outlier"),c("Outlier", "Non-Outlier"))
ukbb_any_coloc$Category = factor(ukbb_any_coloc$Category, levels=c('Non-Outlier', 'Outlier', 'Coloc Outlier'))
fig_5B = ggplot(ukbb_any_coloc, aes(x=Category,y=VPercentile)) + 
  geom_violin(fill='#7EB09B',alpha=0.5) + geom_boxplot(width=0.3) + theme_bw() +
  xlab('') + ylab('Variant effect size percentile') +
  stat_compare_means(method="wilcox.test", comparisons=my_comparisons, 
                     method.args=list(alternative="greater"), label="p.signif", tip.length=0) +
  gtex_v8_figure_theme() + theme(axis.text.x=element_text(hjust=1,angle=35))

### panel C ###
hit_random_table = fread(paste0(data_dir, 'fig5C_input_data.txt'))
hit_random_table$Posterior = factor(hit_random_table$Posterior, levels=c(0.01, 0.05, 0.25, 0.5, 0.75, 0.9))
fig_5C = ggplot(hit_random_table %>% filter(Category == 'Random', Model == 'Watershed'), aes(x=Posterior,y=Prop)) + geom_boxplot() +
  geom_point(data = filter(hit_random_table, Category == 'Actual', Model == 'Watershed'), aes(x=Posterior, y=Prop),color = 'red',size=1) +
  theme_bw() + xlab('Posterior threshold') + ylab('Proportion in top 25%') +
  gtex_v8_figure_theme() 

### panel D ###
asthma_pval_data = fread(paste0(data_dir, 'fig5D_pval_input_data.txt'))

asthma_pval_data$IsVar = factor(asthma_pval_data$IsVar, levels=c('VOI', 'leadSNP', 'Other'))
fig_5D1 = ggplot(asthma_pval_data, aes(x=Pos,y=-log10(pval))) + geom_point(size=0.5) + 
  theme_bw() + guides(color=F) + 
  ggtitle('Asthma') + xlab('Chromosome 9 position') +
  geom_point(data = subset(asthma_pval_data, IsVar == 'VOI'),
             aes(x = Pos, y = -log10(pval)), color = '#CC00B4', size = 1) +
  geom_point(data = subset(asthma_pval_data, IsVar == 'leadSNP'),
             aes(x = Pos, y = -log10(pval)), color = '#2699FF', size = 1) +
  geom_hline(yintercept=-log10(5e-08), color='grey') +
  gtex_v8_figure_theme()

asthma_effect_sizes = fread(paste0(data_dir, 'fig5D_betas_input_data.txt'))
asthma_effect_sizes$VOI = factor(asthma_effect_sizes$VOI, levels=c('Watershed variant','lead SNP','Other'))
fig_5D2 = ggplot(asthma_effect_sizes, aes(x=minor_AF,y=abs(beta_scaled))) +
  geom_point(aes(color=VOI),size=0.5) + theme_bw() + guides(size=F) +
  geom_point(data = subset(asthma_effect_sizes, VOI != 'Other'),
             aes(color = VOI), size=1) +
  ylab('|Effect size|') +
  scale_color_manual(values=c('#CC00B4', '#2699FF', 'black')) + xlab('UKBB MAF') +
  gtex_v8_figure_theme() + guides(color=F)

fig_5D <- plot_grid(fig_5D1, fig_5D2, labels = NULL, nrow=2)

### panel E ###
chol_pval_data = fread(paste0(data_dir, 'fig5E_pval_input_data.txt'))

fig_5E1 = ggplot(chol_pval_data, aes(x=Pos,y=-log10(pval))) + geom_point(size=0.5) + 
  theme_bw() + guides(color=F) + 
  ggtitle('High Cholesterol') + xlab('Chromosome 22 position') +
  geom_point(data = subset(chol_pval_data, IsVar == 'VOI'),
             aes(x = Pos, y = -log10(pval)), color = '#CC00B4', size = 1) +
  geom_point(data = subset(chol_pval_data, IsVar == 'leadSNP'),
             aes(x = Pos, y = -log10(pval)), color = '#2699FF', size = 1) +
  geom_hline(yintercept=-log10(5e-08), color='grey') +
  gtex_v8_figure_theme()

chol_betas_data = fread(paste0(data_dir, 'fig5E_betas_input_data.txt'))
fig_5E2 = ggplot(chol_betas_data, aes(x=minor_AF,y=abs(beta_scaled))) +
  geom_point(aes(color=VOI),size=1) + theme_bw() + guides(size=F) +
  geom_point(data = subset(chol_betas_data, VOI != 'Other'),
             aes(color = VOI), size=1) +
  ylab('|Effect size|') +
  scale_color_manual(values=c('black', '#CC00B4')) + xlab('UKBB MAF') +
  gtex_v8_figure_theme() + guides(color=F)

fig_5E <- plot_grid(fig_5E1, fig_5E2, labels = NULL, nrow=2)

## compile
first_row <- plot_grid(fig_5A, fig_5B, fig_5C, labels = c('A','B','C'), ncol=3, axis='tlbr', align='hv')
second_row1 <- plot_grid(fig_5D, labels = c('D'), ncol=1)
second_row2 <- plot_grid(fig_5E, labels = c('E'), ncol=1)
second_row <- plot_grid(second_row1, second_row2, ncol=2)

fig_5 = plot_grid(first_row, second_row, nrow = 2)

ggsave(fig_5, file=paste0(out_dir, 'fig5.pdf'), width=7.2, height=6, units='in')


