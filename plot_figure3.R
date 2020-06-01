library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(epitools)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figure3/'

dcols = c('#7F5A83', '#0D324D', '#BFCDE0')
names(dcols) = c('aseOutliers', 'sOutliers', 'eOutliers')

### panel A - outlier replication across tissues ###
sharing_med = fread(paste0(data_dir, "fig3A_input_data.txt"))
fig_order = sharing_med %>% filter(type == "aseOutliers") %>% 
  arrange(desc(frac))

sharing_med$tissue_name <- factor(sharing_med$tissue_name, levels = fig_order$tissue_name)
sharing_med$type = factor(sharing_med$type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
fig_3A = ggplot(sharing_med, aes(tissue_name, frac, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") + 
  ylab("Sharing Fraction") + xlab("") +
  scale_fill_manual(values = dcols) +
  geom_errorbar(aes(min = min, max = max), position = position_dodge(.85), width = .4) +
  gtex_v8_figure_theme() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        legend.title=element_blank(),
        legend.position=c(0.9,0.9),
        legend.key.size=unit(0.05,'in')) 

### panel B - single tissue rare variant enrichments ###
risks = fread(paste0(data_dir, 'fig3B_input_data.txt'))
outlier_threshold <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9)
risks$sig <- factor(risks$sig, levels = outlier_threshold)

risks$outlier_type = factor(risks$outlier_type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
fig_3B = ggplot(risks, aes(sig, ratio, fill = outlier_type)) +
  geom_boxplot(size = .3, outlier.size = .3)  +
  scale_fill_manual(values = dcols) + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("Outlier p-value") + ylab("Relative risk") +
  gtex_v8_figure_theme() + guides(fill=F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### panel C - single tissue enrichments by variant type ###
type_risks = fread(paste0(data_dir, "fig3C_input_data.txt"))
type_risks$outlier_type = factor(type_risks$outlier_type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
type_risks$Variant_type = factor(type_risks$Variant_type, levels=c('Splice donor', 'Splice acceptor', 'Stop gained', 'Frameshift', 'Splice', 'Inframe deletion'))
fig_3C = ggplot(type_risks, aes(x = Variant_type, y = ratio, fill = outlier_type)) +
  geom_boxplot(size = .3, outlier.size = .3) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), 
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = dcols) +
  ylab("Relative risk") + guides(fill=F) +
  scale_y_continuous(trans='log2') +
  gtex_v8_figure_theme() 

### panel D - distribution of proportion of tissues with expression changes in median Z-score and correlation outliers ###
mcols = c('#C1CDDE', '#7FCDBB')
names(mcols) = c('MEDZ', 'Correlation')
mfills = c(alpha('#C1CDDE', 0.4), alpha('#7FCDBB',0))
names(mfills) = names(mcols)

outlier_props = fread(paste0(data_dir, 'fig3D_input_data.txt'))

fig_3D = ggplot(outlier_props, aes(x = prop, colour = Method, fill = Method)) +
  geom_density(bw = 0.04, size = 0.75) + ylab('Density') +
  xlab('Proportion tissues with |Z| >= 3') +
  scale_color_manual(values = mcols, breaks = names(mcols)) +
  scale_fill_manual(values=mfills) +
  guides(fill = FALSE, colour = guide_legend(override.aes = list(fill = mfills))) +
  theme(panel.border = element_blank()) +
  gtex_v8_figure_theme() +
  theme(legend.position = c(0.5, 0.85),
        legend.direction = 'vertical',
        legend.title=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size=unit(0.05,'in'))

### panel E - rare variant enrichment in tissue-specific enhancers ###
enh_risks = fread(paste0(data_dir, 'fig3E_input_data.txt'))
enh_risks$Zscore = factor(enh_risks$ZT, levels=c(3,5,7))
enh_risks$Type = factor(enh_risks$Type, levels=c('Matched', 'Unmatched'))
ccols = c('#C1E0BC', '#7FCDBB', '#00ABB8')
fig_3E = ggplot(enh_risks, aes(x=Type,y=Riskratio,Group=Zscore)) +
  geom_point(size=2.5, position=position_dodge(width=0.5),aes(color=Zscore)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0, position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('Enhancer type') +
  labs(color = "Tissue |Z|") +
  scale_color_manual(values=ccols) +
  ggtitle('') + geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  theme(panel.border = element_blank()) +
  gtex_v8_figure_theme() +
  theme(legend.position = c(0.8, 0.85),
        legend.direction = 'vertical',
        legend.key.size=unit(0.05,'in'))

### combine
first_column <- plot_grid(fig_3A, labels = c('A'), nrow=1)
second_column <- plot_grid(fig_3B, fig_3C, fig_3D, fig_3E, labels = c('B', 'C', 'D', 'E'), ncol=4, axis='tblr', align='h')
fig_3 = plot_grid(first_column, second_column, nrow = 2, align='hv', rel_heights = c(1.5,1))
ggsave(fig_3, file=paste0(out_dir, 'fig3.pdf'), width=7.2, height=7.2,units="in")
