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



data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figure1/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'

### panel B - relative risk of rare variant nearby outliers at varying allele frequencies ###
counts = fread(paste0(data_dir, 'fig1B_input_data.txt'))
risks = do.call(rbind, lapply(1:nrow(counts), function(x) 
  epitab(rbind(c(counts$Control_0[x], counts$Control_1[x]), c(counts$Outlier_0[x], counts$Outlier_1[x])),method="riskratio")$tab[2,]))
risks = cbind(counts, risks[,5:8])

risks$MAF = factor(risks$MAF, levels=c('novel', 'single', 'double', 'rare', 'low frequency'))
risks$Type = factor(risks$Type, levels=c('SVs', 'SNVs+indels'))
risks$Method = factor(risks$Method, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
fig_1B = ggplot(risks, aes(x=MAF,y=riskratio,Group=Method)) +
  geom_point(size=3, aes(shape=Method),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('') +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() + scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.x=element_text(size=8)) +
  theme(legend.title=element_blank(), legend.position=c(0.8,0.9)) +
  theme(legend.key.height = unit(0.05, "cm")) +
  facet_wrap(.~Type,scales='free',nrow=2) + theme(strip.background = element_blank())

### panel C - relative risk of rare variant nearby outliers across varying variant annotations ###
plot_cols = readRDS(paste0(data_dir, 'variant_type_colors.rds'))
type_counts = fread(paste0(data_dir, 'fig1C_input_data.txt'))
type_risks = do.call(rbind, lapply(1:nrow(type_counts), function(x) 
  epitab(rbind(c(type_counts$Control_0[x], type_counts$Control_1[x]), c(type_counts$Outlier_0[x], type_counts$Outlier_1[x])),method="riskratio")$tab[2,]))
type_risks = cbind(type_counts, type_risks[,5:8])

type_risks$Variant_category = factor(type_risks$Variant_category, levels=c('no_variant','other_noncoding','TE', 'coding','TSS', 'conserved_noncoding','INV','BND','DEL','CNV','DUP', 'splice_region_variant', 'splice_acceptor_variant','frameshift','splice_donor_variant', 'stop'))
type_risks$Type = factor(type_risks$Type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))

fig_1Ca = ggplot(type_risks, aes(x=Variant_category,y=riskratio,Group=Type)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0,position=position_dodge(width=0.5)) +
  geom_point(size=2, aes(shape=Type,color=Variant_category),position=position_dodge(width=0.5)) +
  theme_bw() + xlab('') + ylab('Relative risk') + ylim(c(0,300)) +
  ggtitle('') + guides(color=F,shape=F) +
  scale_color_manual(values=plot_cols) +
  scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title=element_blank())

cat_keep = as.data.frame(table(filter(type_risks, riskratio < 40)$Variant_category)) %>% filter(Freq == 3)
fig_1Cb = ggplot(type_risks %>% filter(Variant_category %in% cat_keep$Var1), aes(x=Variant_category,y=riskratio,Group=Type)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0,position=position_dodge(width=0.5)) +
  geom_point(size=1, aes(shape=Type,color=Variant_category),position=position_dodge(width=0.5)) +
  theme_bw() + ylab('') + xlab('') + scale_y_continuous(trans='log2') +
  ggtitle('') + guides(color=F,shape=F) +
  scale_color_manual(values=plot_cols) +
  scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA))

fig_1C = fig_1Ca + annotation_custom(ggplotGrob(fig_1Cb), xmin = 0.1, xmax = 12, 
                                     ymin = 30, ymax = 350)


### panel D - proportion of outliers at varying thresholds with nearby rare variant across variant annotations ###
variant_props = fread(paste0(data_dir, 'fig1D_input_data.txt'))
variant_props$Variant_category = factor(variant_props$Variant_category, levels=c('no_variant','other_noncoding','coding', 'conserved_noncoding', 'TSS', 'stop', 'frameshift', 'splice_region_variant', 'splice_donor_variant', 'splice_acceptor_variant', 'TE', 'INV','BND','DEL','CNV','DUP'))
variant_props$pval_bin = factor(variant_props$pval_bin, levels=c('0~1e-07', '1e-07~1e-05', '1e-05~1e-04', '1e-04~1e-03', '1e-03~1e-02', '1e-02~5e-02', 'nonOutlier'))
variant_props$Type = factor(variant_props$Type, levels=c('eOutliers - Over', 'eOutliers - Under', 'aseOutliers', 'sOutliers'))
plot_cols['no_variant'] = 'white'

fig_1D = ggplot(variant_props, aes(x=pval_bin,y=Proportion,Group=Variant_category)) +
  geom_bar(aes(fill=Variant_category),color='black',stat='identity') + 
  theme_bw() + coord_flip() + ylab('') + xlab('') +
  scale_fill_manual(values=plot_cols) + guides(shape=F) +
  ggtitle('') + xlab('P-value bin') + ylab('Proportion with variant') +
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text.x=element_text(size=8)) +
  gtex_v8_figure_theme() + theme(legend.title=element_blank()) +
  theme(legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.05, 'cm'),
        legend.position = "bottom") +
  facet_wrap(.~Type,nrow=4) + theme(strip.background = element_blank())


### panel E - proportion of rare variants of different annotations within 10kb of gene associated with outlier signal ###
dcols = c('#7F5A83', '#0D324D', '#BFCDE0')
names(dcols) = c('aseOutliers', 'sOutliers', 'eOutliers')
abs_data = fread(paste0(data_dir, 'fig1E_input_data.txt'))
abs_data$variable = factor(abs_data$variable, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
fig_1E = ggplot(abs_data,aes(x=Feature,y=Proportion,Group=variable)) +
  geom_bar(aes(fill=variable),stat='identity',color='black',position='dodge') +
  scale_fill_manual(values=dcols) + theme_bw() + 
  xlab('') + ylab('Proportion of variants leading to outlier') + 
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position=c(0.25,0.8),
        legend.title=element_blank(),
        legend.key.height = unit(0.1, "cm"))

## Panel 1A is overview diagram, added in Illustrator
first_col <- plot_grid(NULL, fig_1C, fig_1E, labels = c('A','C', 'E'), nrow=3, align='hv')
second_col <- plot_grid(fig_1B, fig_1D, labels=c('B', 'D'), nrow=2, align='hv', rel_heights=c(1,2))
fig_1 = plot_grid(first_col, second_col, ncol = 2, align='v')
ggsave(fig_1, file=paste0(out_dir, 'fig1.pdf'), width=7.2, height=11, units="in")





