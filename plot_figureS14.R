library(ggplot2)
library(ggthemes)
library(data.table)
library(dplyr)
library(cowplot)
library(epitools)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


data_dir = 'processed_input_data/figureS14/'
out_dir = 'generated_figures/'

calc_risk <- function(data, zt) {
  data = filter(data, Z == zt)
  counttable = rbind(as.numeric(c(data[1,1], data[1,2])), as.numeric(c(data[1,3], data[1,4])))
  rr = epitab(counttable,method="riskratio")$tab
  estims = data.frame(Riskratio = rr[2,5], Lower = rr[2,6], Upper = rr[2,7],
                      Pval = rr[2, 8], Threshold = zt)
  return(estims)
}

### panel A - enrichment of rare SVs near eOutliers split by gene-overlapping and not ###
sv_count_table = fread(paste0(data_dir, 'figS14a_input_data.txt'))
coding_enrich = do.call(rbind, lapply(list(3,4,5,6), function(x) calc_risk(sv_count_table %>% filter(Category == 'Coding'), x)))
noncoding_enrich = do.call(rbind, lapply(list(3,4,5,6), function(x) calc_risk(sv_count_table %>% filter(Category == 'Noncoding'), x)))

both_enrich = rbind(coding_enrich %>% mutate(Category = 'Overlapping gene'),
                    noncoding_enrich %>% mutate(Category = 'Non-coding within 10kb'))
both_enrich$Threshold = factor(both_enrich$Threshold, levels=c(3,4,5,6))
sfig_A = ggplot(both_enrich, aes(x=Threshold,y=Riskratio,Group=Category)) +
  geom_point(size=1.5, aes(color=Category),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('eOutlier |Median Z-score| threshold') +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() + 
  theme(legend.title=element_blank(), legend.position=c(0.325,0.95)) +
  theme(legend.key.height = unit(0.05, "cm"),
        legend.key.width = unit(0.04, "cm"))

## panel B - number outlier-associated rare SVs per person (over and under) ###
count_data = fread(paste0(data_dir, 'figS14b_input_data.txt'))
mean_over = round(mean(filter(count_data, Direction == 'Over eOutliers')$NumSV),2)
mean_under = round(mean(filter(count_data, Direction == 'Under eOutliers')$NumSV),2)
sfig_B = ggplot(count_data, aes(x=Direction, y=NumSV)) + geom_jitter() + geom_violin() + theme_bw() +
  ylab('# outlier-associated rare SVs per person') + xlab('Direction') +
  annotate("text", x=1, y=6, label=paste0('Mean = ', mean_over), cex=3) +
  annotate("text", x=2, y=7.8, label=paste0('Mean = ', mean_under), cex=3) +
  gtex_v8_figure_theme()

## panel C - median Z-scores of gene pairs affected by same rare SV ###
sv_gene_data = fread(paste0(data_dir, 'figS14c_input_data.txt'))
vcols = c('#FCBBA1', '#EA5948', '#9E0142')
names(vcols) = c('BND', 'DEL', 'DUP')
sfig_C = ggplot(sv_gene_data, aes(x=MedZ.x,y=MedZ.y)) + geom_point(size=1.5,aes(color=variant_cat.x)) + theme_bw() +
  xlab('|Median Z| of Gene 1') + ylab('|Median Z| of Gene 2') +
  gtex_v8_figure_theme() + scale_color_manual(values=vcols) +
  geom_hline(yintercept=3, color='grey', linetype='dashed') +
  geom_hline(yintercept=-3, color='grey', linetype='dashed') +
  geom_vline(xintercept=3, color='grey', linetype='dashed') +
  geom_vline(xintercept=-3, color='grey', linetype='dashed') +
  theme(legend.title=element_blank()) +
  theme(legend.key.height = unit(0.05, "cm"),
        legend.position=c(0.9,0.5))

### combine
sfig = plot_grid(sfig_A, sfig_B, sfig_C, labels = c('A', 'B', 'C'), nrow = 1, axis='tlbr')

ggsave(sfig, file=paste0(out_dir, 'figS14.pdf'), width=7.2, height=3,units="in")
