library(ggplot2)
library(data.table)
library(dplyr)
library(epitools)
library(cowplot)

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

### function to calculate promoter motif enrichments ###
calc_motif_enrich <- function(motif, data, direction) {
  motif_data = filter(data, promoter_motif == motif, medz_bin %in% c('Control', direction)) %>% arrange(by=medz_bin)
  if (nrow(motif_data) == 1) {
    return(data.frame(Risk = NA, Lower = NA, Upper = NA, Pval = NA, Motif = motif, Direction = direction))
  } else {
    counttable = rbind(c(motif_data$NumBin[1] - motif_data$NumCat[1], motif_data$NumCat[1]),
                       c(motif_data$NumBin[2] - motif_data$NumCat[2], motif_data$NumCat[2]))
    motif_risk = epitab(counttable, method='riskratio')$tab
    return(data.frame(Risk = motif_risk[2,5], Lower = motif_risk[2,6], Upper = motif_risk[2,7], Pval = motif_risk[2,8], Motif = motif, Direction = direction))
  }
}

### panel A - plot rare variant enrichments at varying windows downstream from gene ###
window_data = fread(paste0(data_dir, 'figS12A_input_data.txt'))
window_data$Method = factor(window_data$Method, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
window_data$Type = factor(window_data$Type, levels=c('SNVs', 'indels', 'SVs'))
sf1 = ggplot(window_data, aes(x=Window, y=Riskratio,Group=Method)) + 
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

### panel B - enrichment of rare variants in promoter motifs nearby outliers ###
tss_data = fread(paste0(data_dir, 'figS12B_input_data.txt')) %>% filter(promoter_motif != 'no_motif')
enrich_over = do.call(rbind, lapply(unique(tss_data$promoter_motif), function(x) calc_motif_enrich(x, tss_data, 'Over')))
enrich_under = do.call(rbind, lapply(unique(tss_data$promoter_motif), function(x) calc_motif_enrich(x, tss_data, 'Under')))
both_enrich = rbind(enrich_over, enrich_under) %>% group_by(Motif) %>%
  mutate(NumNA = length(which(is.na(Risk)))) %>% ungroup() %>%
  filter(NumNA != 2) %>% arrange(by=Risk)
both_enrich$Motif = factor(both_enrich$Motif, levels=unique(both_enrich$Motif))
sf2 = ggplot(both_enrich, aes(x=Motif, y=Risk, Group=Direction)) +
  geom_point(size=2, aes(color=Direction),position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=1, linetype='dashed', color='grey') +
  xlab('') + ylab('Relative risk') + scale_y_log10() +
  scale_color_manual(values=c('#3066BE', '#AABD8C')) +
  gtex_v8_figure_theme() + theme(legend.title=element_blank())
  
### combine
sfig <- plot_grid(sf1, sf2, labels = c('A', 'B'), nrow=2, align='hv', axis='tlbr')
ggsave(sfig, file=paste0(out_dir, 'figS12.png'), width=7.2, height=8, units="in")

