library(gplots)
library(ggplot2)
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


data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figureS26/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'

### panel A - comparison of imputation methods ###
recon_errors = fread(paste0(data_dir, 'figS26A_input_data.txt'))
meth.colors = c(brewer.pal(9, 'YlGnBu')[3:9],'grey')
names(meth.colors) = c('EM','KNN','MEAN','PMD','SOFT','STFZ','MEDZ','RANDOM')
recon_errors$MethodName = factor(recon_errors$MethodName,levels=c('EM', 'KNN', 'MEAN', 'PMD', 'SOFT', 'RANDOM'))
sfig_A = ggplot(recon_errors,aes(x=MethodName,y=Error,Group=MethodName)) + 
  geom_boxplot(aes(fill=MethodName),outlier.size=0.5) + theme_bw() +
  scale_fill_manual(values=meth.colors) + xlab('Imputation method') +
  guides(fill=F) + gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(hjust=1,angle=40))

### panel B - Comparing values of k for KNN imputation ###
recon_error_knn = fread(paste0(data_dir, 'figS26B_input_data.txt'))
recon_error_knn$Cat = factor(recon_error_knn$Cat, levels=c(0,1))
sfig_B = ggplot(recon_error_knn, aes(x = factor(Parameter), y = Error)) +
  geom_boxplot(aes(fill=Cat),outlier.size=0.5) + xlab('Value of k') +
  scale_y_continuous(limits = c(0, 2)) + gtex_v8_figure_theme() +
  scale_fill_manual(values=c('white', '#E5CEDC')) + guides(fill=F) +
  theme(axis.text.x=element_text(hjust=1,angle=90))

### panel C - comparison of enrichments when imputing vs not ###
knn_enrich = fread(paste0(data_dir, 'figS26C_input_data.txt'))
knn_enrich$PT = factor(knn_enrich$PT, levels=pts)
sfig_C = ggplot(knn_enrich, aes(x=PT,y=Riskratio,Group=Method)) + 
  geom_point(size=2,aes(color=Method),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + 
  scale_color_manual(values=c('#7FCDBB', 'black')) +
  ggtitle('') + theme_bw() + xlab('P-value threshold') +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() +
  facet_wrap(.~Type,scales='free') +
  theme(axis.text.x=element_text(hjust=1,angle=40)) +
  theme(strip.background = element_blank(),
        strip.text.x=element_text(size=8),
        legend.title=element_blank())

### combine
first_row <- plot_grid(sfig_A, sfig_B, sfig_C, labels = c('A','B', 'C'), ncol=3, rel_widths=c(1,1,2), axis='tblr', align='h')

ggsave(first_row, file=paste0(out_dir, 'figS26.pdf'), width=7.2, height=3,units="in")





