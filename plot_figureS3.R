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

data_dir = '~/gtex_v8_rare_variant_figure_generation/processed_input_data/figureS3/'
out_dir = '~/gtex_v8_rare_variant_figure_generation/generated_figures/'

### panel A - comparison of Vg estimates in GTEx v7 and v8 in Adipose ###
vg_data = fread(paste0(data_dir, 'figS3A_input_data.txt'))
sfig_A <- ggplot(vg_data, aes(x=vgs_v7, y=vgs_v8)) + geom_point() +
  geom_abline(color='red') + gtex_v8_figure_theme() +
  xlab(bquote('sqrt('*V^G*') in v7')) +
  ylab(bquote('sqrt('*V^G*') in v8')) 

### panel B - distribution of spearman correlations of Vg estimates from v7 to v8 ###
sp_data = fread(paste0(data_dir, 'figS3B_input_data.txt'))
sfig_B <- ggplot(sp_data, aes(x=Cat,y=spearman_correlation)) +  
  geom_violin(fill='lightgrey') + geom_boxplot(width=0.5) +
  gtex_v8_figure_theme() + xlab('') + ylab('Per tissue spearman correlation of v7 and v8') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())

### panel C - number of genes with available Vg estimates ###
num_genes = fread(paste0(data_dir, 'figS3C_input_data.txt'))
num_genes$FullName = factor(num_genes$FullName, levels=rev(unique(num_genes$FullName)))
sfig_C <- ggplot(num_genes, aes(x=FullName, y=value, Group=Version)) +
  geom_bar(stat='identity',aes(fill=Version),position='dodge',alpha=0.8) + gtex_v8_figure_theme() +
  ylab('Number of genes') + xlab('') +
  scale_fill_manual(values=vcols) +
  theme(axis.text.x=element_text(hjust=1,angle=90),
        legend.position=c(0.05,0.9),
        legend.title=element_blank())


first_row <- plot_grid(sfig_A, sfig_B, labels=c('A','B'), ncol=2, align='hv', axis='tblr', rel_widths = c(2,1))
sfig <- plot_grid(first_row, sfig_C, labels=c('', 'C'), nrow=2)

ggsave(sfig, file=paste0(out_dir, 'figS3.pdf'), width=7.2, height=8, units="in")



