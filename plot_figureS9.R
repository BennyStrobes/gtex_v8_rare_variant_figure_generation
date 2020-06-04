library(ggplot2)
library(ggridges)
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

data_dir = 'processed_input_data/figureS9/'
out_dir = 'generated_figures/'

### panel A - sharing aross outlier types ###
sharing_df = fread(paste0(data_dir, 'figS9a_input_data.txt'))
sharing_df$Test = factor(sharing_df$Test, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sharing_df$Discovery = factor(sharing_df$Discovery, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_A = ggplot(sharing_df, aes(y = Test, x = Discovery, fill = value)) +
  geom_tile(color = "grey") +
  xlab('Discovery method') + ylab('Test method') +
  geom_text(aes(label = round(value, 2)),cex=3) +
  scale_fill_gradient("Sharing \nfraction",
                      low = "white", high = "#BC6F8B",
                      na.value = "black", limits = c(0, .52)) +
  gtex_v8_figure_theme()

### panel B - proportion of outliers shared across all three types with nearby rare variant ###
shared_annot = fread(paste0(data_dir, 'figS9b_input_data.txt'))
shared_annot$Var1 = factor(shared_annot$Var1, levels=c('no_variant', 'DEL', 'conserved_noncoding', 'other_noncoding', 'splice'))
sfig_2B = ggplot(shared_annot, aes(x=Var1,y=Freq/35)) + geom_bar(stat='identity') +
  xlab('') + ylab('Proportion with nearby rare variant') +
  gtex_v8_figure_theme() + 
  theme(axis.text.x=element_text(hjust=1,angle=45))

### panel C - Z-score distributions across outlier groups ###
sharing_data = fread(paste0(data_dir, 'figureS9c_input_data.txt'))
sharing_data$Cat = factor(sharing_data$Cat, levels=c('Shared Outlier (n=319)', 'eOutlier (n=305)', 'aseOutlier (n=1890)', 'Non-outlier (n=496536)'))
sfig_C = ggplot(sharing_data, aes(y = Cat, x = abs(MedZ))) + 
  stat_density_ridges(quantiles=2, quantile_lines=T) +
  gtex_v8_figure_theme() + ylab('') + xlab('|Median Z-score|') 

### panel D - proportion of aseOutliers with nearby rare variants split by expression effect ###
ase_summary = fread(paste0(data_dir, 'figS9d_input_data.txt'))
ase_summary$medz_bin = factor(ase_summary$medz_bin,
                               levels=c('-3~-2','-2~-1','-1~1','1~2','2~3'))
ase_summary$variant_cat = factor(ase_summary$variant_cat,
                                  levels=c('no_variant','other_noncoding','coding','conserved_noncoding','TSS','stop','frameshift','splice','TE','INV','BND','DEL','CNV','DUP'))

sfig_D = ggplot(ase_summary, aes(x=medz_bin,y=NProp,Group=variant_cat)) +
  geom_bar(aes(fill=variant_cat),color='black',stat='identity') +
  theme_bw() + coord_flip() + ylab('') + xlab('') +
  scale_fill_manual(values=plot_cols) +
  ggtitle('aseOutliers') + xlab('Median Z-score bin') +
  gtex_v8_figure_theme() +
  theme(legend.title=element_blank(),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width = unit(0.1, "cm"))

### combine
first_row <- plot_grid(sfig_A, sfig_B, labels = c('A','B'), ncol=2, align='h')
second_row <- plot_grid(sfig_C, sfig_D, labels = c('C','D'), ncol=2, align='h')
sfig = plot_grid(first_row, second_row, nrow = 2, align='v')

ggsave(sfig, file=paste0(out_dir, 'figS9.pdf'), width=7.2, height=7.2,units="in")
