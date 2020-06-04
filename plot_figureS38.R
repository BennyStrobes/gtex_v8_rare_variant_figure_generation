library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = 'processed_input_data/figureS38/'
out_dir = 'generated_figures/'

### compare effect sizes between high and low Watershed score rare variants ###
mpra_data = fread(paste0(data_dir, 'figS38_input_data.txt'))
mpval = wilcox.test(abs(filter(mpra_data, WS_Bin == 'High')$log2FoldChange_allele),
            abs(filter(mpra_data, WS_Bin == 'Low')$log2FoldChange_allele),
            alternative='g')$p.value

sf1 = ggplot(mpra_data, aes(x=WS_Bin, y=abs(log2FoldChange_allele))) + 
  geom_boxplot() + theme_bw() + 
  gtex_v8_figure_theme() + scale_y_log10() +
  xlab('Watershed posterior bin') + ylab('|log(allelic fold-change)|') +
  annotate("text", x=1.5, y=1, label=paste0('p= ', round(mpval, 3)), cex=3)

ggsave(sf1, file=paste0(out_dir, 'figS38.pdf'), width=7.2, height=5, units="in")

