# This figure plots the results from the polyclonal ASE based CRISPR assay, for
# experimental validation of high watershed stop-gained variants

# Jonah Einson
# jeinson@nygenome.org

library(tidyverse)
gtex_v8_figure_theme <- function() {
  return(
    theme(
      plot.title = element_text(face="plain",size=8), 
      text = element_text(size=8),
      axis.text=element_text(size=7), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"), 
      legend.text = element_text(size=7), 
      legend.title = element_text(size=8)
    )
  )
}

class_colors <- c(
  rgb(121, 92, 129, maxColorValue = 255), 
  rgb(23, 50, 75, maxColorValue = 255), 
  rgb(193, 205, 222, maxColorValue = 255)
)
names(class_colors) <- c("ASE", "Splicing", "Total Expression")


dat <- read_tsv("processed_input_data/figureS39/figS39_input_data.txt")
dat <- filter(dat, !is.na(mean_aFC))
dat$pbonf <- p.adjust(dat$pvalues, method = "fdr")

# Tidy up dat
dat_tidy <- dat %>%
  pivot_longer(
    cols = ends_with("posterior"), 
    names_to = "method",
    values_to = "posterior_probability"
  )
# Make a column to match to the color key
key <- c("Splicing", "Total Expression", "ASE")
names(key) <- unique(dat_tidy$method)
dat_tidy$method_clean <- key[dat_tidy$method]


# Plot one more with median watershed/ase score

dat <- 
  mutate(dat, med_ase_expr = rowMeans(cbind(ase_watershed_posterior, total_expression_watershed_posterior)))

dat$vtype_sig <- dat$variant_type
dat$vtype_sig[dat$pbonf > .05] <- "stop_gained_nsig"

png("generated_figures/figS39.png", width = 2.5, height = 3, units = "in", res = 300)
ggplot(dat, 
       aes(mean_aFC, med_ase_expr, 
           color = vtype_sig)) + 
  geom_point(pch = 16) +
  labs(color = "Variant annotation") +
  xlab("Mean aFC") + ylab("Mean ASE and Total Expression \nWatershed posterior probability") +
  scale_color_manual(values = c("blue", "red", "pink"), 
                     labels = c(
                       "Non-eQTL control (6)", 
                       "Stop Gained, p < .05 (13)", 
                       "Stop Gained, p > .05 (1)")) +
  geom_hline(yintercept = .9, col = "darkred", lty = 5) +
  geom_vline(xintercept = 0, col = "grey", lty = 2) + 
  coord_flip() +
  scale_x_reverse() + 
  theme(legend.position = c(.4, .82)) +
  gtex_v8_figure_theme()
dev.off()

