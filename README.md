# GTEx version 8 rare variant figure generation
This repository will re-generate all main and supplementary figures in "Transcriptomic signatures across human tissues identify functional rare genetic variation" from a collection of processed input data, which can be found in the processed_input_data directory.

Since the publication of this paper, we found some small mistakes in the code related to the correction of total expression data. Specifically, expression data was corrected for the first 20 genotype PCs rather than the top 3, and genes after chr9 were not corrected for their strongest cis-eQTL. This did not impact any conclusions from our analyses, but for future use, we recommend correcting for the top three genotype PCs and all cis-eQTLs in future work. Updated multi-tissue eOutlier calls are available on the GTEx portal. We thank Taibo Li for helping to identify these inconsistencies.
## Dependencies and Testing Envivonment 1

* R 3.5.1
* R package ggplot2 (3.3.0)
* R package cowplot (1.0.0)
* R package reshape (0.8.8)
* R package RColorBrewer (1.1-2)
* R package ggthemes (4.2.0)
* R package plyr (1.8.6)
* R package ggseqlogo (0.1)

The following scripts were generated using **Dependencies and Testing Environment 1**:
* plot_figure2.R
* plot_figure4.R
* plot_figureS4.R
* plot_figureS5.R
* plot_figureS6.R
* plot_figureS17.R
* plot_figureS18.R
* plot_figureS19.R
* plot_figureS20.R
* plot_figureS24.R
* plot_figureS27.R
* plot_figureS28.R
* plot_figureS29.R
* plot_figureS30.R
* plot_figureS31.R
* plot_figureS32.R
* plot_figureS33.R
* plot_figureS34.R
* plot_figureS35.R
* plot_figureS36.R
* plot_figureS37.R


## Dependencies and Testing Environment 2


## Authors

* **Nicole Ferraro** -- [nmferraro5](https://github.com/nmferraro5) -- nferraro@stanford.edu 

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes) -- bstrober3@gmail.com






