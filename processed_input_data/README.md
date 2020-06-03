# Processed input data required to generate figures

**Below you can find a short description of each input file that we have provided**


## Figure 1
* **fig1b_input_data.txt**
	* Data frame containing counts across outliers and controls of rare variants nearby outlier genes at varying allele frequencies
	* Relevent column desciptors (base 1):
		* Column1: "Control_0": Count of control individuals with no rare variant nearby outlier genes
		* Column2: "Control_1": Count of control individuals with a rare variant nearby outlier genes
		* Column3: "Outlier_0": Count of outlier individuals with no rare variant nearby outlier genes
		* Column4: "Outlier_1": Count of outlier individuals with rare variant nearby outlier genes
		* Column5: "MAF": Allele frequency bin used to define the set of variants included (either novel, single, double, rare, or low frequency)
		* Column6: "Method": Type of outlier (either eOutlier, sOutlier, or aseOutlier)
		* Column7: "Type": Type of variant (either SNVs+indels or SVs)
* **fig1c_input_data.txt**
	* Data frame containing counts across outliers and controls of rare variants nearby outlier genes across variant categories
	* Relevent column desciptors (base 1):
		* Column1: "Control_0": Count of control individuals with no rare variant nearby outlier genes
		* Column2: "Control_1": Count of control individuals with a rare variant nearby outlier genes
		* Column3: "Outlier_0": Count of outlier individuals with no rare variant nearby outlier genes
		* Column4: "Outlier_1": Count of outlier individuals with rare variant nearby outlier genes
		* Column5: "Variant_category": Variant annotation category
		* Column6: "Method": Type of outlier (either eOutlier, sOutlier, or aseOutlier)
* **fig1d_input_data.txt**
	* Data frame containing proportion of outliers in each bin with a nearby rare variant of the given category
	* Relevent column descriptors (base 1):
		* Column1: "Variant_category": Variant annotation category
		* Column2: "pval_bin": Outlier bin based on outlier p-value
		* Column3: "Proportion": Proportion of instances in the given bin with nearby rare variants of the given category
		* Column4: "Type": Type of outlier (either eOutlier, sOutlier, or aseOutlier)
* **fig1e_input_data.txt**
	* Data frame containing proportion of rare variants of different annotations leading to outlier signal in nearby genes
	* Relevent column descriptors (base 1)
		* Column1: "Feature": Variant annotation category 
		* Column2: "Type": Type of outlier (either eOutlier, sOutlier, or aseOutlier)
		* Column3: "Proportion": Proportion of instances of rare variants in a given category leading to outlier signal in nearby genes
* **variant_type_colors.rds**
	* Character object including hex color codes, with names indicating the variant category that color is used for in the above plots
		
## Figure 2
* **fig2a_input_data.txt**
	* Data frame containing riskratio of SNVs, indels, and SVs at varying distances upstream of outlier genes across data types
	* Relevent column desciptors (base 1):
		* Column1: "Riskratio": Relative risk of rare variant lying in outlier gene
		* Column6: "Method": Outlier type (either sOutlier, eOutlier, or aseOutlier)
		* Column7: "Type": The variant class/type (either SVs, indels, or SNVs)
		* Column9: "Window": The distance upstream from a gene (KB)
* **fig2b_input_data.txt**
	* Data frame containing proportion of eOutliers with TSS rare variants in various promotor motifs
	* Relevent column descriptors (base 1):
		* Column3: "promotor_motif": The name of the promotor motif
		* Column10: "ExpBin": The eOutlier bin. (either Under, Control, or Over)
		* Column12: "NumBin": The proportion of rare variants in the given promotor for the given eOutlier bin
* **fig2d_inliers_input_data.txt** and **fig2d_outliers_input_data.txt**
	* Both of these files are formatted in the same way
	* fig2d_inliers_input_data.txt contains positional location of rare variants near inliers splice sites
	* fig2d_inliers_input_data.txt contains positional location of rare variants near outlier splice sites
	* Relevent column descriptors (base 1):
		* Column1: "splice_site_type": Whether it is a donor or acceptor splice site
		* Column2: "distance": Positional location of rare variant around the splice site
		* Column3: "annotated_splice_site": Whether the splice site is annotated or novel
* **fig2e_input_data.txt**
	* Data frame containing types of mutations (the nucleotide changes) occuring at each position around splice sites
	* Relevent column descriptors (base 1)
		* Column1: "mutations": The nucleotide change corresponding to the rare variant 
		* Column2: "distance_to_ss": The position (relative to the nearest splice site) of the rare variant
		* Column3: "odds_ratio": The exp(junction usage of the splice site).. ie. not on log scale
* **fig2f_input_data.txt**
	* Data frame containing junction usage of the closest splice sites to variants that lie in PPT
	* Relevent column descriptors (base 1)
		* Column1: "odds_ratio": The junction usage of the splice site
		* Column2: "type": the type of rare variant. Either "Purine to Pyrimidine" or "Pyrimidine to Purine"

## Figure 3
* **fig3a_input_data.txt**
	* Data frame containing the median proportion of single tissue outliers that replicate across all other tissues
	* Relevent column desciptors (base 1):
		* Column1: "tiss": Short tissue code (i.e. ADPSBQ)
		* Column2: "frac": Median proportion of outliers from given tissue that replicate in other tissues
		* Column3: "min": Minimum proportion of outliers from given tissue that replicate in other tissues
		* Column4: "max": Maximum proportion of outliers from given tissue that replicate in other tissues
		* Column5: "type": Type of outlier (either eOutlier, sOutlier, or aseOutlier)
		* Column6: "tissue_name": Long tissue name (i.e. Adipose - Subcutaneous)
* **fig3b_input_data.txt**
	* Data frame containing estimated relative risk of rare variants nearby single tissue outliers across tissues at varying outlier p-value thresholds
	* Relevent column desciptors (base 1):
		* Column1: "tissue": Tissue name
		* Column2: "ratio": Relative risk estimate of rare variants nearby outliers in the given tissue
		* Column3: "lower.q": Lower bound of 95% confidence interval around estimate
		* Column4: "upper.q": Upper bound of 95% confidence interval around estimate
		* Column10: "sig": Outlier p-value threshold
		* Column11: "outlier_type": Type of outlier (either eOutlier, sOutlier, or aseOutlier)
* **fig3c_input_data.txt**
	* Data frame containing estimated relative risk of rare variants of varying annotation categories nearby single tissue outliers across tissues
	* Relevent column desciptors (base 1):
		* Column1: "tissue": Tissue name
		* Column2: "ratio": Relative risk estimate of rare variants nearby outliers in the given tissue
		* Column3: "lower.q": Lower bound of 95% confidence interval around estimate
		* Column4: "upper.q": Upper bound of 95% confidence interval around estimate
		* Column9: "Variant_type": Variant annotation
		* Column11: "outlier_type": Type of outlier (either eOutlier, sOutlier, or aseOutlier)
* **fig3d_input_data.txt**
	* Data frame containing proportion of tissues supporting outlier calls underlying multi-tissue eOutliers and correlation outliers
	* Relevent column descriptors (base 1)
		* Column3: "Method": Outlier type (MEDZ or correlation) 
		* Column4: "npass": Number of tissues with |Z| > 3 for that gene-individual outlier
		* Column3: "prop": Proportion (range [0,1]) of measured tissue with outlier signal for that outlier
* **fig3e_input_data.txt**
	* Data frame containing relative risk estimates for rare variants occurring in tissue-specific enhancers nearby tissue-specific outliers
	* Relevent column descriptors (base 1)
		* Column1: "Riskratio": Relative risk estimate
		* Column2: "Lower": Lower bound of 95% confidence interval around estimate
		* Column3: "Upper": Upper bound of 95% confidence interval around estimate
		* Column5: "Type": Indicates whether the relative risk was calculated using enhancer regions matched to the tissue driving the outlier signal or not (either Matched or Unmatched)
		* Column6: "Zscore": |Z-score| threshold used to determine tissue-specificity of correlation outliers
		
## Figure 4
* **fig4b_watershed_edge_weights_input_data.txt**
	* File containing Watershed learned edge weights 
	* Relevent column descriptors (base 1)
		* Column1: Edge weight between Splicing and Expression outlier signals
		* Column2: Edge weight between Splicing and ASE outlier signals
		* Column3: Edge weight between Expresssion and ASE outlier signals
* **fig4c_absolute_risk_input_data.txt**
	* Data frame containing proportion of prioritized rare variants leading to outlier
	* Relevent column descriptors (base 1)
		* Column2: "absolute_risk": Proportion of rare variants leading to an outlier
		* Column3: "outlier_type": Type of outlier (either eOutlier, sOutlier, or aseOutlier)
		* Column4: "model_type": Whether using Watershed or GAM
* **fig4d_pr_curve_input_data.txt**
	* Data frame containing precision and recall of Watershed, RIVER, and GAM for each outlier tyep
	* Relevent column descriptors (base 1):
		* Column1: "precision"
		* Column2: "recall"
		* Column3: "outlier_type": either splicing, ase, or expression
		* Column4: "prediction_type": either Watershed, RIVER, or GAM
* **fig4e_tissue_watershed_edge_weights_input_data.txt**
	* File containing learned tissue Watershed edge weights
	* This file is a matrix of dimension NXN where N is the number of tissues
	* Each element of this matrix is the edge weight corresponding to the connection between the two tissues
* **fig4e_gtex_tissue_colors_input_data.txt**
	* File containing mapping from gtex tissue name to a specific color
* **fig4f_tissue_watershed_pr_auc_input_data.txt**
	* File containing area under the precision recall curve for several tissue methods used
	* Each row corresponds to a tissue type
	* Relevent column descriptors (base 1):
		* Column1: "outlier_type": either ASE, Splicing, or Expression
		* Column2: "auc": Area under the precision-recall curve
		* Column3: "method": The model used: either Watershed, RIVER, or GAM

## Figure 5
* **fig5a_input_data.txt**
	* File containing number of outliers per individual across different filters
	* Relevent column descriptors (base 1)
		* Column2: "Type": Type of outlier (either eOutlier, sOutlier, or aseOutlier)
		* Column3: "variable": Filter applied (either Outliers, Outlier RV, Watershed > 0.5, Watershed > 0.9)
		* Column4: "value": Count per individual
* **fig5b_input_data.txt**
	* Data frame containing variant effect size percentiles in UKBB traits across categories
	* Relevent column descriptors (base 1)
		* Column1: "VID": Variant chromosome and position (chr:pos)
		* Column2: "VPercentile": Percentile of variant effect size for the given trait
		* Column3: "Gene": Associated gene 
		* Column4: "Trait": UKBB trait code
		* Column5: "Category": Category of the given variant-gene (either Non-outlier, Outlier, or Coloc Outlier)
* **fig5c_input_data.txt**
	* Data frame containing the proportion of variants falling in the top 25th percentile in co-localized regions for relevant traits split by Watershed posterior
	* Relevent column descriptors (base 1):
		* Column1: "Posterior": Watershed posterior threshold
		* Column2: "Category": Indicates whether the values are from a random permutation or actual
		* Column3: "Prop": Proportion of variants in top 25th percentile at the given threshold
		* Column4: "Model": either Watershed or CADD
* **fig5d_pval_input_data.txt**
	* File containing GWAS p-values for all variants in a given co-localized region for Asthma
	* Relevent column descriptors (base 1):
		* Column1: "Chr": Chromosome
		* Column2: "Pos": Position
		* Column3: "pval": GWAS association p-value
		* Column5: "IsVar": Indicates whether the variant is the outlier-associated variant with a high Watershed score
* **fig5d_betas_input_data.txt**
	* File containing GWAS effect sizes for all variants in a given co-localized region for Asthma
	* Relevent column descriptors (base 1):
		* Column2: "minor_AF": Variant allele frequency in UKBB
		* Column9: "VOI": Indicates whether variant is high-scoring Watershed variant
		* Column10: "beta_scaled": GWAS effect size scaled by case-control ratio
* **fig5e_pval_input_data.txt**
	* File containing GWAS p-values for all variants in a given co-localized region for High cholesterol
	* Relevent column descriptors (base 1):
		* Column1: "Chr": Chromosome
		* Column2: "Pos": Position
		* Column3: "pval": GWAS association p-value
		* Column5: "IsVar": Indicates whether the variant is the outlier-associated variant with a high Watershed score
* **fig5e_betas_input_data.txt**
	* File containing GWAS effect sizes for all variants in a given co-localized region for Asthma
	* Relevent column descriptors (base 1):
		* Column2: "minor_AF": Variant allele frequency in UKBB
		* Column9: "VOI": Indicates whether variant is high-scoring Watershed variant
		* Column10: "beta_scaled": GWAS effect size scaled by case-control ratio

## Figure S4
* **figS4_input_data.txt**
	* Data frame containing information on the number of observed exon-exon junctions, LeafCutter clusters, and genes in each tissue
	* Relevent column descriptors (base 1):
		* Column1: "tissue"
		* Column2: "number_of_junctions": Number of observed exon-exon junctions in this tissue
		* Column3: "number_of_clusters": Number of observed LeafCutter clusters in this tissue
		* Column4: "number_of_genes": Number of observed genes that the LeafCutter clusters map to in this tissue

## Figure S5
* **figS5_input_data.txt**
	* Data frame containing SPOT gene-gene level corrected and uncorrected p-values in Muscle-Skeletal tissue
	* Each row of data frame corresponds to a gene-individual pair
	* Relevent column descriptors (base 1):
		* Column 1: "corrected_pvalue": SPOT gene-level corrected p-value in Muscle-skeletal
		* Column 2: "uncorrected_pvalue": SPOT gene-level uncorrected p-value in Muscle-skeletal
		* Column 3: "number_of_clusters": The number of LeafCutter Clusters assigned to this gene

## Figure S6
* **figS6_input_data.txt**
	* Data frame containing SPOT LeafCutter cluster-level p-values in Muscle-Skeletal tissue for 4 different versions of SPOT
	* Each row of the data frame corresponds to a LeafCutter cluster-individual pair
	* Relevent column descriptors (base 1):
		* Column1: "standard_prior_20K_reads": SPOT LeafCutter cluster p-value using default version of SPOT
		* Column2: "standard_prior_10K_reads": SPOT LeafCutter cluster p-value using version of SPOT with 10,000 simulated reads
		* Column3: "standard_prior_100K_reads": SPOT LeafCutter cluster p-value using version of SPOT with 100,000 simulated reads
		* Column4: "no_prior_20K_reads": SPOT LeafCutter cluster p-value using version of SPOT with no prior placed on the alphas
		* Column5: "fraction": Maximum fraction of reads mapping to a single exon-exon junction in the LeafCutter cluster


## Figure S17
* **figS17a_input_data.txt**
	* Data frame containing data necessary to compute relative risk of rare variants within various window sizes around splice sites for sOutliers
	* Relevent column descriptors (base 1):
		* Column1: "pvalues": p-value threshold used by SPOT to generate sOutliers
		* Column2: "distance": Distance window around splice site
		* Column3: "num_outlier_clusters": The number of LeafCutter clusters that are categorized as sOutliers
		* Column4: "num_outlier_clusters_with_rv": The number of LeafCutter clusterrs that are categorized as sOutliers that have a nearby rare variant
		* Column5: "num_inlier_clusters": The number of LeafCutter clusters that are categorized as not sOutliers
		* Column6: "num_inlier_clusters_with_rv": The number of LeafCutter clusters that are categorized as not sOutliers that have a nearby rare variant
* **figS17b_input_data.txt**
	* Data frame containing junction usage of the closest splice sites to variants that lie in consensus sequence
	* Relevent column descriptors (base 1):
		* Column1: "odds_ratio": The junction usage of the splice site
		* Column2: "type": The type of rare variant in the consensus sequence. Does it create the consensus or destroy the concensus sequence.

## Figure S18
* **figS18_input_data.txt**
	* Data frame containing types of mutations (the nucleotide changes) occuring at each position around splice sites
	* Relevent column descriptors (base 1)
		* Column1: "mutations": The nucleotide change corresponding to the rare variant 
		* Column2: "distance_to_ss": The position (relative to the nearest splice site) of the rare variant
		* Column3: "odds_ratio": The exp(junction usage of the splice site).. ie. not on log scale

## Figure S19
* **figS19_outliers_input_data.txt** and **figS19_inliers_input_data.txt**
	* Both of these files are formatted in the same way
	* figS19_inliers_input_data.txt contains positional location of rare variants near inliers splice sites
	* figS19_outliers_input_data.txt contains positional location of rare variants near outlier splice sites
	* Relevent column descriptors (base 1)
		* Column1: "splice_site_type": Whether nearby splice site to variant is a donor or acceptor splice site
		* Column2: "distance": Positional location of variant around the splice site
		* Column3: "annotated_splice_site": Boolean variable concerning whether splice site is annotated or novel
		* Column4: "variant_allele": Nucleotide of variant allele
		* Column5: "major_allele": Nucleotide of common allele

## Figure S20
* **figS20_outliers_input_data.txt** and **figS20_inliers_input_data.txt**
	* Both of these files are formatted in the same way
	* figS20_inliers_input_data.txt contains positional location of rare variants near inliers splice sites
	* figS20_outliers_input_data.txt contains positional location of rare variants near outlier splice sites
	* Relevent column descriptors (base 1)
		* Column1: "splice_site_type": Whether nearby splice site to variant is a donor or acceptor splice site
		* Column2: "distance": Positional location of variant around the splice site
		* Column3: "annotated_splice_site": Boolean variable concerning whether splice site is annotated or novel
		* Column4: "variant_allele": Nucleotide of variant allele
		* Column5: "major_allele": Nucleotide of common allele


## Figure S24
* **figS24_tissue_names_input_data.txt**
	* File containing ordered list of gtex tissue names
* **figS24_tissue_colors_input_data.txt**
	* File containing mapping from gtex tissue names to colors
* **figS24_tissue_by_tissue_variant_sOutlier_overlap_input_data.txt**
	* Data frame containing contingency table necessary to compute relative risk of rare variants within sOutliers in each tissue
	* Relevent column descriptors (base 1)
		* Column1: "Tissue_name"
		* Column2: "number_of_variants_in_outliers": Number of sOutlier LeafCutter clusters in this tissue with a nearby rare variant
		* Column3: "number_of_outliers": Number of sOutlier LeafCutter clusters in this tissue
		* Column4: "number_of_variants_in_inliers": Number of sInlier LeafCutter clusters in this tissue with a nearby rare variant
		* Column5: "number_of_inliers": Number of sInlier LeafCutter clusters in this tissue
* **figS24_tissue_by_tissue_concensus_jxn_usage_input_data.txt**
	* Data frame containing junction usage of the closest splice sites to variants that lie in consensus sequence
	* Relevent column descriptors (base 1):
		* Column1: "odds_ratio": The junction usage of the splice site
		* Column2: "type": The type of rare variant in the consensus sequence. Does it create the consensus or destroy the concensus sequence.
* **figS24_tissue_by_tissue_ppt_jxn_usage_input_data.txt**
	* Data frame containing junction usage of the closest splice sites to variants that lie in PPT
	* Relevent column descriptors (base 1):
		* Column1: "odds_ratio": The junction usage of the splice site
		* Column2: "type": The type of rare variant in the PPT. Does it create the consensus or destroy the PPT


## Figure S27
* **figS27_input_data.txt**
	* Data frame containing precision and recall of Watershed and CADD for each outlier type
	* Relevent column descriptors (base 1):
		* Column1: "precision"
		* Column2: "recall"
		* Column3: "outlier_type": either splicing, ase, or expression
		* Column4: "prediction_type": either Watershed or CADD


## Figure S28
* **figS28a_input_data.txt**
	* All 6 files (a,b,c,d,e,f) pertaining figure S28 have the same file format
	* Each file is a data frame containing precision and recall of Watershed, RIVER, and GAM for each outlier type
	* Relevent column descriptors (base 1):
		* Column1: "precision"
		* Column2: "recall"
		* Column3: "outlier_type": either splicing, ase, or expression
		* Column4: "prediction_type": either Watershed, RIVER, or GAM

## Figure S29
* **figS29_river_input_data.txt**
	* Confusion matrix for RIVER in jointly predicting outlier status of all three outlier signals (class) using held out pairs of individuals. 
* **figS29_watershed_exact_input_data.txt**
	* Confusion matrix for Watershed (with exact inference) in jointly predicting outlier status of all three outlier signals (class) using held out pairs of individuals. 
* **figS29_watershed_approximate_input_data.txt**
	* Confusion matrix for Watershed (with approximate inference) in jointly predicting outlier status of all three outlier signals (class) using held out pairs of individuals. 
* For all three files, the first element of the binary class abbreviation labels represents median splicing outlier status, the second element of the class abbreviations represents median expression outlier status, and the third element of the class abbreviations represents ASE outlier status.

## Figure S30
* **figS30_input_data.txt**
	* Data frame containing proportion of prioritized rare variants leading to outlier
	* Relevent column descriptors (base 1)
		* Column1: "watershed_threshold": Watershed threshold used to prioritize variants
		* Column2: "absolute_risk": Proportion of rare variants leading to an outlier
		* Column3: "outlier_type": Type of outlier (either eOutlier, sOutlier, or aseOutlier)
		* Column4: "model_type": Whether using Watershed, CADD, or GAM

## Figure S31
* **figS31a_input_data.txt**
	* Data frame containing Watershed genomic annotation coefficients (when model was trained via exact inference and when model was trained via approximate inference)
	* Relevent column descriptors (base 1)
		* Column1: "exact_betas": Genomic annotation coefficient when model was trained via exact inference
		* Column2: "approximate_betas": Genomic annotation coefficient when model was trained via approximate inference
		* Column3: "outlier_class": Outlier type (ase, splicing, or expression) for which genomic annotation coefficent corresponds to
* **figS31b_input_data.txt**
	* Data frame containing precision and recall of Watershed-exact, Watershed-approximate, and RIVER for each outlier type
	* Relevent column descriptors (base 1):
		* Column1: "precision"
		* Column2: "recall"
		* Column3: "outlier_type": either splicing, ase, or expression
		* Column4: "inference": Either exact or approximate
		* Column5: "prediction_type": either Watershed or RIVER

## Figure S32
* **figS32_ase_input_data.txt** and **figS32_splicing_input_data.txt** and **figS32_expression_input_data.txt**
	* One file for each of the three outlier types
	* Each file is formatted identically
		* Each files contains learned tissue Watershed edge weights
		* Each file is a matrix of dimension NXN where N is the number of tissues
		* Each element of this matrix is the edge weight corresponding to the connection between the two tissues

## Figure S33
* **figS33_ase_input_data.txt** and **figS33_splicing_input_data.txt** and **figS33_expression_input_data.txt**
	* One file for each of the three outlier types
	* Each file is formatted identially
	* Relevent column descriptors (base 1):
		* Column1: "auc_watershed": Aread under precision recall for tissue Watershed in a specific tissue
		* Column2: "auc_river": Aread under precision recall for tissue RIVER in a specific tissue
		* Column3: "tissue"

## Figure S34
* **figS34_ase_input_data.txt** and **figS34_splicing_input_data.txt** and **figS34_expression_input_data.txt**
	* One file for each of the three outlier types
	* Each file is formatted identially
	* Relevent column descriptors (base 1):
		* Column1: "tissue"
		* Column3: "delta_auc": Difference in the area under the precision recall curves between tissue-Watershed and tissue-RIVER for a single tissue
		* Column4: "lower_bound": Lower bound of 95% CI on delta_auc generated by non-parametric bootstrapping
		* Column5: "upper_bound": Upper bound of 95% CI on delta_auc generated by non-parametric bootstrapping


## Figure S35
* **figS35_ase_input_data.txt** and **figS35_splicing_input_data.txt** and **figS35_expression_input_data.txt**
	* One file for each of the three outlier types
	* Each file is formatted identially
	* Relevent column descriptors (base 1):
		* Column1: "auc_watershed": Aread under precision recall for tissue Watershed in a specific tissue
		* Column2: "auc_river": Aread under precision recall for RIVER (trained on median outlier signal) evaluated in a specific tissue
		* Column3: "tissue"

## Figure S36
* **figS36_ase_input_data.txt** and **figS36_splicing_input_data.txt** and **figS36_expression_input_data.txt**
	* One file for each of the three outlier types
	* Each file is formatted identially
	* Relevent column descriptors (base 1):
		* Column1: "tissue"
		* Column3: "delta_auc": Difference in the area under the precision recall curves between tissue-Watershed and RIVER (trained on median outlier signal) for a single tissue
		* Column4: "lower_bound": Lower bound of 95% CI on delta_auc generated by non-parametric bootstrapping
		* Column5: "upper_bound": Upper bound of 95% CI on delta_auc generated by non-parametric bootstrapping


## Figure S37
* **figS37_ase_input_data.txt** and **figS37_splicing_input_data.txt** and **figS37_expression_input_data.txt**
	* One file for each of the three outlier types
	* Each file is formatted identially
	* Relevent column descriptors (base 1):
		* Column1: "pval": -log10(outlier_pvalue) of gene nearby rare variant in ASMAD cohort
		* Column2: "type": GTEx Watershed posterior bin of rare variant
		* Column3: "model": Type of model used in gtex. Either Watershed or GAM

