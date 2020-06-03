# Processed input data required to generate figures

**Below you can find a short description of each input file that we have provided**



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


