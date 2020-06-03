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
