# LINDA Example

This repository shows the steps needed to be performed by the user in order to run the LINDA analysis from the gene expression data and exon skipping events. The example here, refers to studying the splicing effects over the protein interaction networks up U2AF2 Knockdown of HepG2 cells (data from [ENCORE project](https://www.encodeproject.org/encore-matrix/?type=Experiment&status=released&internal_tags=ENCORE)). For a better understanding of the LINDA pipeline, please refer to the documentation provided on its [dedicated R-Package GitHub repository](https://github.com/dieterich-lab/LINDA).

The repository is organized with the following directories:

**1. Data:** Here are provided the gene exppression data (*gene_expression_data.csv*) showing the expression of individual genes across two replicates of HepG2-Control (Ctrl) as well as HepG2-U2AF2-KnockDown (KD). Additionally we provide a list of Exon Skipping (ES) events (*SE.MATS.JunctionCountOnly.txt*) obtained from [rMATS](https://www.pnas.org/content/111/51/E5593) for the KDvsCtrl comparison.

**2. Differential Splicing:** Here are provided scripts (*/Differential_Splicing/analysis_script.R*) which have been used to assign an Exon ID to all the genomic coordinates in the identified exon skipping events (*/Data/SE.MATS.JunctionCountOnly.txt*). For more details, please see documentation provided in */Differential_Splicing/README.md*.

**3. Gene Expression:** Here are provided the scripts (*/Gene_Expression/analysis_script.R*) which have been used to estimate the differentially expressed genes (DGE's) for the KDvsCtrl comparison as well as TF activity estimation with [viper](https://www.bioconductor.org/packages/release/bioc/html/viper.html) and [DoRothEA](https://www.bioconductor.org/packages/release/bioc/html/viper.html). For more details, please see documentation provided in */Gene_Expression/README.md*.

**4. LINDA Analysis:** Here are provided the scripts (*/LINDA_Analysis/analysis_script.R*) to perform the splice-aware LINDA network reconstructions for the KD and Ctrl conditions. Here we also provide functionalities (*/LINDA_Analysis/prepare_cytoscape_view.R*) which integrate the two networks and provide tables which can be given to [Cytoscape Application](https://cytoscape.org/) for a nice visualization of the resulting networks and the identified splicing effect mechanisms. For more details, please see documentation provided in */LINDA_Analysis/README.md*.

**IMPORTANT NOTE:** The users need to run the analysis pipeline according to the following order: *Differential Splicing* --> *Gene Expression* --> *LINDA Analysis*.

**TODO:** Add post-hoc network and enrichment analyses over the LINDA solutions that we got.
