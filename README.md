# LINDA Example

This repository shows the steps needed to be performed by the user in order to run the LINDA analysis from the gene and transcript expression data. The example here, refers to studying the splicing effects over the protein interaction networks up U2AF1 Knockdown of HepG2 cells (data from [ENCORE project](https://www.encodeproject.org/encore-matrix/?type=Experiment&status=released&internal_tags=ENCORE)). 

The analysis follows the workflow depicted in the following graph:

<p align="center">
    <img src="https://github.com/enio23/LINDA_Example/blob/main/pipeline.jpeg" alt>
    <em> 1.The Gene and Transcript abundance data for two replicates (of the HepG2 cell-line) for the U2AF1-KD (GEO:GSE88226) and U2AF1-Ctrl (GEO:GSE88002) are downloaded and the expected counts are accesed. 2.Genes and Transcripts are filtered based on their expression. 3.With edgeR we perform differential analysis over Gene and Transcript abundance values for the KDvsCtrl comparison. 4. From the differential Gene Expression data we estimate TF activity scores by using the DoRothEA resource as well as a permutation test to estimate the significance of the activity scores (by randomizing Gene-To-TF relations 1000 times). 5. We filter DIGGER resource to only contain interactions between expressed genes. 6. We map the quantified Transcripts to Domain Pfam ID's and for each Domain we estimate an 'effect' score as the average logFC values of Transcript abundances and a 'significance' score after performing a Fissher aggregation of p-values of the mapping transcripts. 7. We give the three inputs to LINDA (TF's/Domins/DIGGER) in order to contextuaize regulatory signalling networks. </em>
</p>

The repository is organized with the following directories (and please follow the same order when re-analyzing the data):

**1. Gene Expression:** Here are provided the scripts ([/Gene_Expression/analysis_script.R](https://github.com/enio23/LINDA_Example/tree/main/Gene_Expression)) which have been used to estimate the differentially expressed genes (DGE's) for the KDvsCtrl comparison with [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) as well as TF activity estimation with [viper](https://www.bioconductor.org/packages/release/bioc/html/viper.html) and [DoRothEA](https://www.bioconductor.org/packages/release/bioc/html/viper.html). For more details, please see documentation provided in */Gene_Expression/README.md*.

**2. Transcript Expression:** Here are provided scripts ([/Transcript_Expression/analysis_script.R](https://github.com/enio23/LINDA_Example/tree/main/Transcript_Expression)) which have been used to assign an Exon ID to all the genomic coordinates in the identified exon skipping events (*/Data/SE.MATS.JunctionCountOnly.txt*).

**3. LINDA Analysis:** Here are provided the scripts ([/LINDA_Analysis/analysis_script.R](https://github.com/enio23/LINDA_Example/blob/main/LINDA_Analysis/analysis_script.R)) to perform the splice-aware and splice-unaware LINDA network reconstructions. Here we also provide functionalities ([/LINDA_Analysis/prepare_cytoscape_view.R](https://github.com/enio23/LINDA_Example/blob/main/LINDA_Analysis/prepare_cytoscape_visualization.R)) which integrate the two networks and provide tables which can be given to [Cytoscape Application](https://cytoscape.org/) for a nice visualization of the resulting networks and the identified splicing effect mechanisms.

**IMPORTANT NOTE:** For running this example, the users need to obtaine a [CPLEX License](https://www.ibm.com/products/ilog-cplex-optimization-studio) and set the *solverPath* parameter to the path where the CPLEX executable file has been stored.
