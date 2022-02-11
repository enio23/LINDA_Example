**Differential Gene Expression**

Here we perform Differential Gene Expression (DGE) Analysis (*/Data/gene_expression_data.csv*) over gene expression data for the KDvsCtrl comparison. From DGE's
we then estimate differential enzmatic activities of Transcription Factor (TF) targets from [DoRothEA](https://github.com/saezlab/dorothea) and from where we 
estimate the most significantly active TF's. This step is important for the LINDA analyses since it aims to reconstruct regulatory protein interaction networks 
upstream of significantly regulated TF's.

**Analysis Steps:**
Step-by-Step documentation of the *analysis_script.R* script:

**[1-15]:** Calling of the needed libraries.

**[17]:** Sourcing the *estimate_significance.R* function. This is a function used to estimate the significance of the estimated TF activity scores based on a
permutation analysis procedure where the TF-to-Gene relationships from DoRothEA are randomized multiple times.

**[19-20]:** Loading tables which map *ensembl_gene_id* into *external_gene_name* from Ensembl.

**[22-24]:** Loading the TF-to-Gene relations from DoRothEA for the A, B and C classes (for more information please refer [here](https://github.com/saezlab/dorothea)).

**[26-72]:** Downloading the Control-HepG2 and U2AF2-KD-HepG2 shRNA data (two replicates each) from ENCORE and converting Gene ID's to Gene Names.

**[74-116]:** Building the FPKM matrix from the downloaded data.

**[118-123]:** Identifying and filtering Genes with zero counts and save them for filtering purposes in the later steps.

**[125-140]:** Performing DGE analysis over the KDvsCtrl comparison.

**[142-146]:** Estimating TF activity scores and significance from DGE analysis data and saving the results in */Data/*.

