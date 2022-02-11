**Differential Splicing**

Here we prepare the rMATS exon skipping data (*/Data/SE.MATS.JunctionCountOnly.txt*) for the KDvsCtrl comparison into a 
format which can be taken as an input and analyzed by LINDA. In */Data/SE.MATS.JunctionCountOnly.txt* we have been provided an effect size of the exon skipping
event (*IncLevelDifference*) as well as a significance for that event (*PValue/FDR*). However the locations in which the skipping events accur have been provided as
genomic coordinates and we need to transform them into corresponding Exon ID's. For that we rely on the steps provided in the *analysis_script.R* script.

**Analysis Steps:**
Step-by-Step documentation of the *analysis_script.R* script:

**[1-14]:** Calling of the needed libraries.

**[16]:** Reading the rMATS ES data for the KDvsCtrl comparison.

**[18-22]:** Obtaining genomic coordinates for each Ensembl Exon ID's.

**[24-25]:** Initializing the ES input object needed for the LINDA analysis (data-frame with 3 columns: *exon_id*, *effect* and *significance*).

**[27-67]:** Mapping genomic coordinates of */Data/SE.MATS.JunctionCountOnly.txt* into Exon ID's from Ensembl and transformin table into a data-frame. Those 
genomic coordinates that overlap with the coordinates for each individual exon provided by Ensembl will be assigned that specific Exon ID.

**[70-71]:** Inverting the siggne of the *effect* for the CtrlvsKD comparison. Those exons which have a negative *effect* (or rMATS *IncLevelDifference*) in the 
KdvsCtrl comarison are considered as skipped in the KD condition; while those exons which have a negative *effect* (or rMATS *IncLevelDifference*) in the 
CtrlvsKD comarison are considered as skipped in the Ctrl condition.

**[73-74]:** We save both data-frames in the */Data/* directory as we will need them for later analysis.
