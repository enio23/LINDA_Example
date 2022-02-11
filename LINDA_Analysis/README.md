**LINDA Analysis**

Here we perform the network reconstruction analysis with LINDA. LINDA takes as inpu the TF activities data estimated in */Gene_Expression/* as well as the results
from the differential splicing analysis with rMATS (*/Differential_Splicing/*). The TF activities are neeed to infer the protein interaction mechanisms upstream of
the most significantly active TF's while the results from the differential splicing analysis are used to identify the splicing effects on protein interactions. We
perform two splice-aware network inference analysis where we take into consideration the splicing effects at the Control and Knockdown condition.

For that we rely on the steps provided in the *analysis_script.R* script.

**Analysis Steps:**
Step-by-Step documentation of the *analysis_script.R* script:

**[1-13]:** Calling of the needed libraries.

**[18-21]:** Loading the needed inputs for the analysis: the estimated differential TF activities (*diff_tf_act.RData*); differential splicing data at the Ctrl
(*ctrl_vs_kd.RData*) and KD (*kd_vs_ctrl.RData*) conditions; list of missing/non-expressed genes that we have obtained from */Gene_EXpression/analysis_script.R*.

**[23-25]:** We manually assign higher scores those TF's which pass the significance threshold (p-value<=0.05 on this case). However this step is not necessary if
you wish to deifne active TF's only based on their absolute TF activity scores. For this you can just leave the TF activity scores with the same value as they were
inferred by DoRothEA and just define an arbitrary number of the most regulated TF's based on the absolute activity value through the *top* parameter.

**[27-31]:** Filtering from DIGGER those interactions involving proteins whose corresponding genes have not been expressed or are missing. This step is not
mandatory however it is recomended as we filter the background network from redundant interactions and alleviate identifiability issues.

**[33-48]:** Setting the LINDA optimization parameters.

**[50-57]:** Here we perform the KDvsCtrl splice-aware analysis from where we can identify the spliced mechanisms upon U2AF2 Knockdown of HepG2 cells and how they
affecting the protein interaction networks in the Knockdown condition.

**[63-74]:** Here we perform the CtrlvsKD splice-aware analysis from where we can identify the spliced mechanisms on the Control HepG2 cells.

**TODO:** Update something about the results once the *prepare_cytoscape_view.R* function has been impelemnted within the LINDA R-Package.
