library("readxl")
library("readr")
library("XML")
library("igraph")
library("foreach")
library("doParallel")
library("ranger")
library("palmerpenguins")
library("tidyverse")
library("kableExtra")
library("biomaRt")
library("EnsDb.Hsapiens.v86")
library("EnsDb.Hsapiens.v75")
library("ape")

as <- read.delim("../Data/SE.MATS.JunctionCountOnly.txt")

## Map coordinates to Exon ID's
mart=useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="apr2019.archive.ensembl.org")
ensg2symbol=getBM(attributes=c("ensembl_gene_id","external_gene_name", 
                               "ensembl_exon_id", "exon_chrom_start",
                               "exon_chrom_end"), mart=mart)

kd_vs_ctrl <- matrix(data = , nrow = 1, ncol = 3)
colnames(kd_vs_ctrl) <- c("exon_id", "effect", "significance")

for(ii in 1:nrow(as)){
  
  print(paste0("Step ---- ", ii, "/", nrow(as)))
  
  # Find the genomic region for each gene present in the ES table
  idx <- which(ensg2symbol$ensembl_gene_id==strsplit(x = as$GeneID[ii], split = ".", fixed = TRUE)[[1]][1])
  
  temp <- ensg2symbol[idx, ]
  
  # Find which genomic coordinates overlap with the current coordinates under consideration in the ES table
  rng <- matrix(data = , nrow = nrow(temp), ncol = 4)
  rng[, 1] <- temp$exon_chrom_start
  rng[, 2] <- temp$exon_chrom_end
  rng[, 3] <- as$exonStart_0base[ii]
  rng[, 4] <- as$exonEnd[ii]
  
  olap = which((rng[,1] <= rng[,4]) & (rng[,2] >= rng[,3]))
  
  if(length(olap)>0){
    
    toBind <- matrix(data = , nrow = length(olap), ncol = 3)
    for(jj in 1:length(olap)){
      
      toBind[jj, 1] <- temp$ensembl_exon_id[olap[jj]]
      toBind[jj, 2] <- as$IncLevelDifference[ii]
      toBind[jj, 3] <- as$PValue[ii]
      
    }
    
    colnames(toBind) <- c("exon_id", "effect", "significance")
    
    kd_vs_ctrl <- unique(rbind(kd_vs_ctrl, toBind))
    
  }
  
}
kd_vs_ctrl <- kd_vs_ctrl[2:nrow(kd_vs_ctrl), ]

kd_vs_ctrl <- as.data.frame(kd_vs_ctrl)
kd_vs_ctrl$effect <- as.numeric(kd_vs_ctrl$effect)
kd_vs_ctrl$significance <- as.numeric(kd_vs_ctrl$significance)

## Flip the effect signs for the CtrlvsKD comparison and save both files
ctrl_vs_kd <- kd_vs_ctrl
ctrl_vs_kd$effect <- (-1)*ctrl_vs_kd$effect

save(kd_vs_ctrl, file = "../Data/kd_vs_ctrl.RData")
save(ctrl_vs_kd, file = "../Data/ctrl_vs_kd.RData")
