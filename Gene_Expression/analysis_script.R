#! /usr/bin/env Rscript

set.seed(1234)

library("readr")
library("org.Hs.eg.db")
library("biomaRt")
library("stringr")
library("dplyr")
library("edgeR")
library("dorothea")
library("foreach")
library("doParallel")
library("piano")
library("BiRewire")

source("estimate_significance.R")

mart=useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="apr2019.archive.ensembl.org")
ensg2symbol=getBM(attributes=c("ensembl_gene_id","external_gene_name"), mart=mart)

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

## Downloading the HepG2 Control and U2AF2-KD shRNA data from ENCORE
download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE88nnn/GSE88268/matrix/GSE88268_series_matrix.txt.gz", 
              destfile = "../Data/u2af2_kd.txt")
download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE88nnn/GSE88031/matrix/GSE88031_series_matrix.txt.gz", 
              destfile = "../Data/hepg2_ctrl.txt")


kd <- read.table(file = paste0("../Data/u2af2_kd.txt"), sep = " ")
ctrl <- read.table(file = paste0("../Data/hepg2_ctrl.txt"), sep = " ")

# Downloading the U2AF2-KD shRNA
experiments <- kd$V1[which(grepl(pattern = "gene_quantifications_GRCh38.tsv.gz", x = kd$V1, fixed = TRUE))]
experiments <- gsub(pattern = "\t", replacement = " ", x = experiments)
experiments <- unlist(x = strsplit(x = experiments, split = " ", fixed = TRUE))
experiments <- experiments[which(grepl(pattern = ".tsv", x = experiments, fixed = TRUE))]
experiments <- experiments[which(grepl(pattern = "GRCh38", x = experiments, fixed = TRUE))]
cnt <- 1
for(ii in 1:length(experiments)){
  if(grepl(pattern = "gene_quantifications", x = experiments[ii])){
    download.file(url = experiments[ii], destfile = paste0("../Data/temp", ii, ".tsv"))
    temp <- read_tsv(file = paste0("../Data/temp", ii, ".tsv"))
    if("gene_id"%in%colnames(temp)){
      write_tsv(x = temp, file = paste0("../Data/kd", cnt, ".tsv"))
      cnt <- cnt + 1
    }
    file.remove(paste0("../Data/temp", ii, ".tsv"))
  }
}

# Downloading the Ctrl-HepG2 shRNA
experiments <- ctrl$V1[which(grepl(pattern = "gene_quantifications_GRCh38.tsv.gz", x = ctrl$V1, fixed = TRUE))]
experiments <- gsub(pattern = "\t", replacement = " ", x = experiments)
experiments <- unlist(x = strsplit(x = experiments, split = " ", fixed = TRUE))
experiments <- experiments[which(grepl(pattern = ".tsv", x = experiments, fixed = TRUE))]
experiments <- experiments[which(grepl(pattern = "GRCh38", x = experiments, fixed = TRUE))]
cnt <- 1
for(ii in 1:length(experiments)){
  if(grepl(pattern = "gene_quantifications", x = experiments[ii])){
    download.file(url = experiments[ii], destfile = paste0("../Data/temp", ii, ".tsv"))
    temp <- read_tsv(file = paste0("../Data/temp", ii, ".tsv"))
    if("gene_id"%in%colnames(temp)){
      write_tsv(x = temp, file = paste0("../Data/ctrl", cnt, ".tsv"))
      cnt <- cnt + 1
    }
    file.remove(paste0("../Data/temp", ii, ".tsv"))
  }
}

# Building FPKM data matrix
kd1 <- read.table(file = "../Data/kd1.tsv", header = TRUE)
kd2 <- read.table(file = "../Data/kd2.tsv", header = TRUE)
ctrl1 <- read.table(file = "../Data/ctrl1.tsv", header = TRUE)
ctrl2 <- read.table(file = "../Data/ctrl2.tsv", header = TRUE)

commonGenes <- intersect(x = intersect(x = kd1$gene_id, y = kd2$gene_id), 
                         y = intersect(x = ctrl1$gene_id, y = ctrl2$gene_id))

kd1 <- kd1[which(kd1$gene_id%in%commonGenes), ]
kd2 <- kd2[which(kd2$gene_id%in%commonGenes), ]
ctrl1 <- ctrl1[which(ctrl1$gene_id%in%commonGenes), ]
ctrl2 <- ctrl2[which(ctrl2$gene_id%in%commonGenes), ]

kd1 <- kd1[order(kd1$gene_id), ]
kd2 <- kd2[order(kd2$gene_id), ]
ctrl1 <- ctrl1[order(ctrl1$gene_id), ]
ctrl2 <- ctrl2[order(ctrl2$gene_id), ]

idx2rem <- which(duplicated(kd1$gene_id)); if(length(idx2rem)>0){kd1 <- kd1[-idx2rem, ]}
idx2rem <- which(duplicated(kd2$gene_id)); if(length(idx2rem)>0){kd2 <- kd2[-idx2rem, ]}
idx2rem <- which(duplicated(ctrl1$gene_id)); if(length(idx2rem)>0){ctrl1 <- ctrl1[-idx2rem, ]}
idx2rem <- which(duplicated(ctrl2$gene_id)); if(length(idx2rem)>0){ctrl2 <- ctrl2[-idx2rem, ]}

df <- matrix(data = , nrow = nrow(kd1), ncol = 4)
df[, 1] <- kd1$FPKM
df[, 2] <- kd2$FPKM
df[, 3] <- ctrl1$FPKM
df[, 4] <- ctrl2$FPKM
df <- as.data.frame(df)
rownames(df) <- sapply(strsplit(x = kd1$gene_id, split = ".", fixed = TRUE), '[', 1)
colnames(df) <- c(paste0("ko_", 1:2), paste0("ctrl_", 1:2))

uGenes <- unique(ensg2symbol$external_gene_name)
dfMapped <- matrix(data = , nrow = length(uGenes), ncol = ncol(df))
for(ll in 1:length(uGenes)){
  ensg <- ensg2symbol$ensembl_gene_id[which(ensg2symbol$external_gene_name==uGenes[ll])]
  for(mm in 1:ncol(df)){
    dfMapped[ll, mm] <- median(x = df[which(rownames(df)%in%ensg), mm], na.rm = TRUE)
  }
}
rownames(dfMapped) <- uGenes
colnames(dfMapped) <- colnames(df)

## Identify and filter Genes with zero counts and save them for filtering in later steps purposes
missingGenes <- rownames(dfMapped)[which(rowMeans(dfMapped)==0)]
save(missingGenes, file = paste0("../Data/missing_genes.RData"))
df <- dfMapped[-which(rownames(dfMapped)%in%missingGenes), ]
df <- df[complete.cases(df), ]
write_csv(x = as.data.frame(df), file = "../Data/gene_expression_data.csv", col_names = TRUE)

## Differential gene expression analysis
conditions<-factor(c("ko", "ko", "ctrl", "ctrl"))
design <- model.matrix(~ conditions)

y <- DGEList(counts=dfMapped, group=conditions)
keep <- filterByExpr(y, group=conditions)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmQLFit(y,design,robust=TRUE)

res=glmQLFTest(fit, coef=2)
ttop=as.data.frame(topTags(res,n=nrow(df)))
ttop$ID <- rownames(ttop)

## TF Activities
ss <- ttop$logFC
names(ss) <- ttop$ID

input.scores <- estimate_significance(expr = as.matrix(ss), regulons = regulons, nperm = 1000)

save(input.scores, file = "../Data/diff_tf_act.RData")

file.remove("../Data/ctrl1.tsv")
file.remove("../Data/ctrl2.tsv")
file.remove("../Data/kd1.tsv")
file.remove("../Data/kd2.tsv")
