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

dir.create("output")

source("estimate_significance.R")

mart=useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="apr2019.archive.ensembl.org")
ensg2symbol=getBM(attributes=c("ensembl_gene_id","external_gene_name"), mart=mart)

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE88nnn/GSE88226/matrix/GSE88226_series_matrix.txt.gz", 
              destfile = paste0(getwd(), "/temp_kd.txt"))
download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE88nnn/GSE88002/matrix/GSE88002_series_matrix.txt.gz", 
              destfile = paste0(getwd(), "/temp_ctrl.txt"))

temp_kd <- read.table(file = paste0("temp_kd.txt"), sep = " ")
temp_ctrl <- read.table(file = paste0("temp_ctrl.txt"), sep = " ")

# Download KD
experiments <- temp_kd$V1[which(grepl(pattern = "gene_quantifications_GRCh38.tsv.gz", x = temp_kd$V1, fixed = TRUE))]
experiments <- gsub(pattern = "\t", replacement = " ", x = experiments)
experiments <- unlist(x = strsplit(x = experiments, split = " ", fixed = TRUE))
experiments <- experiments[which(grepl(pattern = ".tsv", x = experiments, fixed = TRUE))]
experiments <- experiments[which(grepl(pattern = "GRCh38", x = experiments, fixed = TRUE))]

cnt <- 1
for(jj in 1:length(experiments)){
  if(grepl(pattern = "gene_quantifications", x = experiments[jj])){
    download.file(url = experiments[jj], destfile = paste0(getwd(), "/temp.tsv"))
    temp <- read_tsv(file = paste0(getwd(), "/temp.tsv"))
    if("gene_id"%in%colnames(temp)){
      write_tsv(x = temp, file = paste0("kd", cnt, ".tsv"))
      cnt <- cnt + 1
    }
    file.remove(paste0(getwd(), "/temp.tsv"))
  }
}

# Download Ctrl
experiments <- temp_ctrl$V1[which(grepl(pattern = "gene_quantifications_GRCh38.tsv.gz", x = temp_ctrl$V1, fixed = TRUE))]
experiments <- gsub(pattern = "\t", replacement = " ", x = experiments)
experiments <- unlist(x = strsplit(x = experiments, split = " ", fixed = TRUE))
experiments <- experiments[which(grepl(pattern = ".tsv", x = experiments, fixed = TRUE))]
experiments <- experiments[which(grepl(pattern = "GRCh38", x = experiments, fixed = TRUE))]
cnt <- 1
for(jj in 1:length(experiments)){
  if(grepl(pattern = "gene_quantifications", x = experiments[jj])){
    download.file(url = experiments[jj], destfile = paste0(getwd(), "/temp.tsv"))
    temp <- read_tsv(file = paste0(getwd(), "/temp.tsv"))
    if("gene_id"%in%colnames(temp)){
      write_tsv(x = temp, file = paste0("ctrl", cnt, ".tsv"))
      cnt <- cnt + 1
    }
    file.remove(paste0(getwd(), "/temp.tsv"))
  }
}

# Build FPKM data matrix
kd1 <- read.table(file = paste0("kd1.tsv"), header = TRUE)
kd2 <- read.table(file = paste0("kd2.tsv"), header = TRUE)
ctrl1 <- read.table(file = paste0("ctrl1.tsv"), header = TRUE)
ctrl2 <- read.table(file = paste0("ctrl2.tsv"), header = TRUE)

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
colnames(df) <- c(paste0("kd_", 1:2), paste0("ctrl_", 1:2))

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

dfMapped <- dfMapped[complete.cases(dfMapped), ]

# Do differential gene expression analysis
conditions<-factor(c("kd", "kd", "ctrl", "ctrl"))
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

save(ttop, file = paste0("output/ttop_u2af1_hepg2.RData"))

file.remove(paste0("kd1.tsv"))
file.remove(paste0("kd2.tsv"))
file.remove(paste0("ctrl1.tsv"))
file.remove(paste0("ctrl2.tsv"))
file.remove(paste0("temp_kd.txt"))
file.remove(paste0("temp_ctrl.txt"))

stats <- ttop$logFC
names(stats) <- ttop$ID

tf_act <- estimate_significance(expr = as.matrix(stats), regulons = regulons, nperm = 1000)
save(tf_act, file = "output/tf_act_u2af1_hepg2.RData")
