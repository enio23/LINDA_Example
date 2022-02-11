set.seed(1234)

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
library("LINDA")

dir.create("output")

## Load files
load(file = "../Data/missing_genes.RData")
load(file = "../Data/diff_tf_act.RData")
load(file = "../Data/kd_vs_ctrl.RData")
load(file = "../Data/ctrl_vs_kd.RData")

input.scores$nes <- 0
input.scores$nes[which(input.scores$pval<=0.05)] <- 1
input.scores <- input.scores[, c(1, 2)]

## Load the default DIGGER Background Network from LINDA and filter it from the non-expressed genes
load(file = system.file("extdata", "digger_human_exons.RData", package = "LINDA"))

idx2rem <- unique(c(which(bg$gene_source%in%missingGenes), which(bg$gene_target%in%missingGenes)))
bg <- bg[-idx2rem, ]

##Set the LINDA optimization parameters
lambda1 <- 100
lambda2 <- 1
input.node <- NULL
mipgap = 0
relgap = 0
populate = 500
nSolutions = 100
intensity = 1
timelimit = 3600
process_log = FALSE
replace = 1
solverPath <- "/home/enio/Downloads/cplex"
pValThresh <- 0.1
threads <- 4
top <- length(which(input.scores$nes>=1))

## Perform analysis for the Kd_vs_Ctrl comparison
## Here the spliced mechanisms upon U2AF2 Knockdown of HepG2 cells are identified
res <- runLINDA(input.scores = input.scores, as.input = kd_vs_ctrl, background.network = bg, 
                solverPath = solverPath, input.node = input.node, pValThresh = pValThresh, 
                top = top, lambda1 = lambda1, lambda2 = lambda2, mipgap = mipgap, 
                relgap = relgap, populate = populate, nSolutions = nSolutions, 
                timelimit = timelimit, intensity = intensity, replace = replace, 
                threads = threads, solver = "cplex")

save(res, file = "output/linda_res_u2af2_kd_vs_ctrl.RData")
write.table(x = res$combined_interactions, file = "output/linda_res_u2af2_kd_vs_ctrl.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

## Perform analysis for the Ctrl_vs_KD comparison
## Here the spliced mechanisms on the Control HepG2 cells are identified
res <- runLINDA(input.scores = input.scores, as.input = ctrl_vs_kd, background.network = bg, 
                solverPath = solverPath, input.node = input.node, pValThresh = pValThresh, 
                top = top, lambda1 = lambda1, lambda2 = lambda2, mipgap = mipgap, 
                relgap = relgap, populate = populate, nSolutions = nSolutions, 
                timelimit = timelimit, intensity = intensity, replace = replace, 
                threads = threads, solver = "cplex")

save(res, file = "output/linda_res_u2af2_ctrl_vs_kd.RData")
write.table(x = res$combined_interactions, file = "output/linda_res_u2af2_ctrl_vs_kd.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
