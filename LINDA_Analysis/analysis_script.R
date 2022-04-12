#! /usr/bin/env Rscript

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
library("aggregation")
library("LINDA")
library("fgsea")

dir.create("output")

## Load data
load(file = "../Gene_Expression/output/ttop_u2af1_hepg2.RData")
load(file = "../Gene_Expression/output/tf_act_u2af1_hepg2.RData")
load(file = "../Transcript_Expression/output/splice_effect_u2af1_hepg2.RData")

tf_act$nes <- 0
tf_act$nes[which(tf_act$pval<=0.05)] <- 1
tf_act <- tf_act[, c(1, 2)]

load(file = system.file("extdata", "digger_human_transcripts.RData", package = "LINDA"))

idx2keep <- intersect(x = which(bg$gene_source%in%ttop$ID), 
                      y = which(bg$gene_target%in%ttop$ID))
bg <- bg[idx2keep, ]

##Set the LINDA optimization parameters
lambda1 <- 100
lambda2 <- 1
input.node <- NULL
mipgap = 0
relgap = 0
populate = 500
nSolutions = 100
intensity = 1
timelimit = 7200
process_log = FALSE
replace = 1
solverPath <- "/beegfs/homes/egjerga/cplex"
pValThresh <- 0.05
threads <- 20
top <- length(which(tf_act$nes>=1))



## Perform analysis without considering for splicing effects
## Here the spliced mechanisms upon U2AF2 Knockdown of HepG2 cells are identified
res_ctrl <- runLINDA(input.scores = tf_act, 
                     as.input = NULL, 
                     background.network = bg, 
                     solverPath = solverPath, 
                     input.node = input.node, 
                     pValThresh = pValThresh, 
                     splice_effect_sign = "both", 
                     top = top, 
                     lambda1 = lambda1, 
                     lambda2 = lambda2, 
                     mipgap = mipgap, 
                     relgap = relgap, 
                     populate = populate, 
                     nSolutions = nSolutions, 
                     timelimit = timelimit, 
                     intensity = intensity, 
                     replace = replace, 
                     threads = threads, 
                     solver = "cplex", 
                     condition = 1)

save(res_ctrl, file = "output/res_ctrl.RData")



## Perform analysis while considering for splicing effects
## Here the spliced mechanisms upon U2AF2 Knockdown of HepG2 cells are identified
res_kd <- runLINDA(input.scores = tf_act, 
                   as.input = splice_effect, 
                   background.network = bg, 
                   solverPath = solverPath, 
                   input.node = input.node, 
                   pValThresh = pValThresh, 
                   splice_effect_sign = "both", 
                   top = top, 
                   lambda1 = lambda1, 
                   lambda2 = lambda2, 
                   mipgap = mipgap, 
                   relgap = relgap, 
                   populate = populate, 
                   nSolutions = nSolutions, 
                   timelimit = timelimit, 
                   intensity = intensity, 
                   replace = replace, 
                   threads = threads, 
                   solver = "cplex", 
                   condition = 2)

save(res_kd, file = "output/res_kd.RData")



## Prepare CytoScape visualization
source("prepare_cytoscape_visualization.R")

net <- prepare_cytoscape_visualization(netA = res_ctrl$combined_interactions, 
                                       netB = res_kd$combined_interactions, 
                                       spliceA = NULL, 
                                       spliceB = splice_effect, 
                                       pValThresh = 0.05, 
                                       background.network = bg, 
                                       tf = tf_act$id[which(tf_act$nes==1)], 
                                       sources = "Perturbation")

write.table(x = net$network, file = "output/network.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = net$attributes, file = "output/attributes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
