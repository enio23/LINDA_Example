prepare_cytoscape_visualization <- function(netA = netA,
                                            netB = netB,
                                            spliceA = NULL,
                                            spliceB = NULL,
                                            pValThresh = 0.05,
                                            background.network = background.network,
                                            tf = tf,
                                            sources = "Perturbation",
                                            network_path = NULL,
                                            attributes_path = NULL){
  
  intA <- paste0(netA[, 1], "=", netA[, 3])
  intB <- paste0(netB[, 1], "=", netB[, 3])
  wwA <- as.numeric(netA[, 2])/max(as.numeric(netA[, 2]))
  wwB <- as.numeric(netB[, 2])/max(as.numeric(netB[, 2]))
  ints <- unique(c(intA, intB))
  common <- intersect(x = intA, y = intB)
  
  combined <- matrix(data = , nrow = length(ints), ncol = 6)
  combined[, 1] <- sapply(strsplit(x = ints, split = "=", fixed = TRUE), "[", 1)
  combined[, 2] <- sapply(strsplit(x = ints, split = "=", fixed = TRUE), "[", 2)
  
  for(ii in 1:length(ints)){
    
    idx1 <- which(intA==ints[ii])
    idx2 <- which(intB==ints[ii])
    if(ints[ii]%in%common){
      
      combined[ii, 3] <- "Common"
      reacs <- unique(c(unlist(strsplit(x = netA[idx1, 4], 
                                        split = "; ", 
                                        fixed = TRUE)),
                        unlist(strsplit(x = netB[idx2, 4], 
                                        split = "; ", 
                                        fixed = TRUE))))
      combined[ii, 4] <- paste0(reacs, collapse = "; ")
      combined[ii, 5] <- wwA[idx1] + wwB[idx2]
      combined[ii, 6] <- wwA[idx1] - wwB[idx2]
      
    } else {
      
      if(length(idx1)==1){
        
        combined[ii, 3] <- "Network_A"
        reacs <- unique(unlist(strsplit(x = netA[idx1, 4], 
                                        split = "; ", 
                                        fixed = TRUE)))
        combined[ii, 4] <- paste0(reacs, collapse = "; ")
        combined[ii, 5] <- wwA[idx1]
        combined[ii, 6] <- 1
        
      } else {
        
        combined[ii, 3] <- "Network_B"
        reacs <- unique(unlist(strsplit(x = netB[idx2, 4], 
                                        split = "; ", 
                                        fixed = TRUE)))
        combined[ii, 4] <- paste0(reacs, collapse = "; ")
        combined[ii, 5] <- wwB[idx2]
        combined[ii, 6] <- -1
        
      }
      
    }
    
  }
  
  network <- matrix(data = , nrow = 1, ncol = 6)
  for(ii in 1:nrow(combined)){
    
    reacs <- unique(unlist(strsplit(x = combined[ii, 4], 
                                    split = "; ", 
                                    fixed = TRUE)))
    
    for(jj in 1:length(reacs)){
      
      sGene <- combined[ii, 1]
      tGene <- combined[ii, 2]
      sDomain <- strsplit(x = reacs[jj], split = "=", fixed = TRUE)[[1]][1]
      tDomain <- strsplit(x = reacs[jj], split = "=", fixed = TRUE)[[1]][2]
      
      toBind <- matrix(data = , nrow = 3, ncol = 6)
      
      toBind[, 1] <- c(sGene, 
                       paste0(sGene, "_", sDomain), 
                       paste0(tGene, "_", tDomain))
      
      toBind[, 2] <- c(paste0(sGene, "_", sDomain), 
                       paste0(tGene, "_", tDomain), 
                       tGene)
      
      toBind[, 3] <- combined[ii, 3]
      
      toBind[, 4] <- "normal"
      
      toBind[, 5] <- combined[ii, 5]
      
      toBind[, 6] <- combined[ii, 6]
      
      network <- unique(rbind(network, toBind))
      
    }
    
  }
  network <- network[2:nrow(network), ]
  
  bgReacS <- paste0(background.network$gene_source, 
                    "_", 
                    background.network$pfam_source)
  
  bgReacT <- paste0(background.network$gene_target, 
                    "_", 
                    background.network$pfam_target)
  
  splicedExonsA <- spliceA$id[
    intersect(x = which(spliceA$effect<0), 
              y = which(spliceA$significance<=pValThresh))]
  
  splicedExonsB <- spliceB$id[
    intersect(x = which(spliceB$effect<0), 
              y = which(spliceB$significance<=pValThresh))]
  
  splicedExons <- unique(c(splicedExonsA, splicedExonsB))
  
  nodes <- unique(c(network[, 1], network[, 2]))
  attributes <- matrix(data = , nrow = length(nodes), ncol = 4)
  attributes[, 1] <- nodes
  attributes[, 2] <- "P"
  attributes[which(grepl(pattern = "_", x = nodes, fixed = TRUE)), 2] <- "D"
  attributes[which(nodes%in%tf), 2] <- "T"
  attributes[, 3] <- "normal"
  
  for(ii in 1:nrow(attributes)){
    
    if(attributes[ii, 2]=="D"){
      
      exons <- unique(c(unlist(
        strsplit(
          x = background.network$exon_source[
            which(bgReacS==attributes[ii, 1])], 
          split = "_", 
          fixed = TRUE)),
        unlist(
          strsplit(
            x = background.network$exon_target[
              which(bgReacT==attributes[ii, 1])], 
            split = "_", 
            fixed = TRUE))))
      
      ss <- intersect(x = exons, y = splicedExons)
      
      if(length(ss)>0){
        
        attributes[ii, 3] <- "spliced"
        attributes[ii, 4] <- paste0(ss, collapse = "_")
        
      }
      
    }
    
  }
  
  idx <- which(attributes[, 4]!="")
  if(length(idx)>0){
    
    sDomains <- attributes[idx, 1]
    
    for(jj in 1:length(sDomains)){
      
      idx1 <- which(attributes[, 1]==sDomains[jj])
      attributes[idx1, 1] <- paste0(attributes[idx1, 1], 
                                    "(", 
                                    attributes[idx1, 4], 
                                    ")")
      
      idx2 <- which(network[, 1]==sDomains[jj])
      if(length(idx2)>0){
        
        network[idx2, 1] <- attributes[idx1, 1]
        
        network[idx2, 4] <- "spliced"
        
      }
      
      idx3 <- which(network[, 2]==sDomains[jj])
      if(length(idx3)>0){
        
        network[idx3, 2] <- attributes[idx1, 1]
        
        network[idx3, 4] <- "spliced"
        
      }
      
    }
    
  }
  
  idx <- intersect(x = which(network[, 4]=="spliced"), 
                   y = c(which(!grepl(pattern = "_", x = network[, 1], fixed = TRUE)),
                         which(!grepl(pattern = "_", x = network[, 2], fixed = TRUE))))
  
  if(length(idx)>0){
    
    network[idx, 4] <- "normal"
    
  }
  
  colnames(network) <- c("source", "target", "interaction", 
                         "splice_status", "weight", "shade")
  
  attributes <- attributes[, 1:3]
  colnames(attributes) <- c("species", "type", "status")
  
  reac <- paste0(network[, 1], "=", network[, 2])
  ind <- which(duplicated(reac))

  if(length(ind)>0){

    uReac <- unique(reac)
    net <- matrix(data = , nrow = length(uReac), ncol = ncol(network))
    colnames(net) <- colnames(network)

    for(ii in 1:length(uReac)){

      idx1 <- which(reac==uReac[ii])
      idx2 <- idx1[which(as.numeric(network[idx1, 5])==max(as.numeric(network[idx1, 5])))[1]]

      net[ii, 1] <- strsplit(x = uReac[ii], split = "=", fixed = TRUE)[[1]][1]
      net[ii, 2] <- strsplit(x = uReac[ii], split = "=", fixed = TRUE)[[1]][2]

      if(length(intersect(x = "Common", y = network[idx, 3]))>0 ||
         length(intersect(x = c("Network_A", "Network_B"), y = network[idx, 3]))==2){

        net[ii, 3] <- "Common"

      } else {

        net[ii, 3] <- network[idx2, 3]

      }

      net[ii, 4:6] <- network[idx2, 4:6]

    }

  }
  
  status <- rep("solid", nrow(net))
  idx <- intersect(x = which(grepl(pattern = "_", x = net[, 1], fixed = TRUE)), 
                   y = which(grepl(pattern = "_", x = net[, 2], fixed = TRUE)))
  status[idx] <- "arrow"
  status <- as.matrix(status)

  net <- cbind(net, status)
  colnames(net)[ncol(net)] <- "status"
  
  idx <- intersect(x = which(!grepl(pattern = "_", x = net[, 1], fixed = TRUE)), 
                   y = which(grepl(pattern = "_", x = net[, 2], fixed = TRUE)))
  
  for(ii in 1:length(idx)){
    
    ind <- intersect(x = which(net[, 1]==net[idx[ii], 2]), 
                     y = which(net[, 2]==net[idx[ii], 1]))
    
    if(length(ind)==1){
      
      net[ind, 1:ncol(net)] <- net[idx[ii], 1:ncol(net)]
      
    }
    
  }
  
  idx <- which(attributes[, 1]%in%sources)
  if(length(idx)>0){
    
    attributes[idx, 2] <- "S"
    
  }
  
  res <- list()
  res[[length(res)+1]] <- unique(net)
  res[[length(res)+1]] <- attributes
  names(res) <- c("network", "attributes")
  
  if(!is.null(network_path)){
    
    write.table(x = net, file = network_path, quote = FALSE, 
                sep = "\t", row.names = FALSE, col.names = TRUE)
    
  }
  
  if(!is.null(attributes_path)){
    
    write.table(x = attributes, file = attributes_path, quote = FALSE, 
                sep = "\t", row.names = FALSE, col.names = TRUE)
    
  }
  
  return(res)
  
}
