################################################################
#### Script to store functions used throughout the analysis ####
################################################################

library(plyr)

# Function to collect the total proportion of a chain across te whole population
# of cells

chain.prop <- function(files, keywords, annotation = NULL){
  cur_files <- lapply(as.list(files), function(n){
    if(all(sapply(keywords, grepl, n))){
      read.table(n, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    }
  })
  # Remove empty slots
  library.names <- tail(sapply(files, function(n){unlist(strsplit(n, "\\/"))}), n=1)
  library.names <- library.names[!unlist(lapply(cur_files, is.null))]
  cur_files <- cur_files[!unlist(lapply(cur_files, is.null))]
  
  # Collect the clone fractions for variable clones
  cur_files <- lapply(cur_files, function(n){
    n <- data.frame(n[,c("cloneFraction", "allVHitsWithScore")])
    n$allVHitsWithScore <- as.factor(as.character(sapply(n$allVHitsWithScore, 
                                                         function(x){unlist(strsplit(x, "\\*"))[1]})))
    n
  })
  
  # Sum across clone faction
  sum.list <- lapply(cur_files, function(n){ddply(n, "allVHitsWithScore", summarise, sum = sum(cloneFraction))})
  
  # Sum across similar annotations if annotation df is givern
  if(!is.null(annotation)){
    sum.list <- lapply(sum.list, function(n){
        n$chain <- annotation[as.character(n$allVHitsWithScore),1]
        n <- ddply(n, "chain", summarise, sum = sum(sum))
        colnames(n) <- c("allVHitsWithScore", "sum")
        n
        })
  }
  
  # Merge into one dataframe
  mat <- matrix(data = NA, ncol = length(cur_files), 
                nrow = length(unique(unlist(lapply(sum.list, function(n){n$allVHitsWithScore})))))
  rownames(mat) <- as.character(unique(unlist(lapply(sum.list, function(n){n$allVHitsWithScore}))))
  colnames(mat) <- library.names
  
  for(i in 1:ncol(mat)){
    cur_n <- sum.list[[i]]
    mat[as.character(cur_n$allVHitsWithScore),i] <- cur_n$sum
  }
  mat[is.na(mat)] <- 0
  
  mat
}

#### CDR3 analysis
topCDR3 <- function(files, top, keywords, select = NULL,
                    annotation = NULL){
  cur_files <- lapply(as.list(files), function(n){
    if(all(sapply(keywords, grepl, n))){
      read.table(n, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    }
  })
  # Remove empty slots
  library.names <- tail(sapply(files, function(n){unlist(strsplit(n, "\\/"))}), n=1)
  library.names <- library.names[!unlist(lapply(cur_files, is.null))]
  cur_files <- cur_files[!unlist(lapply(cur_files, is.null))]
  
  # Collect clone fractions
  cur_files <- lapply(cur_files, function(n){
    n <- data.frame(n[,c("cloneFraction", "allVHitsWithScore", "clonalSequence", "aaSeqCDR3")])
    n$allVHitsWithScore <- as.factor(as.character(sapply(n$allVHitsWithScore, function(x){unlist(strsplit(x, "\\*"))[1]})))
    n
  })
  
  # Exclude certain V chains if select != NULL
  if(!is.null(select)){
    cur_files <- lapply(cur_files, function(n){
      n[n$allVHitsWithScore == select,]
    })
  }
  
  # Return the top expanded clones
  # Collect 
  cur_out.na <- unlist(lapply(cur_files, function(n){
    if(length(n$clonalSequence) > 1){
      n$clonalSequence[1:top]
    }
    else{
      n$clonalSequence
    }
  }))
  cur_out.na.unique <- unique(cur_out.na)
  
  cur_out.aa <- unlist(lapply(cur_files, function(n){
    if(length(n$clonalSequence) > 1){
      n$aaSeqCDR3[1:top]
    }
    else{
      n$aaSeqCDR3
    }
  }))
  
  cur_out.chain <- unlist(lapply(cur_files, function(n){
    if(length(n$clonalSequence) > 1){
      as.character(n$allVHitsWithScore[1:top])
    }
    else{
      as.character(n$allVHitsWithScore)
    }
  }))
  
  # Reannotate chains
  if(!is.null(annotation)){
    cur_out.chain <- as.character(annotation[cur_out.chain,])
  }
  
  
  data.frame(na = cur_out.na.unique[!is.na(cur_out.na.unique)],
             aa = cur_out.aa[match(cur_out.na.unique[!is.na(cur_out.na.unique)],
                                         cur_out.na)],
             chain = cur_out.chain[match(cur_out.na.unique[!is.na(cur_out.na.unique)],
                                         cur_out.na)], stringsAsFactors = FALSE)
}

# Function to collect the clones fraction for given clones
clones.fraction <- function(files, name, clones){
  keywords <- unlist(strsplit(name, "_"))[1:4]
  select <- unlist(strsplit(name, "_"))[5]
  replicate <- unlist(strsplit(name, "_"))[6]
  
  cur_files <- lapply(as.list(files), function(n){
    if(all(sapply(keywords, grepl, n))){
      read.table(n, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    }
  })
  # Remove empty slots
  library.names <- tail(sapply(files, function(n){unlist(strsplit(n, "\\/"))}), n=1)
  library.names <- library.names[!unlist(lapply(cur_files, is.null))]
  cur_files <- cur_files[!unlist(lapply(cur_files, is.null))]
  
  # Collect clone fractions
  cur_files <- lapply(cur_files, function(n){
    n <- data.frame(n[,c("cloneFraction", "allVHitsWithScore", "clonalSequence")])
    n$allVHitsWithScore <- as.factor(as.character(sapply(n$allVHitsWithScore, function(x){unlist(strsplit(x, "\\*"))[1]})))
    n
  })
  
  # Exclude certain V chains if select != X
  if(!(select == "X")){
    cur_files <- lapply(cur_files, function(n){
        n[grepl(select, n$allVHitsWithScore),]
    })
  }
  
  cur_rep <- cur_files[[as.numeric(replicate)]]
  
  out_vector <- rep(0, length(clones))
  m <- match(cur_rep$clonalSequence, clones)
  out_vector[m[!is.na(m)]] <- cur_rep$cloneFraction[!is.na(m)]
  out_vector
}

#### Function to compute shared clones
shared_clones <- function(files, in.names, all.reads = FALSE, subsample = NULL,
                          select = NULL){
  cur_files <- lapply(as.list(files), function(n){
      read.table(n, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  })
  
  # Exclude certain V chains if select != NULL
  if(!is.null(select)){
    cur_files <- lapply(cur_files, function(n){
      n[grepl(select, n$allVHitsWithScore),]
    })
  }
  
  mat.out <- matrix(data = NA, ncol = length(cur_files), 
                    nrow = length(cur_files))
  colnames(mat.out) <- rownames(mat.out) <- in.names
  
  if(all.reads){
    for(i in 1:length(cur_files)){
      for(j in 1:length(cur_files)){
        if(j <= i){
          x <- rep(cur_files[[i]]$clonalSequence, 
                   times = cur_files[[i]]$cloneCount)
          y <- rep(cur_files[[j]]$clonalSequence,
                   times = cur_files[[j]]$cloneCount)
          if(!is.null(subsample)){
            set.seed(12345)
            x <- x[sample(1:length(x), subsample)]
            set.seed(12345)
            y <- y[sample(1:length(y), subsample)]
          }
          
          mat.out[i,j] <- length(rep(sort(intersect(x, y)), 
                                       pmin(table(x[x %in% y]), table(y[y %in% x]))))

        }
        else{
          x <- rep(cur_files[[i]]$clonalSequence, 
                   times = cur_files[[i]]$cloneCount)
          y <- rep(cur_files[[j]]$clonalSequence,
                   times = cur_files[[j]]$cloneCount)
          if(!is.null(subsample)){
            set.seed(12345)
            x <- x[sample(1:length(x), subsample)]
            set.seed(12345)
            y <- y[sample(1:length(y), subsample)]
          }
          cur_inter <- length(rep(sort(intersect(x, y)), 
                                     pmin(table(x[x %in% y]), table(y[y %in% x]))))

          mat.out[i,j] <- cur_inter/(length(x) + 
                                       length(y) -
                                       cur_inter)
        }
      }
    }
  }
  else{
    for(i in 1:length(cur_files)){
      for(j in 1:length(cur_files)){
        if(j <= i){
          mat.out[i,j] <- length(intersect(cur_files[[i]]$clonalSequence,
                                           cur_files[[j]]$clonalSequence))
        }
        else{
          cur_inter <- length(intersect(cur_files[[i]]$clonalSequence,
                                        cur_files[[j]]$clonalSequence))
          mat.out[i,j] <- cur_inter/(length(cur_files[[i]]$clonalSequence) + 
                                       length(cur_files[[j]]$clonalSequence) -
                                       cur_inter)
        }
      }
    }
  }
  
  mat.out
}

#### Clonal diversity analysis
clone_diversity <- function(files, in.names, all.reads = FALSE, 
                            subsample = NULL, select = NULL){
  cur_files <- lapply(as.list(files), function(n){
    read.table(n, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  })
  
  # Exclude certain V chains if select != NULL
  if(!is.null(select)){
    cur_files <- lapply(cur_files, function(n){
      n[grepl(select, n$allVHitsWithScore),]
    })
  }
  
  cur_files <- lapply(cur_files, function(n){
    n <- data.frame(n[,c("allVHitsWithScore", "cloneCount", "cloneFraction", "clonalSequence", "aaSeqCDR3")])
    n$allVHitsWithScore <- sapply(n$allVHitsWithScore, function(x){unlist(strsplit(x, "\\*"))[1]})
    n
  })
  
  mat <- matrix(data = NA, nrow = 5, ncol = length(cur_files))
  colnames(mat) <- in.names
  rownames(mat) <- c("Highly expanded: > 2%", "Expanded: 1% - 1.99%",
                     "Frequent: 0.5% - 0.99", "Infrequent: < 0.5%",
                     "Singletons")
  
  if(all.reads){
    mat["Singletons",] <- unlist(lapply(cur_files, function(n){
      sum(n$cloneCount[n$cloneCount == 1])/sum(n$cloneCount)
    }))
    mat["Infrequent: < 0.5%",] <- unlist(lapply(cur_files, function(n){
      sum(n$cloneCount[n$cloneFraction < 0.005 & n$cloneCount > 1])/sum(n$cloneCount)
    }))
    mat["Frequent: 0.5% - 0.99",] <- unlist(lapply(cur_files, function(n){
      sum(n$cloneCount[n$cloneFraction >= 0.005 & n$cloneFraction < 0.01])/sum(n$cloneCount)
    }))
    mat["Expanded: 1% - 1.99%",] <- unlist(lapply(cur_files, function(n){
      sum(n$cloneCount[n$cloneFraction >= 0.01 & n$cloneFraction < 0.02])/sum(n$cloneCount)
    }))
    mat["Highly expanded: > 2%",] <- unlist(lapply(cur_files, function(n){
      sum(n$cloneCount[n$cloneFraction >= 0.02])/sum(n$cloneCount)
    }))
  }
  else{
    mat["Singletons",] <- unlist(lapply(cur_files, function(n){
      sum(n$cloneCount == 1)/nrow(n)
    }))
    mat["Infrequent: < 0.5%",] <- unlist(lapply(cur_files, function(n){
      sum(n$cloneFraction < 0.005 & n$cloneCount > 1)/nrow(n)
    }))
    mat["Frequent: 0.5% - 0.99",] <- unlist(lapply(cur_files, function(n){
      sum(n$cloneFraction >= 0.005 & n$cloneFraction < 0.01)/nrow(n)
    }))
    mat["Expanded: 1% - 1.99%",] <- unlist(lapply(cur_files, function(n){
      sum(n$cloneFraction >= 0.01 & n$cloneFraction < 0.02)/nrow(n)
    }))
    mat["Highly expanded: > 2%",] <- unlist(lapply(cur_files, function(n){
      sum(n$cloneFraction >= 0.02)/nrow(n)
    }))
  }
  mat
}

