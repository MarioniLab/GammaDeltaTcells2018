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
topCDR3 <- function(files, top, keywords, select = NULL){
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
      n$clonalSequence
    }
  }))
  data.frame(na = cur_out.na.unique[!is.na(cur_out.na.unique)],
                   aa = cur_out.aa[match(cur_out.na.unique[!is.na(cur_out.na.unique)],
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
  out_vector[clones %in% cur_rep$clonalSequence] <- cur_rep$cloneFraction[cur_rep$clonalSequence %in% clones[clones %in% cur_rep$clonalSequence]]
  out_vector
}

#### Function to compute shared clones
shared_clones <- function(files, in.names){
  cur_files <- lapply(as.list(files), function(n){
      read.table(n, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  })
  
  mat.out <- matrix(data = NA, ncol = length(cur_files), 
                    nrow = length(cur_files))
  colnames(mat.out) <- rownames(mat.out) <- in.names
  
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
  mat.out
}

