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