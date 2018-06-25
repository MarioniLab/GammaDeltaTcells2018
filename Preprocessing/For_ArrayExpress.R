library(openxlsx)

# Load all libraries

raw <- read.table("Dropbox (Cambridge University)/GammaDelta/Analysis/Data/GammaDelta_all.txt", sep = "\t")
norm <- read.table("Dropbox (Cambridge University)/GammaDelta/Analysis/Data/NormCounts.txt", sep = "\t")

meta <- read.csv("/Users/nils/Dropbox (Cambridge University)/GammaDelta/Analysis/Data/Metadata.csv", stringsAsFactors = FALSE)
meta <- meta[1:48,]

rownames(meta) <- paste("do", sub("DO", "", meta$DO.number), sep = "")

sdrf <- data.frame(libraries = colnames(raw),
                   age = sub(" mo", "", meta[colnames(raw),"age"]),
                   cell.type = sapply(meta[colnames(raw),"condtition"], function(n){unlist(strsplit(n, "-"))[1]}),
                   high.quality = meta$DO.number %in% sapply(colnames(norm), function(n){unlist(strsplit(n, "_"))[4]}))

write.xlsx(sdrf, "Dropbox (Cambridge University)/GammaDelta/Analysis/For_ArrayExpress/sdrf.xlsx")
