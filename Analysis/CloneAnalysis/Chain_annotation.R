all.files <- list.files("Dropbox (Cambridge University)/GammaDelta/Analysis/Clones/Txt_files/", full.names = TRUE)

# Collect all unique names

names <- unique(unlist(lapply(as.list(all.files), function(n){
  cur_table <- read.table(n, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cur_names <- cur_table[,"allVHitsWithScore"]
  as.character(sapply(cur_names, function(n){unlist(strsplit(n, "\\*"))[1]}))
})))

# Order
names <- names[order(names)]

# Write out
library(openxlsx)
write.xlsx(data.frame(names = names), "Dropbox (Cambridge University)/GammaDelta/Analysis/Clones/Variable_sequences.xlsx")
