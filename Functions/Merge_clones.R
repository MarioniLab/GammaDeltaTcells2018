library(openxlsx)

files <- list.files("Dropbox (Cambridge University)/GammaDelta/Analysis/Clones/Txt_files/", full.names = TRUE)
files.2 <- list.files("Dropbox (Cambridge University)/GammaDelta/Analysis/Clones/Txt_files/", full.names = FALSE)


for(i in seq(1, length(files), 6)){
  cur_data = list("TRD" = read.table(files[i+4], sep = "\t", header = TRUE),
                  "TRG" = read.table(files[i+5], sep = "\t", header = TRUE),
                  "TRA" = read.table(files[i+2], sep = "\t", header = TRUE),
                  "TRB" = read.table(files[i+3], sep = "\t", header = TRUE),
                  "IGH" = read.table(files[i], sep = "\t", header = TRUE),
                  "IGL" = read.table(files[i+1], sep = "\t", header = TRUE))
  file.name <- paste(unlist(strsplit(files.2[i], "_"))[1:11], collapse = "_")
  write.xlsx(cur_data, file = paste("Dropbox (Cambridge University)/GammaDelta/Analysis/Clones/Xlsx_files/", 
                                    file.name, ".xlsx", sep = ""))
}
