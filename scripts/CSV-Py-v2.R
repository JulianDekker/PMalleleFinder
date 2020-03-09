#setwd("C:/Users/julian/Documents/bafstu/SHARK")
args = commandArgs(trailingOnly = T)
path= args[1]
name= args[2]
setwd(path)

myfun <- function(filename) {
  d1 <- read.table(filename, header = TRUE, check.names=FALSE, sep="\t")
  d1[,sapply(d1,class) == "logical"] <-
    sapply(d1[,sapply(d1,class) == "logical"],
           function(i) substr(as.character(i),1,1)) #read.table reads "T" as "TRUE" so this will convert it back to just "T"
  print(d1)
  d2<-d1[,-c(1,7:33)]
  d3=data.frame(d2)
  colnames(d3) <- sub(" .*", "", colnames(d3))
  newfile <- gsub(".xls", ".csv", filename)
  newfile1 <- gsub("Hap", "py", newfile)
  write.csv(d3[,-c(1:5)],newfile1, row.names = FALSE, quote=FALSE)
  d4 = as.data.frame(d3)
  colnames(d4) <- sub("_.*", "", colnames(d4))
  colnames(d4) <- sub("X", "", colnames(d4))
  print(d4[,-c(1:5), drop=F])
  write.csv(d4[,-c(1:5), drop=F],newfile1, row.names = FALSE, quote=FALSE)
}

files = list.files(pattern=paste0(name, "-Hap.xls"))
#files = list.files(pattern="TRAV22-201-exon1-Hap.xls")
invisible(lapply(files, myfun))

