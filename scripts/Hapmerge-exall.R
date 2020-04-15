library(dplyr, warn.conflicts = FALSE)

args <- commandArgs(trailingOnly = TRUE)
popfile = args[1]
thresh= as.numeric(args[2])
ped1 = args[3]
pedfiles = args[4:length(args)]
population=read.table(popfile, header=TRUE)

d1=read.table(paste0(ped1 ,".vcf.ped"), colClasses = "character")
if (!is.na(pedfiles[1])) {
  for (ped in pedfiles) {
  d=read.table(paste0(ped ,".vcf.ped"), colClasses = "character")
  d1 = cbind(d1,d[,-c(1:2)])
  }
}
data = d1
Maternal=data[ , !c(TRUE,FALSE) ]
Paternal=data[ , !c(FALSE,TRUE) ]

colnames(Maternal)[1] <- "sample"
Maternalp=cbind(population,Maternal)		
colnames(Paternal)[1] <- "sample"
Paternalp=cbind(population,Paternal)		
Maternalp=merge(population[, c(1:3)],Maternal)
Paternalp=merge(population[, c(1:3)],Paternal)
Maternalp$mat <- c("_M")
Paternalp$pat <- c("_P")
vcf1=read.table(paste0(ped1, ".vcf"))

if (!is.na(pedfiles[1])) {
  for (ped in pedfiles) {
    vcf2=read.table(paste0(ped, ".vcf"))
    vcf1=rbind(vcf1,vcf2)
  }
}

vcf = vcf1
vcf.ref=t(vcf[,2:4])
vcf.ref=apply(format(vcf.ref), 2, paste, collapse="_")
colnames(Maternalp) <- c("pop", "superpop","sample",vcf.ref)
colnames(Paternalp) <- c("pop", "superpop","sample",vcf.ref)

Maternalp<- Maternalp[ , !grepl( "esv" , names( Maternalp ) ) ]
Maternalp.dd=Maternalp[,-c(1:3)]
if (!is.null(dim(Maternalp.dd))){
  dd <- apply( Maternalp.dd ,1, paste , collapse = "" )
  dd <- as.data.frame(dd)
  matern=cbind(Maternalp.dd=Maternalp,dd)
}else{
  matern=cbind(Maternalp.dd=Maternalp,Maternalp.dd)
}

Paternalp<- Paternalp[ , !grepl( "esv" , names( Paternalp ) ) ]
Paternalp.dd=Paternalp[,-c(1:3)]
if (!is.null(dim(Maternalp.dd))){
  dd <- apply( Paternalp.dd ,1, paste , collapse = "" )
  dd <- as.data.frame(dd)
  patern=cbind(Paternalp.dd=Paternalp,dd)
} else {
  patern=cbind(Paternalp.dd=Paternalp,Paternalp.dd)
}
vcf.ref=t(vcf[,2:4])
vcf.ref<- vcf.ref[ , !grepl( "esv" , vcf.ref[2,]) ]
if (!is.null(dim(vcf.ref))){
  vcf.ref=apply(format(vcf.ref), 2, paste, collapse="_")  
} else{
  vcf.ref=paste(as.array(format(vcf.ref)), collapse='_')
}

colnames(matern) <- c("sample", "superpop","pop",vcf.ref,"phase","dd")
colnames(patern) <- c("sample", "superpop","pop",vcf.ref,"phase","dd")


combined=rbind(matern,patern)
combined <- combined[, !duplicated(colnames(combined), fromLast = TRUE)] 
combined1 <- combined %>%
  dplyr::group_by(dd) %>%
  dplyr::summarise(test = toString(sample)) %>%
  dplyr::ungroup()

colnames(combined1) <- c ("dd","sample")

freq.all <- data.frame(table(combined$dd)) ##Retreive haplotype numbers for 2504 individuals all together
freq.pop <- as.data.frame.matrix(table(combined$dd,combined$pop)) ##Retreive haplotype numbers for populations
freq.suppop <- as.data.frame.matrix(table(combined$dd,combined$superpop)) ##Retreive haplotype numbers for superpopulations

######Representing the data
mergehapfreq <- function (x){
  dd <- rownames(x)
  rownames(x) <- NULL
  x1 <- cbind(dd, x)
  x1$dd <- gsub("_.", "", x1$dd)
  x1_1 <- dplyr::group_by(x1, dd) %>% dplyr::summarise_all(sum)
  df<-data.frame(x1_1)
  rownames(df) <- df[,1]
  df[,1] <- NULL
  freq.pop1=df[apply(df,1,function(x) !all(x<thresh)),] ###removing rows with all populations <0.005 AF
  dd <- rownames(freq.pop1)
  rownames(freq.pop1) <- NULL
  x2 <- cbind(dd, freq.pop1)
  return (x2)
}
uniqfreq <- mergehapfreq(freq.pop)
dd <- as.data.frame(rownames(freq.pop))
colnames(dd) <- c("dd")
dd1 <- dd %>% 
  tidyr::separate(dd, c("dd", "a"), extra='drop')
freq.pop <- cbind(dd1, freq.pop)
rownames(freq.pop) <- NULL

freq.pop1 <- subset(freq.pop, dd %in% uniqfreq$dd)
mp_col <- paste(freq.pop1$dd, freq.pop1$a, sep='_')
rownames(freq.pop1) <- mp_col
freq.pop1[, c(1:2)] <- NULL

op1=merge(freq.pop1,freq.suppop,by="row.names",all.x=TRUE)
op1_1=merge(op1,combined1,by.x="Row.names", by.y="dd",all.x=TRUE)
op2=merge(op1_1,freq.all, by.x="Row.names", by.y="Var1", all.x = TRUE)
comb = combined[!duplicated(combined$dd),]
op3=merge(op2,comb[,-c(1:3)],by.x="Row.names", by.y="dd",all.x=TRUE)

op4 = op3[, -c(33:length(op3))]
op5 = op3[, -c(1:32)]

if (!ncol(op5[vapply(op5, function(x) length(unique(x)) > 1, logical(1L))]) == 3){
  op6 <- cbind(op4, op5[vapply(op5, function(x) length(unique(x)) > 1, logical(1L))])
} else{
  op6 <- cbind(op4, op5)
}

aa <- dplyr::distinct(op6, Row.names, .keep_all = TRUE)
if (!is.null(aa)){
  write.table(format(op6,digits=3),paste0(sapply(strsplit(ped1, '-exon'), `[`, 1), "-Hap.xls"), quote=FALSE, sep="\t", row.names = FALSE)
  tryCatch({
    aa=aa[order(aa$Freq, decreasing = TRUE), ]
    write.table(format(aa,digits=3),paste0(sapply(strsplit(ped1, '-exon'), `[`, 1), "-Hap.xls"), quote=FALSE, sep="\t", row.names = FALSE)
  },error = function(e) {})
} else{
  write.table(format(op6,digits=3),paste0(sapply(strsplit(ped1, '-exon'), `[`, 1), "-Hap.xls"), quote=FALSE, sep="\t", row.names = FALSE)
}


