args <- commandArgs(trailingOnly = TRUE)

popfile = args[1]
ped1 = args[2]
ped2 = args[3]
ped3 = args[4]
thresh = args[5]
population=read.table(popfile, header=TRUE)

d1=read.table(paste0(ped1 ,".vcf.ped"), colClasses = "character")
d2=read.table(paste0(ped2, ".vcf.ped"), colClasses = "character")
d3=read.table(paste0(ped3, ".vcf.ped"), colClasses = "character")
data=cbind(d1,d2[,-c(1:2)])
data=cbind(data, d3[,-c(1:2)])

Maternal=data[ , !c(TRUE,FALSE) ]
Paternal=data[ , !c(FALSE,TRUE) ]

#Maternalp=cbind(population,Maternal)
#Paternalp=cbind(population,Paternal)
colnames(Maternal)[1] <- "sample"
Maternalp=cbind(population,Maternal)		
colnames(Paternal)[1] <- "sample"
Paternalp=cbind(population,Paternal)		
Maternalp=merge(population[, c(1:3)],Maternal)
Paternalp=merge(population[, c(1:3)],Paternal)

vcf1=read.table(paste0(ped1, ".vcf"))
vcf2=read.table(paste0(ped2, ".vcf"))
vcf3=read.table(paste0(ped3, ".vcf"))
vcf=rbind(vcf1,vcf2)
vcf=rbind(vcf, vcf3)

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
#colnames(matern) <- c("pop", "superpop","sampleId",vcf.ref,"dd")
#colnames(patern) <- c("pop", "superpop","sampleId",vcf.ref,"dd")
colnames(matern) <- c("sample", "superpop","pop",vcf.ref,"dd")
colnames(patern) <- c("sample", "superpop","pop",vcf.ref,"dd")


combined=rbind(matern,patern)
#print(combined)
freq.all <- data.frame(table(combined$dd)) ##Retreive haplotype numbers for 2504 individuals all together
freq.pop <- as.data.frame.matrix(table(combined$dd,combined$pop)) ##Retreive haplotype numbers for populations
freq.suppop <- as.data.frame.matrix(table(combined$dd,combined$superpop)) ##Retreive haplotype numbers for superpopulations

######Representing the data
freq.pop1=freq.pop[apply(freq.pop,1,function(x) !all(x<thresh)),] ###removing rows with all populations <0.005 AF

op1=merge(freq.pop1,freq.suppop,by="row.names",all.x=TRUE)
op2=merge(op1,freq.all, by.x="Row.names", by.y="Var1", all.x = TRUE)
comb = combined[!duplicated(combined$dd),]
op3=merge(op2,comb[,-c(1:3)],by.x="Row.names", by.y="dd",all.x=TRUE)

aa=op3[vapply(op3, function(x) length(unique(x)) > 1, logical(1L))]
if (!is.null(aa)){
  write.table(format(op3,digits=3),paste0(sapply(strsplit(ped1, '-exon'), `[`, 1), "-Hap.xls"), quote=FALSE, sep="\t", row.names = FALSE)
  #  write.table(format(op3,digits=3),paste0(gsub("-exon1-", "", ped1), "-Hap.xls"), quote=FALSE, sep="\t", row.names = FALSE)
  tryCatch({
    aa=aa[order(aa$Freq, decreasing = TRUE), ]
    write.table(format(aa,digits=3),paste0(sapply(strsplit(ped1, '-exon'), `[`, 1), "-Hap.xls"), quote=FALSE, sep="\t", row.names = FALSE)
  },error = function(e) {})
    #write.table(format(aa,digits=3),paste0(gsub("-exon1-", "", ped1), "-Hap.xls"), quote=FALSE, sep="\t", row.names = FALSE)
} else{
  write.table(format(op3,digits=3),paste0(sapply(strsplit(ped1, '-exon'), `[`, 1), "-Hap.xls"), quote=FALSE, sep="\t", row.names = FALSE)
  #write.table(format(op3,digits=3),paste0(gsub("-exon1-", "", ped1), "-Hap.xls"), quote=FALSE, sep="\t", row.names = FALSE)
}

