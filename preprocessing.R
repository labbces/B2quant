library(tximport)
rm(list=ls())

wd<-"/data/diriano/B2quant/"
setwd(wd)

targets<-read.delim("targets.txt",header=T)
rownames(targets)<-targets$SampleName
colnames(targets)
targets$Individual<-as.factor(targets$Individual)
targets$DevStage<-as.factor(targets$DevStage)
targets$LibraryName<-as.factor(targets$LibraryName)
targets
targets<-targets[which(! targets$SampleName %in% c('SRR1979663')),]

#Getting tx2gene files, we have four collections of genes two for BRAKER and two for GALBA. There are two of each of these, because one includes duplicates and the other do not.
#here we will carry analyses for each of these.
tx2gene_BRAKER3<-read.delim("data/BRAKER3.tx2gene",header=F)
tx2gene_BRAKER3_dups<-read.delim("data/BRAKER3_dups.tx2gene",header=F)
tx2gene_GALBA<-read.delim("data/GALBA.tx2gene",header=F)
tx2gene_GALBA_dups<-read.delim("data/GALBA_dups.tx2gene",header=F)

#Preparing quantification files to load
myFiles_BRAKER3<-paste(wd, 'data/', targets$SampleName,"_salmon_BRAKER3/quant.sf",sep="")
names(myFiles_BRAKER3)<-targets$LibraryName
all(file.exists(myFiles_BRAKER3))
txi.salmon_BRAKER3<-tximport(files = myFiles_BRAKER3, type = 'salmon', tx2gene = tx2gene_BRAKER3, txIn = TRUE, txOut = FALSE)
head(txi.salmon_BRAKER3$counts)
dataTPM_BRAKER3<-txi.salmon_BRAKER3$abundance
# dataMeanTPM_BRAKER3<-cbind(rowMeans(dataTPM_BRAKER3[,which(targets$DevStage =='B')]),
#           rowMeans(dataTPM_BRAKER3[,which(targets$DevStage =='B0')]),
#           rowMeans(dataTPM_BRAKER3[,which(targets$DevStage =='M')]),
#           rowMeans(dataTPM_BRAKER3[,which(targets$DevStage =='P')])
# )

saveRDS(dataTPM_BRAKER3, 'dataTPM_BRAKER3.rds')

myFiles_BRAKER3_dups<-paste(wd, 'data/', targets$SampleName,"_salmon_dups_BRAKER3/quant.sf",sep="")
names(myFiles_BRAKER3_dups)<-targets$LibraryName
all(file.exists(myFiles_BRAKER3_dups))
txi.salmon_BRAKER3_dups<-tximport(files = myFiles_BRAKER3_dups, type = 'salmon', tx2gene = tx2gene_BRAKER3_dups, txIn = TRUE, txOut = FALSE)
head(txi.salmon_BRAKER3_dups$counts)
dataTPM_BRAKER3_dups<-txi.salmon_BRAKER3_dups$abundance
saveRDS(dataTPM_BRAKER3_dups, 'dataTPM_BRAKER3_dups.rds')




myFiles_GALBA<-paste(wd, 'data/', targets$SampleName,"_salmon_GALBA/quant.sf",sep="")
names(myFiles_GALBA)<-targets$LibraryName
all(file.exists(myFiles_GALBA))
txi.salmon_GALBA<-tximport(files = myFiles_GALBA, type = 'salmon', tx2gene = tx2gene_GALBA, txIn = TRUE, txOut = FALSE)
head(txi.salmon_GALBA$counts)
dataTPM_GALBA<-txi.salmon_GALBA$abundance
saveRDS(dataTPM_GALBA, 'dataTPM_GALBA.rds')

myFiles_GALBA_dups<-paste(wd, 'data/', targets$SampleName,"_salmon_dups_GALBA/quant.sf",sep="")
names(myFiles_GALBA_dups)<-targets$LibraryName
all(file.exists(myFiles_GALBA_dups))
txi.salmon_GALBA_dups<-tximport(files = myFiles_GALBA_dups, type = 'salmon', tx2gene = tx2gene_GALBA_dups, txIn = TRUE, txOut = FALSE)
head(txi.salmon_GALBA_dups$counts)
dataTPM_GALBA_dups<-txi.salmon_GALBA_dups$abundance
saveRDS(dataTPM_GALBA_dups, 'dataTPM_BRAKER3_dups.rds')
