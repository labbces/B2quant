library(tximport)
library(reshape2)
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

geneIDs_GALBA_dups<-read.delim("geneIDs_GALBA_dups.txt",header=T)
saveRDS(geneIDs_GALBA_dups$Name, 'geneIDs_GALBA_dups.rds')

geneIDs_BRAKER3_dups<-read.delim("geneIDs_BRAKER3_dups.txt",header=T)
saveRDS(geneIDs_BRAKER3_dups$Name, 'geneIDs_BRAKER3_dups.rds')
#targets<-targets[which(! targets$SampleName %in% c('SRR1979663')),]

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
dataTPM_BRAKER3_melt<-melt(dataTPM_BRAKER3[which(rowSums(dataTPM_BRAKER3)>0),])
colnames(dataTPM_BRAKER3_melt)<-c('Gene','Sample','TPM')
dataTPM_BRAKER3_melt$DevStage<-NA
dataTPM_BRAKER3_melt[which(dataTPM_BRAKER3_melt$Sample %in% c('SC_163P','SC_138P','SC_235P','SC_234P')),'DevStage']<-'P'
dataTPM_BRAKER3_melt[which(dataTPM_BRAKER3_melt$Sample %in% c('SC_234M','SC_163M','SC_138M','SC_235M')),'DevStage']<-'M'
dataTPM_BRAKER3_melt[which(dataTPM_BRAKER3_melt$Sample %in% c('SC_234B0','SC_235B0','SC_138B0')),'DevStage']<-'B0'
dataTPM_BRAKER3_melt[which(dataTPM_BRAKER3_melt$Sample %in% c('SC_234B','SC_163B','SC_138B')),'DevStage']<-'B'

# dataMeanTPM_BRAKER3<-cbind(rowMeans(dataTPM_BRAKER3[,which(targets$DevStage =='B')]),
#           rowMeans(dataTPM_BRAKER3[,which(targets$DevStage =='B0')]),
#           rowMeans(dataTPM_BRAKER3[,which(targets$DevStage =='M')]),
#           rowMeans(dataTPM_BRAKER3[,which(targets$DevStage =='P')])
# )

saveRDS(dataTPM_BRAKER3_melt, 'dataTPM_BRAKER3.rds')

myFiles_BRAKER3_dups<-paste(wd, 'data/', targets$SampleName,"_salmon_dups_BRAKER3/quant.sf",sep="")
names(myFiles_BRAKER3_dups)<-targets$LibraryName
all(file.exists(myFiles_BRAKER3_dups))
txi.salmon_BRAKER3_dups<-tximport(files = myFiles_BRAKER3_dups, type = 'salmon', tx2gene = tx2gene_BRAKER3_dups, txIn = TRUE, txOut = FALSE)
head(txi.salmon_BRAKER3_dups$counts)
dataTPM_BRAKER3_dups<-txi.salmon_BRAKER3_dups$abundance
dataTPM_BRAKER3_dups_melt<-melt(dataTPM_BRAKER3_dups[which(rowSums(dataTPM_BRAKER3_dups)>0),])
colnames(dataTPM_BRAKER3_dups_melt)<-c('Gene','Sample','TPM')
dataTPM_BRAKER3_dups_melt$DevStage<-NA
dataTPM_BRAKER3_dups_melt[which(dataTPM_BRAKER3_dups_melt$Sample %in% c('SC_163P','SC_138P','SC_235P','SC_234P')),'DevStage']<-'P'
dataTPM_BRAKER3_dups_melt[which(dataTPM_BRAKER3_dups_melt$Sample %in% c('SC_234M','SC_163M','SC_138M','SC_235M')),'DevStage']<-'M'
dataTPM_BRAKER3_dups_melt[which(dataTPM_BRAKER3_dups_melt$Sample %in% c('SC_234B0','SC_235B0','SC_138B0')),'DevStage']<-'B0'
dataTPM_BRAKER3_dups_melt[which(dataTPM_BRAKER3_dups_melt$Sample %in% c('SC_234B','SC_163B','SC_138B')),'DevStage']<-'B'
saveRDS(dataTPM_BRAKER3_dups_melt, 'dataTPM_BRAKER3_dups.rds')




myFiles_GALBA<-paste(wd, 'data/', targets$SampleName,"_salmon_GALBA/quant.sf",sep="")
names(myFiles_GALBA)<-targets$LibraryName
all(file.exists(myFiles_GALBA))
txi.salmon_GALBA<-tximport(files = myFiles_GALBA, type = 'salmon', tx2gene = tx2gene_GALBA, txIn = TRUE, txOut = FALSE)
head(txi.salmon_GALBA$counts)
dataTPM_GALBA<-txi.salmon_GALBA$abundance
dataTPM_GALBA_melt<-melt(dataTPM_GALBA[which(rowSums(dataTPM_GALBA)>0),])
colnames(dataTPM_GALBA_melt)<-c('Gene','Sample','TPM')
dataTPM_GALBA_melt$DevStage<-NA
dataTPM_GALBA_melt[which(dataTPM_GALBA_melt$Sample %in% c('SC_163P','SC_138P','SC_235P','SC_234P')),'DevStage']<-'P'
dataTPM_GALBA_melt[which(dataTPM_GALBA_melt$Sample %in% c('SC_234M','SC_163M','SC_138M','SC_235M')),'DevStage']<-'M'
dataTPM_GALBA_melt[which(dataTPM_GALBA_melt$Sample %in% c('SC_234B0','SC_235B0','SC_138B0')),'DevStage']<-'B0'
dataTPM_GALBA_melt[which(dataTPM_GALBA_melt$Sample %in% c('SC_234B','SC_163B','SC_138B')),'DevStage']<-'B'
saveRDS(dataTPM_GALBA_melt, 'dataTPM_GALBA.rds')

myFiles_GALBA_dups<-paste(wd, 'data/', targets$SampleName,"_salmon_dups_GALBA/quant.sf",sep="")
names(myFiles_GALBA_dups)<-targets$LibraryName
all(file.exists(myFiles_GALBA_dups))
txi.salmon_GALBA_dups<-tximport(files = myFiles_GALBA_dups, type = 'salmon', tx2gene = tx2gene_GALBA_dups, txIn = TRUE, txOut = FALSE)
head(txi.salmon_GALBA_dups$counts)
dataTPM_GALBA_dups<-txi.salmon_GALBA_dups$abundance
dataTPM_GALBA_dups_melt<-melt(dataTPM_GALBA_dups[which(rowSums(dataTPM_GALBA_dups)>0),])
colnames(dataTPM_GALBA_dups_melt)<-c('Gene','Sample','TPM')
dataTPM_GALBA_dups_melt$DevStage<-NA
dataTPM_GALBA_dups_melt[which(dataTPM_GALBA_dups_melt$Sample %in% c('SC_163P','SC_138P','SC_235P','SC_234P')),'DevStage']<-'P'
dataTPM_GALBA_dups_melt[which(dataTPM_GALBA_dups_melt$Sample %in% c('SC_234M','SC_163M','SC_138M','SC_235M')),'DevStage']<-'M'
dataTPM_GALBA_dups_melt[which(dataTPM_GALBA_dups_melt$Sample %in% c('SC_234B0','SC_235B0','SC_138B0')),'DevStage']<-'B0'
dataTPM_GALBA_dups_melt[which(dataTPM_GALBA_dups_melt$Sample %in% c('SC_234B','SC_163B','SC_138B')),'DevStage']<-'B'
saveRDS(dataTPM_GALBA_dups_melt, 'dataTPM_GALBA_dups.rds')