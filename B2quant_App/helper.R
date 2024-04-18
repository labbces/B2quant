getExpressionProfile<-function(data, minCPM, libGroup, pctLibs, libraries){
  keep<-(rowSums(data != 0)/ncol(data)*100)>=pctLibs
  expresedpct<-data[keep,]
  txpExpresed<-expresedpct[rowSums(as.data.frame(expresedpct[,libraries[which(libraries$LibraryTissue == libGroup),
                                                           'Library']] > minCPM))>0,]
  return(txpExpresed)
}
getNumberTranscripts<-function(data, minCPM, libGroup, pctLibs, libraries){
  keep<-(rowSums(data != 0)/ncol(data)*100)>=pctLibs
  expresedpct<-data[keep,]
  ntxpExpresed<-nrow(expresedpct[rowSums(as.data.frame(expresedpct[,libraries[which(libraries$LibraryTissue == libGroup),
                                                                'Library']] > minCPM))>0,])
  return(ntxpExpresed)
}

getTranscriptNames<-function(data, minCPM, libGroup, pctLibs, libraries){
  keep<-(rowSums(data != 0)/ncol(data)*100)>=pctLibs
  expresedpct<-data[keep,]
  txpNames<-rownames(expresedpct[rowSums(as.data.frame(expresedpct[,libraries[which(libraries$LibraryTissue == libGroup),
                                                                'Library']] > minCPM))>0,])
  return(txpNames)
}

getGeneNames<-function(data){
  geneNames<-unique(data)
  return(geneNames)
}