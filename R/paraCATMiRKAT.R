paraCATMiRKAT<-function(testIter,testList,otutable,taxonomy,metric,metaData,outcomeVar,numBS,origR2,origBS,tree){
  otutemp<-otutable
  otuUnder<-rownames(otutable)[apply(taxonomy, 1, function(x) sum(any(x %in% testList[testIter])))>0] 
  otutemp[otuUnder,]<-0
  Ks<-compDistList(otutemp,metric,tree)
  taxaR2<-MiRKATR2(metaData,outcomeVar,Ks)
  taxaBS<-rep(NA,numBS)
  for(BSiter in 1:numBS){
    indi<-sample(1:nrow(metaData),nrow(metaData),replace = TRUE)
    outcomeBS<-metaData[indi,outcomeVar]
    KsBS<-Ks
    temp<-NA
    for(k in 1:length(Ks)){
      for(i in 1:length(outcomeBS)){
        for(j in 1:length(outcomeBS)){
          KsBS[[k]][i,j]<-Ks[[k]][indi[i],indi[j]]
        }
      }
    }
    taxaBS[BSiter]<-MiRKATR2(outcomeBS,outcomeVar,KsBS)
  }
  BSvar<-var(taxaBS-origBS)
  BSpvalue<-pnorm((taxaR2-origR2)/sqrt(BSvar/nrow(metaData)))
  BSpvalue2Sided<-2*pnorm(-abs(taxaR2-origR2)/sqrt(BSvar/nrow(metaData)))
  c(BSpvalue,BSpvalue2Sided)
}