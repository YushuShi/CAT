paraCAT<-function(j,testList,otutable,taxonomy,metric,metaData,outcomeVar,nperm,origpValue,tree){
  otutemp<-otutable
  otuUnder<-rownames(otutable)[apply(taxonomy, 1, function(x) sum(any(x %in% testList[j])))>0] 
  otutemp[otuUnder,]<-0
  distMat<-compDist(otutemp,metric,tree)
  distResult<-adonis(distMat~metaData[,outcomeVar],permutations = nperm)$aov.tab
  catpValue<-distResult[1,6]
  prop.test(c(origpValue*(nperm+1),catpValue*(nperm+1)),rep(nperm+1,2))$p.value
}