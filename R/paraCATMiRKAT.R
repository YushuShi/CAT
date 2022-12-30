paraCATMiRKAT<-function(j,testList,otutable,taxonomy,metric,metaData,outcomeVar,nperm,origpValue,tree){
  otutemp<-otutable
  otuUnder<-rownames(otutable)[apply(taxonomy, 1, function(x) sum(any(x %in% testList[j])))>0] 
  otutemp[otuUnder,]<-0
  Ks<-compDistList(otutemp,metric,tree)
  temp<-NA
  if(length(outcomeVar)>1){
    temp<-MiRKATS(obstime = metaData[,outcomeVar[1]], delta = metaData[,outcomeVar[2]], Ks = Ks,
                  perm = TRUE, omnibus="permutation", nperm=nperm)
  }else{
    if(length(table(metaData[,outcomeVar]))==2){
      temp<-MiRKAT(metaData[,outcomeVar], Ks = Ks, method="permutation",
                   out_type="D",
                          omnibus="permutation", nperm=nperm)        
    }else{
      temp<-MiRKAT(metaData[,outcomeVar], Ks = Ks,
                   omnibus="permutation", nperm=nperm)  
    }
  }
  catpValue<-temp$omnibus_p 
  prop.test(c(origpValue*(nperm+1),catpValue*(nperm+1)),rep(nperm+1,2))$p.value
}