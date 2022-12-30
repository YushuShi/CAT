# source("compDist.R")
# source("compDistList.R")
# source("paraCAT.R")
# source("paraCATMiRKAT.R")
CAT<-function(testList,otutable,taxonomy,metric="Weighted UniFrac",metaData,outcomeVar,tree=NULL,method="PERMANOVA",nperm=9999,parallel=TRUE,nCore=2){
  testResult<-rep(NA,length(testList))

  if(is.null(tree)){
    if(sum(metric%in%c("WeightUniFrac","Unweighted UniFrac","robust"))>0){
      stop("UniFrac distance needs a tree!")
    }
    if((sum(rownames(otutable)%in% rownames(taxonomy))!=nrow(taxonomy))|
       (sum(rownames(taxonomy) %in% rownames(otutable))!=nrow(otutable))){
      if(nrow(otutable)<nrow(taxonomy)){
        otutable<-otutable[rownames(otutable)%in% rownames(taxonomy),]
        print(paste0("Before padding: number of rows in OTU table: ",nrow(otutable),", and number of rows in taxonomy table: ",nrow(taxonomy)))
        padding<-matrix(0,ncol=ncol(otutable),nrow=nrow(taxonomy)-sum(rownames(taxonomy)%in% rownames(otutable)))
        rownames(padding)<-rownames(taxonomy)[!(rownames(taxonomy)%in% rownames(otutable))]
        colnames(padding)<-colnames(taxonomy)
        otutable<-rbind(otutable,padding)
        print(paste0("After padding: number of rows in OTU table: ",nrow(otutable),", and number of rows in taxonomy table: ",nrow(taxonomy)))
      }else{
        taxonomy<-taxonomy[rownames(taxonomy)%in% rownames(otutable),]
        print(paste0("Before padding: number of rows in OTU table: ",nrow(otutable),", and number of rows in taxonomy table: ",nrow(taxonomy)))
        padding<-matrix(0,ncol=ncol(taxonomy),nrow=nrow(otutable)-sum(rownames(otutable)%in% rownames(taxonomy)))
        rownames(padding)<-rownames(otutable)[!(rownames(otutable)%in% rownames(taxonomy))]
        colnames(padding)<-colnames(otutable)
        taxonomy<-rbind(taxonomy,padding)
        print(paste0("After padding: number of rows in OTU table: ",nrow(otutable),", and number of rows in taxonomy table: ",nrow(taxonomy)))
      }
      }
    }else{
    tree$root.edge<-0
    if((sum(rownames(otutable)%in% tree$tip.label)!=length(tree$tip.label))|
      (sum(tree$tip.label %in% rownames(otutable))!=nrow(otutable))){
      otutable<-otutable[rownames(otutable)%in% tree$tip.label,]
      print(paste0(nrow(otutable)," OTUs and ",length(tree$tip.label)," tip nodes."))
      if(nrow(otutable)<length(tree$tip.label)){
        padding<-matrix(0,ncol=ncol(otutable),nrow=length(tree$tip.label)-nrow(otutable))
        rownames(padding)<-tree$tip.label[!(tree$tip.label%in% rownames(otutable))]
        colnames(padding)<-colnames(otutable)
        otutable<-rbind(otutable,padding)
      }
    }
    if((sum(rownames(taxonomy)%in% tree$tip.label)!=length(tree$tip.label))|
       (sum(tree$tip.label %in% rownames(taxonomy))!=nrow(taxonomy))){
      taxonomy<-taxonomy[rownames(taxonomy)%in% tree$tip.label,]
      print(paste0(nrow(taxonomy)," rows in taxonomy table, and ",length(tree$tip.label)," tip nodes."))
      if(nrow(taxonomy)<length(tree$tip.label)){
        padding<-matrix(0,ncol=ncol(taxonomy),nrow=length(tree$tip.label)-sum(tree$tip.label%in% rownames(taxonomy)))
        rownames(padding)<-tree$tip.label[!(tree$tip.label%in% rownames(taxonomy))]
        colnames(padding)<-colnames(taxonomy)
        taxonomy<-rbind(taxonomy,padding)
      }
    }
    }
  otutable<-otutable[,rownames(metaData)]
  taxonomy<-taxonomy[rownames(otutable),]

  if(method=="PERMANOVA"){
    distMat<-compDist(otutable,metric,tree)
    suppressMessages(distResult<-adonis(distMat~metaData[,outcomeVar],permutations = nperm)$aov.tab)
    origpValue<-distResult[1,6]
    
    if(parallel){
      registerDoParallel(nCore)
      testResult<-foreach(seedNum=1:length(testList),.combine='c',
                          .packages=c("ape","vegan","GUniFrac")
      ) %dopar% paraCAT(seedNum,testList,otutable,taxonomy,metric,metaData,outcomeVar,nperm,
                        origpValue,tree)
      stopImplicitCluster()
    }else{
      for(j in 1:length(testList)){
        otutemp<-otutable
        otuUnder<-rownames(otutable)[apply(taxonomy, 1, function(x) sum(any(x %in% testList[j])))>0] 
        otutemp[otuUnder,]<-0
        distMat<-compDist(otutemp,metric,tree)
        suppressMessages(distResult<-adonis(distMat~metaData[,outcomeVar],permutations = nperm)$aov.tab)
        catpValue<-distResult[1,6]
        testResult[j]<-prop.test(c(origpValue*(nperm+1),catpValue*(nperm+1)),rep(nperm+1,2))$p.value
      }
    }
  }else{
    Ks<-compDistList(otutable,metric,tree)
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

    origpValue<-temp$omnibus_p 
    if(parallel){
      registerDoParallel(nCore)
      testResult<-foreach(seedNum=1:length(testList),.combine='c',
                          .packages=c("ape","vegan","GUniFrac","MiRKAT")
      ) %dopar% paraCATMiRKAT(seedNum,testList,otutable,taxonomy,metric,metaData,outcomeVar,nperm,
                        origpValue,tree)
      stopImplicitCluster()
    }else{
      for(j in 1:length(testList)){
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
        testResult[j]<-prop.test(c(origpValue*(nperm+1),catpValue*(nperm+1)),rep(nperm+1,2))$p.value
      }
    }
  }
  testResult
}
  


