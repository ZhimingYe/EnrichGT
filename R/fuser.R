setGeneric("cpres_internal_getter",function(x,...) standardGeneric("cpres_internal_getter"))
setMethod("cpres_internal_getter", signature(x = "enrichResult"),function(x,...){
  y<-x@result
  if(dim(y)[1]==0){
    stop("no enrichment result contains")
  }
  if(sum(colnames(y)%in%c("ID","Description","GeneRatio","pvalue","p.adjust","geneID","Count"))!=7){
    stop("At list contains needed columns")
  }
  y<-y |> dplyr::select(ID,Description,GeneRatio,pvalue,p.adjust,geneID,Count)
  return(y)
})
setMethod("cpres_internal_getter", signature(x = "compareClusterResult"),function(x,...){
  y<-x@compareClusterResult
  if(dim(y)[1]==0){
    stop("no enrichment result contains")
  }
  if(sum(colnames(y)%in%c("ID","Description","NES","pvalue","p.adjust","core_enrichment"))!=6){
    stop("At list contains needed columns")
  }
  y<-y |> dplyr::select(ID,Description,NES,pvalue,p.adjust,core_enrichment)
  return(y)
})
setMethod("cpres_internal_getter", signature(x = "data.frame"),function(x,...){
  if(sum(colnames(x)=="NES")==0){
    if(sum(colnames(y)%in%c("ID","Description","GeneRatio","pvalue","p.adjust","geneID","Count"))!=7){
      stop("At list contains needed columns")
    }
    if(dim(y)[1]==0){
      stop("no enrichment result contains")
    }
    y<-y |> dplyr::select(ID,Description,GeneRatio,pvalue,p.adjust,geneID,Count)
    return(y)
  }
  else{
    if(sum(colnames(y)%in%c("ID","Description","NES","pvalue","p.adjust","core_enrichment"))!=6){
      stop("At list contains needed columns")
    }
    if(dim(y)[1]==0){
      stop("no enrichment result contains")
    }
    y<-y |> dplyr::select(ID,Description,NES,pvalue,p.adjust,core_enrichment)
    return(y)
  }

})

