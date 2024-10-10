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
setMethod("cpres_internal_getter", signature(x = "gseaResult"),function(x,...){
  y<-x@result
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

#' 2-Group Comparison of enrichment results and further clustering and visualizing
#'
#' @param obj.test the enriched object from tested group. WARNING: `obj.test` and `obj.ctrl` should come from same database (e.g. GO Biological Process(GOBP)).
#' @param obj.ctrl the enriched object from control group. WARNING: `obj.test` and `obj.ctrl` should come from same database (e.g. GO Biological Process(GOBP)).
#' @param name.test optional, the name of the testing group. If is `NULL`, the object name of `obj.test` will be used.
#' @param name.ctrl optional, the name of the control group. If is `NULL`, the object name of `obj.ctrl` will be used.
#' @param ClusterNum how many cluster will be clustered
#' @param P.adj p.adjust cut-off. To avoid slow visualization, you can make stricter p-cut off.
#' @param force ignore all auto-self-checks, which is useful
#' @param nTop keep n top items according to p-adj in each cluster.
#' @param ... Others options.
#'
#' @details
#' Execute `obj.test` VS `obj.ctrl` tests, showing pathway overlaps (or differences) and meta-gene modules of test group and control group. Supports ORA and GSEA results (enriched object or data.frame). !WARNING!: `obj.test` and `obj.ctrl` should come from same database (e.g. GO Biological Process(GOBP)).
#' @return `List` containing multiple `EnrichGT_obj` objects. The `List` contains objects with overlapped enriched terms, unique enrich terms. In each slot, slot `enriched_result` contains a data.frame with enriched results. `gt_object` contains `gt` object. you can use `obj@gt_object` to get it and use functions from `gt` like `gtsave`. `gene_modules` is a list containing meta-gene modules of each cluster. `pathway_clusters` contains pathways names in each cluster.
#' @export
#'
#' @author Zhiming Ye
compareGT<-function(obj.test,obj.ctrl,name.test=NULL,name.ctrl=NULL,ClusterNum=15,P.adj=0.05,force=F,nTop=10,...){
  testdf<-cpres_internal_getter(obj.test)
  ctrldf<-cpres_internal_getter(obj.ctrl)
  overlapped <- testdf$ID[testdf$ID%in%ctrldf$ID]
  testOnly <- testdf$ID[!testdf$ID%in%ctrldf$ID]
  ctrlOnly <- ctrldf$ID[!ctrldf$ID%in%testdf$ID]
  ctrlOnlyDF<-ctrldf[ctrldf$ID%in%ctrlOnly,]
  testOnlyDF<-testdf[testdf$ID%in%testOnly,]
  ctrlOverlapDF<-ctrldf[!ctrldf$ID%in%ctrlOnly,]
  testOverlapDF<-testdf[!testdf$ID%in%testOnly,]
  Overlap_Control<-universalGT(ctrlOverlapDF,objname2="Overlap_Control")
  Overlap_Test<-universalGT(testOverlapDF,objname2="Overlap_Test")
  Control_Only<-universalGT(ctrlOnlyDF,objname2="Control_Only")
  Test_Only<-universalGT(testOnlyDF,objname2="Test_Only")
  resultlist<-list(`Overlap_Control`=Overlap_Control,`Overlap_Test`=Overlap_Test,`Control_Only`=Control_Only,`Test_Only`=Test_Only)
  return(resultlist)
}

universalGT<-function(x,...){
  if("NES"%in%colnames(x)){
    x1<-.genGSEAGT(x,...)
  }
  else{
    x1<-.genGT(x,...)
  }
  return(x1)
}
