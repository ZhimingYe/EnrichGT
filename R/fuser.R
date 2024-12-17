setGeneric("cpres_internal_getter",function(x,...) standardGeneric("cpres_internal_getter"))
setMethod("cpres_internal_getter", signature(x = "enrichResult"),function(x,...){
  y<-x@result
  if(dim(y)[1]==0){
    cli::cli_abort("no enrichment result contains")
  }
  if(sum(colnames(y)%in%c("ID","Description","GeneRatio","pvalue","p.adjust","geneID","Count"))!=7){
    cli::cli_abort("At list contains needed columns")
  }
  y<-y |> dplyr::select(ID,Description,GeneRatio,pvalue,p.adjust,geneID,Count)
  return(y)
})
setMethod("cpres_internal_getter", signature(x = "gseaResult"),function(x,...){
  y<-x@result
  if(dim(y)[1]==0){
    cli::cli_abort("no enrichment result contains")
  }
  if(sum(colnames(y)%in%c("ID","Description","NES","pvalue","p.adjust","core_enrichment"))!=6){
    cli::cli_abort("At list contains needed columns")
  }
  y<-y |> dplyr::select(ID,Description,NES,pvalue,p.adjust,core_enrichment)
  return(y)
})
setMethod("cpres_internal_getter", signature(x = "data.frame"),function(x,...){
  if(sum(colnames(x)=="NES")==0){
    if(sum(colnames(y)%in%c("ID","Description","GeneRatio","pvalue","p.adjust","geneID","Count"))!=7){
      cli::cli_abort("At list contains needed columns")
    }
    if(dim(y)[1]==0){
      cli::cli_abort("no enrichment result contains")
    }
    y<-y |> dplyr::select(ID,Description,GeneRatio,pvalue,p.adjust,geneID,Count)
    return(y)
  }
  else{
    if(sum(colnames(y)%in%c("ID","Description","NES","pvalue","p.adjust","core_enrichment"))!=6){
      cli::cli_abort("At list contains needed columns")
    }
    if(dim(y)[1]==0){
      cli::cli_abort("no enrichment result contains")
    }
    y<-y |> dplyr::select(ID,Description,NES,pvalue,p.adjust,core_enrichment)
    return(y)
  }

})

#' 2-Group Comparison of enrichment results and further clustering and visualizing
#'
#' @description
#' See `?EnrichGT()`
#'
#' @param obj.test the enriched object from tested group. WARNING: `obj.test` and `obj.ctrl` should come from same database (e.g. GO Biological Process(GOBP)).
#' @param obj.ctrl the enriched object from control group. WARNING: `obj.test` and `obj.ctrl` should come from same database (e.g. GO Biological Process(GOBP)).
#' @param name.test optional, the name of the testing group. If is `NULL`, the object name of `obj.test` will be used.
#' @param name.ctrl optional, the name of the control group. If is `NULL`, the object name of `obj.ctrl` will be used.
#' @param ClusterNum how many cluster will be clustered
#' @param P.adj p.adjust cut-off. To avoid slow visualization, you can make stricter p-cut off.
#' @param force ignore all auto-self-checks, which is useful
#' @param nTop keep n top items according to p-adj in each cluster.
#' @param method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param ... Others options.
#'
#' @details
#' Execute `obj.test` VS `obj.ctrl` tests, showing pathway overlaps (or differences) and meta-gene modules of test group and control group.
#'
#' Supports ORA and GSEA results (enriched object or data.frame).
#'
#' !WARNING!: `obj.test` and `obj.ctrl` should come from same database (e.g. GO Biological Process(GOBP)).
#' @return
#' `List` containing multiple `EnrichGT_obj` objects.
#'
#' The `List` contains objects with overlapped enriched terms, unique enrich terms.
#' @export
#'
#' @author Zhiming Ye
egt_compare_groups<-function(obj.test,obj.ctrl,name.test=NULL,name.ctrl=NULL,ClusterNum=15,P.adj=0.05,force=F,nTop=10,method="ward.D2",...){
  testdf<-cpres_internal_getter(obj.test)
  ctrldf<-cpres_internal_getter(obj.ctrl)
  overlapped <- testdf$ID[testdf$ID%in%ctrldf$ID]
  testOnly <- testdf$ID[!testdf$ID%in%ctrldf$ID]
  ctrlOnly <- ctrldf$ID[!ctrldf$ID%in%testdf$ID]
  ctrlOnlyDF<-ctrldf[ctrldf$ID%in%ctrlOnly,]
  testOnlyDF<-testdf[testdf$ID%in%testOnly,]
  ctrlOverlapDF<-ctrldf[!ctrldf$ID%in%ctrlOnly,]
  testOverlapDF<-testdf[!testdf$ID%in%testOnly,]
  Overlap_Control<-tryCatch({
    universalGT(ctrlOverlapDF,objname2="Overlap_Control",ClusterNum=ClusterNum,P.adj=P.adj,force=force,nTop=nTop,method=method,...)},error=function(e){
      message_egt("Error in parsing Overlap_Control, EnrichGT will return the raw data frame.",Type=1)
      return(ctrlOverlapDF)
    }
  )
  Overlap_Test<-tryCatch({
    universalGT(testOverlapDF,objname2="Overlap_Test",ClusterNum=ClusterNum,P.adj=P.adj,force=force,nTop=nTop,method=method,...)},error=function(e){
      message_egt("Error in parsing Overlap_Test, EnrichGT will return the raw data frame.",Type=1)
      return(testOverlapDF)
    }
  )
  Control_Only<-tryCatch({
    universalGT(ctrlOnlyDF,objname2="Control_Only",ClusterNum=ClusterNum,P.adj=P.adj,force=force,nTop=nTop,method=method,...)},error=function(e){
      message_egt("Error in parsing Control_Only, EnrichGT will return the raw data frame.",Type=1)
      return(ctrlOnlyDF)
    }
  )
  Test_Only<-tryCatch({
    universalGT(testOnlyDF,objname2="Test_Only",ClusterNum=ClusterNum,P.adj=P.adj,force=force,nTop=nTop,method=method,...)},error=function(e){
      message_egt("Error in parsing Test_Only, EnrichGT will return the raw data frame.",Type=1)
      return(testOnlyDF)
    }
  )
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
