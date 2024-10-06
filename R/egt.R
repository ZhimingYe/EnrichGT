
setGeneric("doEnrichGT",function(x,ClusterNum,P.adj=0.05,...) standardGeneric("doEnrichGT"))
setMethod("doEnrichGT", signature(x = "enrichResult"),function(x,...){
  if(sum(grepl("^GO",names(x@geneSets)))>5){
    nsimp()
  }
  y<-.genGT(x=x@result,ClusterNum,P.adj=0.05,nTop=nTop,...)
  return(y)
})
setMethod("doEnrichGT", signature(x = "compareClusterResult"),function(x,...){
  x<-x@compareClusterResult
  y<-.cprres(x,nTop=nTop)
  if(sum(grepl("^GO",names(y$ID)))>5){
    nsimp()
  }
  return(y)
})
setMethod("doEnrichGT", signature(x = "gseaResult"),function(x,...){
  y<-.genGSEAGT(x=x@result,ClusterNum,P.adj=0.05,nTop=nTop,...)
  return(y)
})
setMethod("doEnrichGT", signature(x = "data.frame"),function(x,...){
  if("NES"%in%colnames(x)){
    y<-.genGSEAGT(x=x,ClusterNum,P.adj=0.05,nTop=nTop,...)
    return(y)
  }
  else if("Cluster"%in%colnames(x)){
    y<-.cprres(x,nTop=nTop)
    if(sum(grepl("^GO",names(y$ID)))>5){
      nsimp()
    }
    return(y)
  }
  else{
    if(sum(grepl("^GO",names(y$ID)))>5){
      nsimp()
    }
    y<-.genGT(x=x,ClusterNum,P.adj=0.05,nTop=nTop,...)
    return(y)
  }
})

.genGT<-function(x,ClusterNum,P.adj=0.05,force=F,objname,nTop,...){
  InnerDF<-x
  .checkRowNames(x,"ORA")
  ClusterNum0<-ClusterNum
  ClusterNum<-.genClusterNum(x=x,ClusterNum = ClusterNum0,force = force) |> round()
  InnerDF<-InnerDF |> dplyr::left_join(.enrichpws(InnerDF$ID,InnerDF$geneID,ClusterNum)) # Merge according to "ID"
  InnerDF<-InnerDF |> dplyr::filter(pvalue<0.05,Count>=5,p.adjust<P.adj) |> dplyr::select(ID,Description,GeneRatio,`p.adjust`,geneID,Cluster,Count) # Need Fix
  .checkNrows(InnerDF,force = force)
  obj<-InnerDF |>
    dplyr::mutate(PCT=sapply(InnerDF$GeneRatio,function(x)eval(parse(text = x)))*100) |>
    dplyr::mutate(Padj = signif(p.adjust, 2),PCT=signif(PCT, 2)) |>
    dplyr::select(Description,ID,Count,Cluster,PCT,Padj,geneID) |>
    dplyr::mutate(geneID=gsub("/",", ",geneID)) # |>
  obj0 <-obj |> gt_ora()
  obj2 <-obj
  obj3 <-obj2 |> genMetaGM(type="ORA")
  obj3_1 <- obj3[[1]]
  obj3_2 <- obj3[[2]]
  return(obj)
}

.genGSEAGT<-function(x,ClusterNum,P.adj=0.05,force=F,objname,nTop,...){
  InnerDF<-x
  .checkRowNames(x,"GSEA")
  ClusterNum0<-ClusterNum
  ClusterNum<-.genClusterNum(x=x,ClusterNum = ClusterNum0,force = force) |> round()
  InnerDF<-InnerDF |> dplyr::left_join(.enrichpws(InnerDF$ID,InnerDF$core_enrichment,ClusterNum)) # Merge according to "ID"
  InnerDF<-InnerDF |> dplyr::filter(pvalue<0.05,abs(NES)>=0.9,p.adjust<P.adj) |> dplyr::select(ID,Description,NES,`p.adjust`,Cluster,core_enrichment) # Need Fix
  .checkNrows(InnerDF,force = force)
  obj<-InnerDF |>
    dplyr::mutate(absNES=abs(NES)) |>
    dplyr::mutate(Reg=ifelse(NES>0,"red","forestgreen")) |>
    dplyr::mutate(Padj = signif(p.adjust, 2),absNES=signif(absNES, 4)) |>
    dplyr::select(Description,ID,Reg,absNES,Cluster,Padj,core_enrichment) |>
    dplyr::mutate(core_enrichment=gsub("/",", ",core_enrichment))#|> need fix colors
  obj <- obj |>
    dplyr::group_by(Cluster) |>
    dplyr::arrange(Padj) |>
    dplyr::slice_head(n = nTop) |>
    dplyr::ungroup()
  obj0 <-obj |> gt_gsea()
  obj2 <- obj |> dplyr::mutate(Reg=ifelse(Reg=="red","UpReg","DownReg"))
  obj3 <-obj2 |> genMetaGM(type="GSEA")
  obj3_1 <- obj3[[1]]
  obj3_2 <- obj3[[2]]
  return(obj)
}
# attachment::att_amend_desc()
