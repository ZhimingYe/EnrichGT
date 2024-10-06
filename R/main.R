
setGeneric("doEnrichGT",function(x,ClusterNum,P.adj=0.05,...) standardGeneric("doEnrichGT"))
setMethod("doEnrichGT", signature(x = "enrichResult"),function(x,...){
  if(sum(grepl("^GO",names(x@geneSets)))>5){
    message("\n=====[SUGGESTION]=====\nYou are passing an object from GO Enrichment.\nPlease ensure that `obj |> clusterProfiler::simplify()` is executed, to pre-simplify result,\nFor better enriched result.\n")
  }
  y<-.genGT(x=x@result,ClusterNum,P.adj=0.05,...)
  return(y)
})
setMethod("doEnrichGT", signature(x = "compareClusterResult"),function(x,...){
  x<-x@compareClusterResult
  y<-.cprres(x)
  if(sum(grepl("^GO",names(y$ID)))>5){
    message("\n=====[SUGGESTION]=====\nYou are passing an object from GO Enrichment.\nPlease ensure that `obj |> clusterProfiler::simplify()` is executed, to pre-simplify result,\nFor better enriched result.\n")
  }
  return(y)
})
setMethod("doEnrichGT", signature(x = "gseaResult"),function(x,...){
  y<-.genGSEAGT(x=x@result,ClusterNum,P.adj=0.05,...)
  return(y)
})
setMethod("doEnrichGT", signature(x = "data.frame"),function(x,...){
  if("NES"%in%colnames(x)){
    y<-.genGSEAGT(x=x,ClusterNum,P.adj=0.05,...)
    return(y)
  }
  else if("Cluster"%in%colnames(x)){
    y<-.cprres(x)
    if(sum(grepl("^GO",names(y$ID)))>5){
      message("\n=====[SUGGESTION]=====\nYou are passing an object from GO Enrichment.\nPlease ensure that `obj |> clusterProfiler::simplify()` is executed, to pre-simplify result,\nFor better enriched result.\n")
    }
    return(y)
  }
  else{
    if(sum(grepl("^GO",names(y$ID)))>5){
      message("\n=====[SUGGESTION]=====\nYou are passing an object from GO Enrichment.\nPlease ensure that `obj |> clusterProfiler::simplify()` is executed, to pre-simplify result,\nFor better enriched result.\n")
    }
    y<-.genGT(x=x,ClusterNum,P.adj=0.05,...)
    return(y)
  }
})
.cprres<-function(x,...){

  if("cluster"%in%colnames(x)){
    colnames(x)[colnames(x)=="cluster"]<-"zzz"
  }
  x<-x |> dplyr::rename(cluster=Cluster)
  if(sum(table(x$cluster)>10)>1){
    x<-x |> dplyr::filter(cluster%in%names(table(x$cluster)[table(x$cluster)>10]))
  }
  else{
    stop("Error.")
  }
  x<-split(x,x$cluster)
  tryCatch({y<-lapply(x, function(x2)try({.genGT(x2,ClusterNum,P.adj=0.05,...)}))},error=function(e){
    stop("[Message]Error: might be too few columns. ")
  })
  return(y)
}
is_numeric_string <- function(x) {
  grepl("^-?\\d+(\\.\\d+)?$", x)
}
.enrichpws<-function(ID,geneID,k,sep="/"){
  require(proxy)
  require(text2vec)
  tokens_list <- strsplit(geneID,sep)
  if(sum(is_numeric_string(unlist(tokens_list)))>length(unlist(tokens_list))*0.5){
    message("[Message]Please run setReadable first!")
  }
  tokens_list<-lapply(tokens_list,function(x){
    x[is_numeric_string(x)]<-paste0("Gene_",x)
  })

  names(tokens_list)<-ID
  tokens_iter <- text2vec::itoken(tokens_list, progressbar = FALSE)
  vocab <- text2vec::create_vocabulary(tokens_iter)
  vectorizer <- text2vec::vocab_vectorizer(vocab)
  dtm <- text2vec::create_dtm(tokens_iter, vectorizer)
  distance_matrix <- proxy::dist(dtm |> as.matrix(), method = "cosine")
  hc <- hclust(distance_matrix, method = "ward.D2")
  plot(hc)
  clusters <- cutree(hc, k = k)
  clusters<-clusters |> as.data.frame() |> tibble::rownames_to_column(var="ID")
  colnames(clusters)[2]<-"Cluster"
  clusters$Cluster<-paste0("Cluster_",clusters$Cluster)
  return(clusters)
}


.genClusterNum<-function(x,ClusterNum,force){
  if(!force){
    ClusterNum0<-ClusterNum
    ClusterNum<-dplyr::case_when(dim(x)[1]<10~1,
                                 dim(x)[1]<20~2,
                                 dim(x)[1]<30~2,
                                 dim(x)[1]<40~3,
                                 dim(x)[1]<50~3,
                                 dim(x)[1]>=50~ClusterNum0)
    if(ClusterNum>60){
      message("[Message]Too many clusters! Try with max as 50...\nuse force=T to forbid the self-check")
      ClusterNum<-60
    }
    if(ClusterNum>dim(x)[1]/10 & dim(x)[1]>=50){
      message("[Message]Too many clusters! Try with max as ncol/10...\nuse force=T to forbid the self-check")
      ClusterNum<-dim(x)[1]/11
    }
  }
  else{
  }
  return(ClusterNum)
}
.checkNrows<-function(x,force){
  if(nrow(x)>1000&!force){
    stop("[Message]Too many rows!(>1000), please subset! use force=T to forbid the self-check")
  }
  if(nrow(x)>750){
    message("[Message]Too many rows! It might be slow...\nWorking, but please consider increase P.adj ...")
  }
}
.checkRowNames<-function(x,Type){
  if(Type=="ORA"){
    judge<-.TypeChecker(x,c("ID","Description","GeneRatio","pvalue","p.adjust","geneID","Count"))
    return(judge)
  }
  else if(Type=="GSEA"){
    judge<-.TypeChecker(x,c("ID","Description","NES","pvalue","p.adjust","core_enrichment"))
    return(judge)
  }
}
.TypeChecker<-function(x,vec){
  judgeinner<-F
  if(sum(colnames(x)%in%vec)==length(vec)){
    judgeinner<-T
    return(judgeinner)
  }
  else{
    stop(paste0("cols: ",paste(vec[!vec%in%colnames(x)],sep = ", "), " not found!\n"))
  }
}

.genGT<-function(x,ClusterNum,P.adj=0.05,force=F,objname,...){
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
    dplyr::mutate(geneID=gsub("/",", ",geneID))|>
    gt::gt(groupname_col = "Cluster") |>
    gt_hulk_col_numeric2(PCT,pal = RColorBrewer::brewer.pal(8,"PiYG") |> rev()) |>
    gt_hulk_col_numeric2(Padj,pal = RColorBrewer::brewer.pal(8,"Spectral")) |>
    gt_hulk_col_numeric2(Count,pal = RColorBrewer::brewer.pal(8,"PuBuGn") |> rev()) |>
    gt_merge_stack2(col2 = ID, col1 = Description) |>
    gt::tab_style(
      style = cell_text(size = px(13)),
      locations = cells_body()
    ) |>
    gt::tab_header(title =paste0("Parse form: ",objname),subtitle = paste0("Split into ",ClusterNum, " Clusters. Generated by github@zhimingye/EnrichGT"))
  return(obj)
}

.genGSEAGT<-function(x,ClusterNum,P.adj=0.05,force=F,objname,...){
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
    dplyr::mutate(core_enrichment=gsub("/",", ",core_enrichment))|>
    gt::gt(groupname_col = "Cluster") |>
    gt_hulk_col_numeric2(Padj,pal = RColorBrewer::brewer.pal(8,"Spectral")) |>
    gt_hulk_col_numeric2(absNES,pal = RColorBrewer::brewer.pal(8,"PuBuGn") |> rev()) |>
    gt_merge_stack2(col2 = ID, col1 = Description) |>
    # gt::cols_add(dir = ifelse(Reg=="Up", "red", "forestgreen")) |>
    gt::cols_label(Reg = "") |>
    gt::text_case_match(
      "red" ~ fontawesome::fa("arrow-up"),
      "forestgreen" ~ fontawesome::fa("arrow-down")
    ) |>
    gt::tab_style(
      style = cell_text(color = from_column("Reg")),
      locations = cells_body(columns = Reg)
    ) |>
    gt::tab_style(
      style = cell_text(size = px(13)),
      locations = cells_body()
    ) |>
    gt::tab_header(title =paste0("Parse form: ",objname),subtitle = paste0("Split into ",ClusterNum, " Clusters. Generated by github@zhimingye/EnrichGT"))
  return(obj)
}
# attachment::att_amend_desc()
