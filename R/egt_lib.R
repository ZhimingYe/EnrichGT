genMetaGM<-function(x,type){
  if(type=="ORA"){
    x <- x |> dplyr::rename(genelists=geneID)
  }
  else{
    x <- x |> dplyr::rename(genelists=core_enrichment)
  }
  x$IDDs<-paste0(x$ID,"|",x$Description)
  x2 <- split(x,x$Cluster)
  x3 <- lapply(x2,function(y){
    q <- strsplit(y$genelists,", ") |> unlist() |> table() |> names()
    return(q)
  })
  x4 <- lapply(x2,function(y){
    q <- y$IDDs |> table() |> names()
    return(q)
  })
  return(list(`gene_modules`=x3,`pathway_clusters`=x4))
}

.cprres<-function(x,...){
  if("cluster"%in%colnames(x)){
    colnames(x)[colnames(x)=="cluster"]<-"zzz"
  }
  x<-x |> dplyr::rename(cluster=Cluster)
  if(sum(table(x$cluster)>10)>1){
    x<-x |> dplyr::filter(cluster%in%names(table(x$cluster)[table(x$cluster)>10]))
  }
  else{
    cli::cli_abort("Error.")
  }
  x<-split(x,x$cluster)
  tryCatch({y<-lapply(x, function(x2)try({.genGT(x2,...)}))},error=function(e){
    cli::cli_abort("[EnrichGT]Error: might be too few columns. ")
  })
  return(y)
}
is_numeric_string <- function(x) {
  grepl("^-?\\d+(\\.\\d+)?$", x)
}
.enrichpws<-function(ID,geneID,k,method,sep="/"){

  require(text2vec)
  tokens_list <- strsplit(geneID,sep)
  if(sum(is_numeric_string(unlist(tokens_list)))>length(unlist(tokens_list))*0.5){
    message_egt("Please run setReadable first!")
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
  hc <- hclust(distance_matrix, method = method)
  plot(hc)
  clusters <- cutree(hc, k = k)
  clusters<-clusters |> as.data.frame() |> tibble::rownames_to_column(var="ID")
  colnames(clusters)[2]<-"Cluster"
  clusters$Cluster<-paste0("Cluster_",clusters$Cluster)
  return(list(`clusters`=clusters,`hc`=hc,`dtm`=dtm))
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
      message_egt("Too many clusters! Try with max as 50...\nuse force=T to forbid the self-check")
      ClusterNum<-60
    }
    if(ClusterNum>dim(x)[1]/10 & dim(x)[1]>=50){
      message_egt("Too many clusters! Try with max as ncol/10...\nuse force=T to forbid the self-check")
      ClusterNum<-dim(x)[1]/11
    }
  }
  else{
  }
  return(ClusterNum)
}
.checkNrows<-function(x,force){
  if(nrow(x)>10000&!force){
    cli::cli_abort("[EnrichGT]Too many rows!(>10000), please subset! use force=T to forbid the self-check")
  }
  if(nrow(x)>750){
    message_egt("Too many rows! It might be slow...\nWorking, but please consider increase P.adj ...")
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
    cli::cli_abort(paste0("cols: ",paste(vec[!vec%in%colnames(x)],sep = ", "), " not found!\n"))
  }
}




#' @importFrom umap umap
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_classic
#' @importFrom ggrepel geom_text_repel
.egtUMAP <- function(x){
  mat<-x@document_term_matrix
  umap_result <- umap::umap(mat)
  umap_df <- data.frame(ID=rownames(umap_result[["layout"]]),
                        UMAP1 = umap_result$layout[, 1],
                        UMAP2 = umap_result$layout[, 2])
  udf<-x@enriched_result |> left_join(umap_df,by="ID")
  fig<-ggplot(udf, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    geom_point(size = 2) +
    geom_text_repel(aes(label = Description),
                    size = 3,
                    max.overlaps = 20,
                    box.padding = 0.3,
                    point.padding = 0.2) +
    labs(title = "Enrichment Results",
         x = "UMAP1", y = "UMAP2") +
    theme_classic()
  return(fig)
}
