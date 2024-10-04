
#' Parse enrichment results and further clustering and visualizing
#'
#' @param x an enrichment result from `clusterProfiler`, or a `data.frame` containing result from `clusterProfier`
#' @param ClusterNum how many cluster will be clustered
#' @param P.adj p.adjust cut-off. To avoid slow visualization, you can make stricter p-cut off.
#' @param ... Others options.
#' For an ORA result, c("ID","Description","GeneRatio","pvalue","p.adjust","geneID","Count") should be contained; For GSEA, c("ID","Description","NES","pvalue","p.adjust","core_enrichment") should be contain. For `compareClusterResult`, a `compareClusterResult` object or a data-frame with additional `Cluster` column should be contained, others similar to ORA result.
#' @return a `gt` great table object or a `List` containing multiple `gt` objects for `compareClusterResult`.
#' @export
#'
#' @author Zhiming Ye
EnrichGT<-function(x,ClusterNum=15,P.adj=0.05,...){
  res<-doEnrichGT(x,ClusterNum,P.adj,...)
  return(res)
}
