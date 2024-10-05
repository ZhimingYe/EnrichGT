
#' Parse enrichment results and further clustering and visualizing
#'
#' @param x an enrichment result from `clusterProfiler`, or a `data.frame` containing result from `clusterProfier`
#' @param ClusterNum how many cluster will be clustered
#' @param P.adj p.adjust cut-off. To avoid slow visualization, you can make stricter p-cut off.
#' @param ... Others options.
#'
#'
#' @details
#'  For an ORA result, c("ID","Description","GeneRatio","pvalue","p.adjust","geneID","Count") should be contained; For GSEA, c("ID","Description","NES","pvalue","p.adjust","core_enrichment") should be contain. For `compareClusterResult`, a `compareClusterResult` object or a data-frame with additional `Cluster` column should be contained, others similar to ORA result.
#' @return a `gt` great table object or a `List` containing multiple `gt` objects for `compareClusterResult`.
#' @export
#'
#' @author Zhiming Ye
EnrichGT<-function(x,ClusterNum=15,P.adj=0.05,force=F,...){
  objname<-deparse(substitute(x))
  if(objname=="."){
    objname<-"`Magrittr` pipe conveyed object"
  }
  res<-doEnrichGT(x,ClusterNum,P.adj,force,objname=objname,...)
  return(res)
}

.onAttach <- function(libname, pkgname) {
  required_packages <- c(
    "dplyr", "fontawesome", "glue", "gt", "proxy",
    "RColorBrewer", "rlang", "scales", "text2vec",
    "tibble", "clusterProfiler", "ReactomePA"
  )
  missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

  if (length(missing_packages) > 0) {
    package_message <- paste("The following packages are not installed but are required: ",
                             paste(missing_packages, collapse = ", "))
    package_message <- paste(package_message, "\nPlease install them using install.packages() or BiocManager::install() for Bioconductor packages.")
    packageStartupMessage(package_message)
  } else {
    packageStartupMessage("View your enrichment result by entring `EnrichGT(result)`\nby Zhiming Ye, https://github.com/ZhimingYe/EnrichGT")
  }
}
