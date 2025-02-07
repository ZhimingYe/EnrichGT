
#' Parse enrichment results and further clustering and visualizing
#'
#' @description
#' Cluster enrichment results based on hit genes for ORA (e.g, typical GO enrichment) or core enrichment from GSEA using term frequency analysis. This provides a clearer view of biological relevance by focusing on the genes that matter most.
#'
#' Gene enrichment analysis can often be misleading due to the redundancy within gene set databases and the limitations of most enrichment tools. Many tools, by default, only display a few top results and fail to filter out redundancy. This can result in both biological misinterpretation and valuable information being overlooked.
#'
#' For instance, high expression of certain immune genes can cause many immune-related gene sets to appear overrepresented. However, a closer look often reveals that these gene sets are derived from the same group of genes, which might represent only a small fraction. Less than 1/10 of the differentially expressed genes (DEGs). What about the other 9/10?  Do they hold no biological significance?
#'
#' The main purpose of developing this package is to provide a lightweight and practical solution to the problems mentioned above.
#'
#' @param x an enrichment result from `EnrichGT` or `clusterProfiler`. To perform fusing multi-database enrichment results, please give a `list` object.
#' @param ClusterNum how many cluster will be clustered
#' @param P.adj p.adjust cut-off. To avoid slow visualization, you can make stricter p-cut off.
#' @param force ignore all auto-self-checks, which is useful
#' @param nTop keep n top items according to p-adj in each cluster.
#' @param method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param ... Others options.
#'
#' @details
#'  For an ORA result, c("ID","Description","GeneRatio","pvalue","p.adjust","geneID","Count") should be contained;
#'
#'  For GSEA, c("ID","Description","NES","pvalue","p.adjust","core_enrichment") should be contain.
#'
#'  For `compareClusterResult`, a `compareClusterResult` object or a data-frame with additional `Cluster` column should be contained, others similar to ORA result.
#'
#'  To perform fusing multi-database enrichment results, please give a `list` object.
#'
#' @return an `EnrichGT_obj` object.
#'
#' slot `enriched_result` contains a data.frame with enriched results. `gt_object` contains `gt` object.
#'
#' you can use `obj@gt_object` to get it and use functions from `gt` like `gtsave`.
#'
#' `gene_modules` is a list containing meta-gene modules of each cluster.
#'
#' `pathway_clusters` contains pathways names in each cluster.
#'
#' `clustering_tree` contains the clustering tree object from `hclust()`, you can use other packages like `ggtree` for further visualization and analysis.
#'
#' `raw_enriched_result` contains raw table without selecting `nTop`.
#'
#' @export
#'
#' @author Zhiming Ye
egt_recluster_analysis<-function(x,ClusterNum=10,P.adj=0.05,force=F,nTop=10,method="ward.D2",...){
  objname<-deparse(substitute(x))
  if(objname=="."){
    objname<-"`Magrittr` pipe conveyed object"
  }
  res<-doEnrichGT(x,ClusterNum,P.adj,force,objname=objname,nTop=nTop,method,...)
  cli::cli_alert_success("re-enrichment done.")
  cli::cli_alert_info("You can adjust the param of egt_recluster_analysis() for better results. Please refer to the help page. ")
  return(res)
}

.onAttach <- function(libname, pkgname) {
  required_packages <- c(
    "dplyr", "fontawesome", "glue", "gt", "proxy",
    "RColorBrewer", "rlang", "scales", "text2vec",
    "tibble","forcats","ggplot2", "Matrix", "xfun"
  )
  missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

  if (length(missing_packages) > 0) {
    package_message <- paste("The following packages are not installed but are required: ",
                             paste(missing_packages, collapse = ", "))
    package_message <- paste(package_message, "\nPlease install them using install.packages() or BiocManager::install() for Bioconductor packages.")
    packageStartupMessage(package_message)
  } else {
    suppressPackageStartupMessages({
      requireNamespace("dplyr")
      requireNamespace("tibble")
      requireNamespace("ggplot2")
      requireNamespace("Matrix")
      requireNamespace("gt")
      requireNamespace("cli")
    })
    cli::cli_h1("EnrichGT Version 0.9")
    cli::cli_alert_info("See help on https://zhimingye.github.io/EnrichGT/")
    cli::cli_alert("by Zhiming Ye")
  }
}
