#' Cluster and re-enrichment enrichment results
#'
#' @description
#' Performs hierarchical clustering on enrichment results (ORA or GSEA) based on gene-term associations to reduce redundancy and improve biological interpretation. The function helps identify coherent groups of related terms while preserving important but less significant findings.
#'
#' @param x Enrichment result from `EnrichGT` or `clusterProfiler`. For multi-database results, provide a `list`.
#' @param ClusterNum Number of clusters to create (default: 10).
#' @param P.adj Adjusted p-value cutoff (default: 0.05). Stricter values improve performance.
#' @param force Logical to bypass validation checks (default: FALSE).
#' @param nTop Number of top terms to keep per cluster by p-value (default: 10).
#' @param method Hierarchical clustering method (default: "ward.D2"). One of:
#' "ward.D", "ward.D2", "single", "complete", "average" (UPGMA),
#' "mcquitty" (WPGMA), "median" (WPGMC), or "centroid" (UPGMC).
#' @param ... Additional arguments passed to clustering functions.
#'
#' @details
#' Input requirements by analysis type:
#'
#' ORA results:
#'   Required columns: "ID", "Description", "GeneRatio", "pvalue",
#'   "p.adjust", "geneID", "Count"
#'
#' GSEA results:
#'   Required columns: "ID", "Description", "NES", "pvalue",
#'   "p.adjust", "core_enrichment"
#'
#' compareClusterResult:
#'   Either the compareClusterResult object or a data frame with:
#'   - All ORA columns listed above
#'   - Additional "Cluster" column
#'
#' Multi-database:
#'   Provide as a named list of the above result types
#'
#' @return An `EnrichGT_obj` containing:
#' \describe{
#'   \item{enriched_result}{Filtered results data frame}
#'   \item{tinytable_obj}{Formatted `tinytable` table object}
#'   \item{gene_modules}{List of gene modules per cluster}
#'   \item{pathway_clusters}{Pathway names by cluster}
#'   \item{clustering_tree}{`hclust` object for visualization}
#'   \item{raw_enriched_result}{Unfiltered results table}
#' }
#'
#' @examples
#' \dontrun{
#' # ORA example
#' res <- egt_recluster_analysis(ora_result, ClusterNum=8)
#' plot(res@clustering_tree)
#'
#' # GSEA example
#' gsea_res <- egt_recluster_analysis(gsea_result, method="average")
#' gsea_res
#' }
#'
#' @importFrom dplyr group_by arrange slice_head ungroup filter
#' @importFrom cli cli_alert_info cli_alert_warning cli_abort
#' @export
#'
#' @author Zhiming Ye
egt_recluster_analysis <- function(
  x,
  ClusterNum = 10,
  P.adj = 0.05,
  force = F,
  nTop = 10,
  method = "ward.D2",
  ...
) {
  objname <- deparse(substitute(x))
  if (objname == ".") {
    objname <- "`Magrittr` pipe conveyed object"
  }
  res <- doEnrichGT(
    x,
    ClusterNum,
    P.adj,
    force,
    objname = objname,
    nTop = nTop,
    method,
    ...
  )
  cli::cli_alert_success("re-enrichment done.")
  cli::cli_alert_info(
    "You can adjust the param of egt_recluster_analysis() for better results. Please refer to the help page. "
  )
  return(res)
}

.onAttach <- function(libname, pkgname) {
  required_packages <- c(
    "dplyr",
    "glue",
    "proxy",
    "RColorBrewer",
    "scales",
    "text2vec",
    "tibble",
    "forcats",
    "ggplot2",
    "Matrix",
    "xfun"
  )
  missing_packages <- required_packages[
    !(required_packages %in% installed.packages()[, "Package"])
  ]

  if (length(missing_packages) > 0) {
    package_message <- paste(
      "The following packages are not installed but are required: ",
      paste(missing_packages, collapse = ", ")
    )
    package_message <- paste(
      package_message,
      "\nPlease install them using install.packages() or BiocManager::install() for Bioconductor packages."
    )
    packageStartupMessage(package_message)
  } else {
    suppressPackageStartupMessages({
      requireNamespace("dplyr")
      requireNamespace("tibble")
      requireNamespace("ggplot2")
      requireNamespace("Matrix")
      requireNamespace("cli")
    })
    cli::cli_h1("EnrichGT")
    cli::cli_alert_info("See help on https://zhimingye.github.io/EnrichGT/")
    cli::cli_alert("by Zhiming Ye")
  }
}
