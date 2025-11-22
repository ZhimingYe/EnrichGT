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
#'   \item{gt_object}{Formatted `gt_tbl` table object}
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
  suppressPackageStartupMessages({
    requireNamespace("dplyr")
    requireNamespace("tibble")
    requireNamespace("ggplot2")
    requireNamespace("Matrix")
    requireNamespace("cli")
  })
  message("EnrichGT, See help on https://zhimingye.github.io/EnrichGT/")
}


#' Export Quarto Report
#'
#' @param re_enrichment_results The `EnrichGT_obj`, AI summarized result is more recommanded.
#' @param output_path Path of the output qmd file (e.g., `test.qmd`)
#' @param type export pdf or html. In default output is HTML format.
#'
#' @returns A quarto document
#' @export
#'
egt_generate_quarto_report <- function(
  re_enrichment_results,
  output_path = "Report.qmd",
  type = "html"
) {
  rds_path <- file.path(
    dirname(output_path),
    paste0(output_path, "dependency.rds")
  )
  HTMLtitle <- c(
    '---',
    'title: "Re-enrich Enriched Results"',
    'format:',
    '  html:',
    '    toc: true',
    '    self-contain: true',
    '---',
    '',
    '# Enrichment Summary',
    '',
    '```{r}',
    '#| echo: false',
    '#| message: false',
    '#| warning: false',
    'library(EnrichGT)',
    paste0('re_enrichment_results <- readRDS("', basename(rds_path), '")'),
    '```',
    ''
  )
  PDFtitle <- c(
    '---',
    'title: "Re-enrich Enriched Results"',
    'format:',
    '  pdf:',
    '    toc: true',
    '    number-sections: true',
    '---',
    '',
    '# Enrichment Summary',
    '',
    '```{r}',
    '#| echo: false',
    '#| message: false',
    '#| warning: false',
    'library(EnrichGT)',
    paste0('re_enrichment_results <- readRDS("', basename(rds_path), '")'),
    '```',
    ''
  )
  HTMLGTrmd <- c(
    '# Enrichment Full results',
    '',
    '```{r, comment = ""}',
    '#| echo: false',
    're_enrichment_results@gt_object',
    '```'
  )
  PDFGTrmd <- c(
    '# Enrichment Full results',
    '',
    '```{r, comment = ""}',
    '#| echo: false',
    're_enrichment_results@gt_object_noHTML',
    '```'
  )

  if (!class(re_enrichment_results) == "EnrichGT_obj") {
    cli::cli_abort(
      "Please provide `EnrichGT_obj`. You should perform re-clustering"
    )
  }

  saveRDS(re_enrichment_results, file = rds_path)
  cluster_names <- names(re_enrichment_results)
  if (type == "pdf") {
    qmd_content <- PDFtitle
  } else {
    qmd_content <- HTMLtitle
  }
  qmd_content <- qmd_content
  for (cluster in cluster_names) {
    qmd_content <- c(
      qmd_content,
      paste0('## ', cluster),
      '',
      '```{r, comment = ""}',
      '#| echo: false',
      '#| output: asis',
      paste0(
        "a <- EnrichGT:::egt_summary(re_enrichment_results, '",
        cluster,
        "')"
      ),
      'a',
      '```',
      '',
      ''
    )
  }
  if (type == "pdf") {
    qmd_content <- c(
      qmd_content,
      PDFGTrmd
    )
  } else {
    qmd_content <- c(
      qmd_content,
      HTMLGTrmd
    )
  }

  writeLines(qmd_content, con = output_path)
  cli::cli_alert_success(paste("Report generated at:", output_path))
  cli::cli_alert_success(paste("RDS file saved at:", rds_path))
}
