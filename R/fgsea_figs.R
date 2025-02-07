#' Generate GSEA Enrichment Plots
#'
#' This function creates graphical representations of Gene Set Enrichment Analysis (GSEA) results, 
#' including either a multi-panel GSEA table plot for multiple pathways or a single pathway 
#' enrichment plot. The visualization leverages the \code{fgsea} package's plotting functions.
#'
#' @usage 
#' 
#' egt_plot_gsea(resGSEA$Description[1],genes = genes_with_weights(genes = DEGexample$...1, weights = DEGexample$log2FoldChange),database = database_GO_BP(org.Hs.eg.db))
#' 
#' egt_plot_gsea(resGSEA[1:8,],genes = genes_with_weights(genes = DEGexample$...1, weights = DEGexample$log2FoldChange),database = database_GO_BP(org.Hs.eg.db))
#' 
#' @param x A GSEA result object. Can be either:
#' \itemize{
#'   \item A data frame containing GSEA results (requires columns: pvalue, p.adjust, Description)
#'   \item A character string specifying a single pathway name
#' }
#' @param genes A named numeric vector from \code{genes_with_weights()}. These should match the gene identifiers
#' used in the GSEA analysis.
#' @param database A database \code{data.frame}, You can obtain it from \code{database_xxx()} functions.
#' This should correspond to the database used in the original GSEA analysis.
#'
#' @return A ggplot object:
#' \itemize{
#'   \item When \code{x} is a data frame: Returns a multi-panel plot showing normalized enrichment
#'         scores (NES), p-values, and leading edge plots for top pathways
#'   \item When \code{x} is a pathway name: Returns an enrichment plot showing the running
#'         enrichment score for the specified pathway
#' }
#'
#' @author Zhiming Ye, warpped from \code{fgsea}
#' @importFrom fgsea plotGseaTable
#' @importFrom fgsea plotEnrichment
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @export
#' 
egt_plot_gsea <- function(x, genes, database){
  suppressPackageStartupMessages({
    requireNamespace("fgsea")
  })
  ref <- egt_gsea_analysis_internal(genes=genes,database=database,for_figures=T)
  if(is.data.frame(x)){
    x$pval <- x$pvalue
    x$padj <- x$p.adjust
    x$pathway <- x$Description
    figure<-fgsea::plotGseaTable(ref[[1]][x$pathway], ref[[2]], x, gseaParam=0.5)
  }else{
    figure<-fgsea::plotEnrichment(ref[[1]][[x]], ref[[2]]) # not sure why 2[[]]?
    figure <- figure+ ggplot2::theme_bw() + ggplot2::labs(title = x)
    figure$layers[[1]]$aes_params$colour <- "#4975ae"
    figure$layers[[3]]$aes_params$colour <- "grey"
    figure$layers[[4]]$aes_params$colour <- "grey" # some hack of colors. 
  }
  return(figure)
}