#' Generate Biological Theme Wordcloud from Enrichment Results
#'
#' Creates a wordcloud visualization from enrichment results, either from a data.frame
#' or an EnrichGT_obj. For data.frames, filters by p.adjust < 0.05 by default.
#'
#' @param x Input object - either a data.frame with enrichment results or an EnrichGT_obj
#' @param cluster Cluster name/number (required if x is EnrichGT_obj)
#' @param skip_filtering Logical, if TRUE skips filtering of data.frame by p.adjust
#' @param ... Additional arguments passed to wordcloud_generator2
#' @return ggwordcloud object
#' @author Zhiming Ye
#' @export
egt_fetch_biological_theme <- function(
  x,
  cluster = NULL,
  skip_filtering = F,
  ...
) {
  if (is.data.frame(x)) {
    if (!skip_filtering) {
      x <- x |> dplyr::filter(p.adjust < 0.05)
    }
    wordcloud_generator2(x$Description, ...)
  } else if (class(x) == "EnrichGT_obj") {
    if (is.null(cluster)) {
      cli::cli_abort("Please provide cluster name or number")
    }
    name <- getCLSNAME(cluster)
    wordcloud_generator2(sapply(
      strsplit(x@pathway_clusters[[name]], "[|]"),
      function(q) q[2]
    ))
  } else NULL
}
