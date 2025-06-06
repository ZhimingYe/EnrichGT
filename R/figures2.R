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


#' Fetch and visualize termwise relationships in enrichment analysis results
#'
#' This function analyzes the relationships between enriched terms and creates a
#' hierarchical clustering visualization.
#'
#' @param x A character vector of terms to analyze
#' @param database A data frame containing the database information
#' @param according_to Character string specifying which column to use for filtering, default is "Description"
#' @param ClusterNum Integer specifying the number of clusters, default is 4
#' @param method Character string specifying the hierarchical clustering method, default is "ward.D2"
#'
#' @return A list containing:
#'   \item{Figure}{A ggplot2 object showing the hierarchical clustering}
#'   \item{Cluster_DF}{A data frame with cluster assignments}
#'   \item{EnrichedHC}{The hierarchical clustering results}
#'
#' @importFrom dplyr filter left_join
#' @importFrom stringr str_wrap
#' @importFrom ggplot2 ggplot geom_segment geom_text coord_flip scale_y_continuous theme_void theme element_text margin expansion aes
#' @importFrom ggdendro dendro_data
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- egt_fetch_termwise_relationship(
#'   x = c("pathway1", "pathway2"),
#'   database = my_pathway_db,
#'   according_to = "Description"
#' )
#' result$Figure
#' }
#'
#' @export

egt_fetch_termwise_relationship <- function(
  x,
  database,
  according_to = "Description",
  ClusterNum = 4,
  method = "ward.D2"
) {
  if (according_to == "Description") {
    according_to <- "pathway"
  }
  stopifnot(is.data.frame(database))
  res <- prepare_database(database, db0_name = "pathway")
  ddd0 <- res$db0 |> dplyr::filter(!!sym(according_to) %in% x)
  ddd <- res$database |> dplyr::filter(term %in% ddd0$pathway)
  terms2 <- split(ddd$gene, ddd$term)
  kkk <- genclusterscore(terms2, names(terms2), ClusterNum, method)
  hc <- kkk$hc
  dd <- ggdendro::dendro_data(hc, type = "rectangle")
  seg_df <- dd$segments
  lab_df <- dd$labels
  cdf <- kkk$clusters
  colnames(cdf) <- c("label", "cluster")
  lab_df <- lab_df |> dplyr::left_join(cdf)
  lab_df$label <- stringr::str_wrap(lab_df$label, 40)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = seg_df,
      mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.5,
      colour = "grey40"
    ) +
    ggplot2::geom_text(
      data = lab_df,
      mapping = ggplot2::aes(
        x = x,
        y = -0.03 * max(seg_df$y),
        label = label,
        colour = cluster
      ),
      hjust = 1,
      size = 3,
      lineheight = 0.9
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.10))
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = .5, face = "bold"),
      plot.margin = ggplot2::margin(
        t = 10,
        r = 15,
        b = 10,
        l = 200,
        unit = "pt"
      )
    )
  print(p)
  invisible(list(Figure = p, Cluster_DF = kkk$clusters, EnrichedHC = kkk))
}
