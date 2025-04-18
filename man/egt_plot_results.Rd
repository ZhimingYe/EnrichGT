% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/figures.R
\name{egt_plot_results}
\alias{egt_plot_results}
\title{Visualize enrichment results using simple plot}
\usage{
egt_plot_results(
  x,
  ntop = NULL,
  showIDs = F,
  max_len_descript = 40,
  keepAll = F,
  ...,
  P.adj = NULL
)
}
\arguments{
\item{x}{a data frame form enriched result like `egt_enrichment_analysis()` or `egt_gsea_analysis()`, or an re-clustered `EnrichGT` object}

\item{ntop}{Show top N in each cluster. In default, for origin enriched result, showing top 15; for re-clustered object, showing top 5 in each cluster.}

\item{showIDs}{bool, show pathway IDs or not. Default is FALSE}

\item{max_len_descript}{the label format length, default as 40.}

\item{keepAll}{Do filtering to avoid overlap of same genes or not}

\item{...}{Other param}

\item{P.adj}{If pass an origin data.frame from original enriched result, you can specify the P-adjust value cut off. If is null, default is 0.05. When passing `EnrichGT_obj`, this filter is previously done by `egt_recluster_analysis`.}

\item{low.col}{the color for the lowest}

\item{hi.col}{the color for the highest}
}
\value{
a ggplot2 object
}
\description{
This plot is the most widely like `enrichplot::dotplot()`used method to visualize enriched terms. It shows the enrichment scores (e.g. p values) and gene ratio or NES as dot size and color / or bar height. Users can specify the number of terms using `ntop` or selected terms to color via the `low.col` and `hi.col`.
}
\author{
Zhiming Ye
}
