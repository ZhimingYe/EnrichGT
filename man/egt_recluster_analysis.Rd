% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{egt_recluster_analysis}
\alias{egt_recluster_analysis}
\title{Cluster and re-enrichment enrichment results}
\usage{
egt_recluster_analysis(
  x,
  ClusterNum = 10,
  P.adj = 0.05,
  force = F,
  nTop = 10,
  method = "ward.D2",
  ...
)
}
\arguments{
\item{x}{Enrichment result from `EnrichGT` or `clusterProfiler`. For multi-database results, provide a `list`.}

\item{ClusterNum}{Number of clusters to create (default: 10).}

\item{P.adj}{Adjusted p-value cutoff (default: 0.05). Stricter values improve performance.}

\item{force}{Logical to bypass validation checks (default: FALSE).}

\item{nTop}{Number of top terms to keep per cluster by p-value (default: 10).}

\item{method}{Hierarchical clustering method (default: "ward.D2"). One of:
"ward.D", "ward.D2", "single", "complete", "average" (UPGMA),
"mcquitty" (WPGMA), "median" (WPGMC), or "centroid" (UPGMC).}

\item{...}{Additional arguments passed to clustering functions.}
}
\value{
An `EnrichGT_obj` containing:
\describe{
  \item{enriched_result}{Filtered results data frame}
  \item{gt_object}{Formatted `gt_tbl` table object}
  \item{gene_modules}{List of gene modules per cluster}
  \item{pathway_clusters}{Pathway names by cluster}
  \item{clustering_tree}{`hclust` object for visualization}
  \item{raw_enriched_result}{Unfiltered results table}
}
}
\description{
Performs hierarchical clustering on enrichment results (ORA or GSEA) based on gene-term associations to reduce redundancy and improve biological interpretation. The function helps identify coherent groups of related terms while preserving important but less significant findings.
}
\details{
Input requirements by analysis type:

ORA results:
  Required columns: "ID", "Description", "GeneRatio", "pvalue",
  "p.adjust", "geneID", "Count"

GSEA results:
  Required columns: "ID", "Description", "NES", "pvalue",
  "p.adjust", "core_enrichment"

compareClusterResult:
  Either the compareClusterResult object or a data frame with:
  - All ORA columns listed above
  - Additional "Cluster" column

Multi-database:
  Provide as a named list of the above result types
}
\examples{
\dontrun{
# ORA example
res <- egt_recluster_analysis(ora_result, ClusterNum=8)
plot(res@clustering_tree)

# GSEA example
gsea_res <- egt_recluster_analysis(gsea_result, method="average")
gsea_res
}

}
\author{
Zhiming Ye
}
