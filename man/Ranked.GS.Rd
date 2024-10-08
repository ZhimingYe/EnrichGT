% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_helpers.R
\name{Ranked.GS}
\alias{Ranked.GS}
\title{Return ranked gene list which is use for "GSEA" and "cnetplot" in package "clusterProfiler"}
\usage{
Ranked.GS(Gene, log2FC, FromType = "SYMBOL", OrgDB = org.Hs.eg.db)
}
\arguments{
\item{Gene}{A vector containing genes}

\item{log2FC}{A vector containg log2FC}

\item{FromType}{default is "SYMBOL"}

\item{OrgDB}{default as org.Hs.eg.db}
}
\value{
A ranked named numeric vector. Names of the numbers is the ENTREZID.
}
\description{
Return ranked gene list which is use for "GSEA" and "cnetplot" in package "clusterProfiler"
}
\author{
Zhiming Ye
}
