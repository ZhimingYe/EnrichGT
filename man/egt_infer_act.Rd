% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infAct.R
\name{egt_infer_act}
\alias{egt_infer_act}
\title{Infering Pathway or Transcript Factors activity from EnrichGT meta-gene modules}
\usage{
egt_infer_act(x, DB = "collectri", species = "human")
}
\arguments{
\item{x}{an EnrichGT_obj object.}

\item{DB}{can be "progeny" (the Pathway activity database), or "collectri" (TF activity database)}

\item{species}{can be "human" or "mouse"}
}
\value{
an ORA result list
}
\description{
Only supports gene symbols.

[PROGENy](https://saezlab.github.io/progeny/) is a comprehensive resource containing a curated collection of pathways and their target genes, with weights for each interaction.

[CollecTRI](https://github.com/saezlab/CollecTRI) is a comprehensive resource containing a curated collection of TFs and their transcriptional targets compiled from 12 different resources. This collection provides an increased coverage of transcription factors and a superior performance in identifying perturbed TFs compared to our previous.

If when doing re-enrichment, you select a high number of clusters, that may cause low gene number in each meta-gene module, and then can't be infered sucessfully. So if result is empty, please increase the number of re-clustering when doing it.
}
\author{
Zhiming Ye, Saez-Rodriguez Lab (The decoupleR package, https://saezlab.github.io/decoupleR/)
}
