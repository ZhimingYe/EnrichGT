% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{egt_generate_quarto_report}
\alias{egt_generate_quarto_report}
\title{Export Quarto Report}
\usage{
egt_generate_quarto_report(
  re_enrichment_results,
  output_path = "Report.qmd",
  type = "html"
)
}
\arguments{
\item{re_enrichment_results}{The `EnrichGT_obj`, AI summarized result is more recommanded.}

\item{output_path}{Path of the output qmd file (e.g., `test.qmd`)}

\item{type}{export pdf or html. In default output is HTML format.}
}
\value{
A quarto document
}
\description{
Export Quarto Report
}
