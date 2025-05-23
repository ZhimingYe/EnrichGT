% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objFrame.R
\name{\%-delete->\%}
\alias{\%-delete->\%}
\title{Filter Enrichment Results by Description Pattern}
\usage{
x \%-delete->\% y
}
\arguments{
\item{x}{Either an EnrichGT_obj object or data.frame containing enrichment results}

\item{y}{Regular expression pattern to match against Description field}
}
\value{
For EnrichGT_obj input: A new EnrichGT_obj with filtered and re-clustered results.
For data.frame input: A filtered data.frame.
}
\description{
Infix operator to filter enrichment results by matching against Description field.
For EnrichGT_obj objects, re-runs clustering analysis after filtering.
}
\details{
This operator helps refine enrichment results by removing terms matching
the given pattern from the Description field. When applied to EnrichGT_obj, it
preserves all original parameters and re-runs the clustering analysis on the
filtered results.
}
\examples{
\dontrun{
# Filter out "ribosome" related terms
filtered_results <- reenrichment_obj \%-delete->\% "ribosome"

# Filter data.frame directly
filtered_df <- df \%-delete->\% "metabolism"
}

}
