% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multicpr.R
\name{comparison_reactor_ora}
\alias{comparison_reactor_ora}
\title{ORA Comparison Reactor}
\description{
Typically created via \code{\link{egt_comparison_reactor}("ORA")}
}
\details{
This reactor is optimized for comparing ORA results across multiple groups,
with methods tailored for p-value based comparisons.
}
\examples{

## ------------------------------------------------
## Method `comparison_reactor_ora$make_plans`
## ------------------------------------------------

\dontrun{
ora_reactor$make_plans(group = c("group1", "group2"), use_value = "padj")
}
}
\seealso{
\code{\link{comparison_reactor_base}} for inherited methods
A specialized reactor for comparing Over-Representation Analysis (ORA) results.
Inherits from comparison_reactor_base and provides ORA-specific functionality.
}
\section{Super class}{
\code{\link[EnrichGT:comparison_reactor_base]{EnrichGT::comparison_reactor_base}} -> \code{comparison_reactor_ora}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-comparison_reactor_ora-new}{\code{comparison_reactor_ora$new()}}
\item \href{#method-comparison_reactor_ora-make_plans}{\code{comparison_reactor_ora$make_plans()}}
\item \href{#method-comparison_reactor_ora-clone}{\code{comparison_reactor_ora$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="EnrichGT" data-topic="comparison_reactor_base" data-id="append_enriched_result"><a href='../../EnrichGT/html/comparison_reactor_base.html#method-comparison_reactor_base-append_enriched_result'><code>EnrichGT::comparison_reactor_base$append_enriched_result()</code></a></span></li>
<li><span class="pkg-link" data-pkg="EnrichGT" data-topic="comparison_reactor_base" data-id="do_recluster"><a href='../../EnrichGT/html/comparison_reactor_base.html#method-comparison_reactor_base-do_recluster'><code>EnrichGT::comparison_reactor_base$do_recluster()</code></a></span></li>
<li><span class="pkg-link" data-pkg="EnrichGT" data-topic="comparison_reactor_base" data-id="fetch_biological_theme"><a href='../../EnrichGT/html/comparison_reactor_base.html#method-comparison_reactor_base-fetch_biological_theme'><code>EnrichGT::comparison_reactor_base$fetch_biological_theme()</code></a></span></li>
<li><span class="pkg-link" data-pkg="EnrichGT" data-topic="comparison_reactor_base" data-id="fetch_relationship"><a href='../../EnrichGT/html/comparison_reactor_base.html#method-comparison_reactor_base-fetch_relationship'><code>EnrichGT::comparison_reactor_base$fetch_relationship()</code></a></span></li>
<li><span class="pkg-link" data-pkg="EnrichGT" data-topic="comparison_reactor_base" data-id="find_relationship"><a href='../../EnrichGT/html/comparison_reactor_base.html#method-comparison_reactor_base-find_relationship'><code>EnrichGT::comparison_reactor_base$find_relationship()</code></a></span></li>
<li><span class="pkg-link" data-pkg="EnrichGT" data-topic="comparison_reactor_base" data-id="get_recluster_result"><a href='../../EnrichGT/html/comparison_reactor_base.html#method-comparison_reactor_base-get_recluster_result'><code>EnrichGT::comparison_reactor_base$get_recluster_result()</code></a></span></li>
<li><span class="pkg-link" data-pkg="EnrichGT" data-topic="comparison_reactor_base" data-id="get_splited_list"><a href='../../EnrichGT/html/comparison_reactor_base.html#method-comparison_reactor_base-get_splited_list'><code>EnrichGT::comparison_reactor_base$get_splited_list()</code></a></span></li>
<li><span class="pkg-link" data-pkg="EnrichGT" data-topic="comparison_reactor_base" data-id="prefilter_by_p_adj"><a href='../../EnrichGT/html/comparison_reactor_base.html#method-comparison_reactor_base-prefilter_by_p_adj'><code>EnrichGT::comparison_reactor_base$prefilter_by_p_adj()</code></a></span></li>
<li><span class="pkg-link" data-pkg="EnrichGT" data-topic="comparison_reactor_base" data-id="print"><a href='../../EnrichGT/html/comparison_reactor_base.html#method-comparison_reactor_base-print'><code>EnrichGT::comparison_reactor_base$print()</code></a></span></li>
<li><span class="pkg-link" data-pkg="EnrichGT" data-topic="comparison_reactor_base" data-id="split_by_cluster"><a href='../../EnrichGT/html/comparison_reactor_base.html#method-comparison_reactor_base-split_by_cluster'><code>EnrichGT::comparison_reactor_base$split_by_cluster()</code></a></span></li>
<li><span class="pkg-link" data-pkg="EnrichGT" data-topic="comparison_reactor_base" data-id="summarize"><a href='../../EnrichGT/html/comparison_reactor_base.html#method-comparison_reactor_base-summarize'><code>EnrichGT::comparison_reactor_base$summarize()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-comparison_reactor_ora-new"></a>}}
\if{latex}{\out{\hypertarget{method-comparison_reactor_ora-new}{}}}
\subsection{Method \code{new()}}{
Create a new ORA comparison reactor
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{comparison_reactor_ora$new()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-comparison_reactor_ora-make_plans"></a>}}
\if{latex}{\out{\hypertarget{method-comparison_reactor_ora-make_plans}{}}}
\subsection{Method \code{make_plans()}}{
Create comparison plans between specified ORA groups
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{comparison_reactor_ora$make_plans(group = "auto", use_value = "p")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{group}}{Character vector of group names to compare or "auto" for all groups}

\item{\code{use_value}}{Which value to use for comparison ("p" for p-value or "padj" for adjusted p-value)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The reactor object (invisible) for method chaining
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
ora_reactor$make_plans(group = c("group1", "group2"), use_value = "padj")
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-comparison_reactor_ora-clone"></a>}}
\if{latex}{\out{\hypertarget{method-comparison_reactor_ora-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{comparison_reactor_ora$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
