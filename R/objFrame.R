setClass(
  "egt_llm",
  slots = list(
    pathways = "list",
    genes_and_title = "list",
    llm_model_info = "character"
  )
)

setClass(
  "egt_llm_comparison",
  slots = list(
    llm_results = "list",
    model_names = "character",
    comparison_summary = "list"
  )
)
setOldClass("gt_tbl")
setOldClass("hclust")

setClass(
  "EnrichGT_obj",
  slots = list(
    enriched_result = "data.frame",
    gt_object = "gt_tbl",
    gt_object_noHTML = "gt_tbl",
    gene_modules = "list",
    pathway_clusters = "list",
    document_term_matrix = "dgCMatrix",
    clustering_tree = "hclust",
    raw_enriched_result = "data.frame",
    fused = "logical",
    param = "list",
    LLM_Annotation = "egt_llm",
    LLM_Comparison = "egt_llm_comparison"
  )
)

setMethod("show", "EnrichGT_obj", function(object) {
  if (isTRUE(getOption('knitr.in.progress'))) {
    cli::cli_alert_warning(
      "Please use `object@gt_object` to get the table when inside knitr(Quarto) environment, instead of simply calling this object. "
    )
    cli::cli_alert_success("`object@gt_object` is correct.")
    cli::cli_alert_danger("`object` only is wrong. ")
  } else {
    print(object@gt_object)
  }
})


setMethod("names", "EnrichGT_obj", function(x) {
  if (!is.data.frame(x@enriched_result)) return
  NULL
  if (nrow(x@enriched_result) > 3) {
    return(names(table(x@enriched_result$Cluster)))
  } else {
    return(NULL)
  }
})

setMethod("$", "EnrichGT_obj", function(x, name) {
  name <- getCLSNAME(name)

  # Check for multi-LLM comparison results first
  if (length(x@LLM_Comparison@llm_results) >= 2) {
    summary_use_llm_comparison(x, name)
  } else if (length(x@LLM_Comparison@llm_results) == 1) {
    # Single LLM comparison result - display it directly
    summary_use_single_llm_comparison(x, name)
  } else if ((x@LLM_Annotation@pathways |> length()) >= 2) {
    summary_use_llm(x, name)
  } else {
    summary_use_local(x, name)
  }
})

#' Generate HTML Summary for EnrichGT Object
#'
#' This function generates an HTML-formatted summary of enrichment results
#' that can be displayed in RStudio viewer or web browsers.
#'
#' @param x An EnrichGT_obj object
#' @param name The cluster name to summarize
#' @return HTML formatted summary that can be viewed with htmltools::html_print() or similar
#' @export
#' @examples
#' \dontrun{
#' # Generate HTML summary for cluster 1
#' html_summary <- egt_summary(obj, "Cluster_1")
#' htmltools::html_print(html_summary)
#'
#' # Or view in RStudio viewer
#' if (interactive()) {
#'   html_summary <- egt_summary(obj, "Cluster_1")
#'   print(html_summary)
#' }
#' }
#' @importFrom htmltools html_print
#' @importFrom htmltools HTML
egt_summary <- function(x, name) {
  if (!inherits(x, "EnrichGT_obj")) {
    stop("Input must be an EnrichGT_obj object")
  }

  name <- getCLSNAME(name)

  # Check for multi-LLM comparison results first
  if (length(x@LLM_Comparison@llm_results) >= 2) {
    r1 <- summary_use_llm_comparison_html(x, name)
  } else if (length(x@LLM_Comparison@llm_results) == 1) {
    # Single LLM comparison result - display it directly
    r1 <- summary_use_single_llm_comparison_html(x, name)
  } else if ((x@LLM_Annotation@pathways |> length()) >= 2) {
    r1 <- summary_use_llm_html(x, name)
  } else {
    r1 <- summary_use_local_html(x, name)
  }
  htmltools::html_print(r1)
}

getCLSNAME <- function(name){
  if (nchar(name) <= 5 & !(as.character(name) |> is_numeric_string())) {
    first_char <- substr(name, 1, 1)
    name <- gsub(first_char, "Cluster_", name)
  }
  if (as.character(name) |> is_numeric_string()) {
    name <- as.character(paste0("Cluster_", name))
  }
  name
}

new.egt <- function(x1, x2, x3, x4, x5, x6, x7, x8, x9) {
  flag0 <- F
  tryCatch(
    {
      objegt <- new("EnrichGT_obj")
      objegt@enriched_result <- x1
      objegt@gt_object <- x2
      objegt@gene_modules <- x3
      objegt@pathway_clusters <- x4
      objegt@document_term_matrix <- x5
      objegt@clustering_tree <- x6
      objegt@raw_enriched_result <- x7
      objegt@param <- x8
      objegt@fused <- F
      objegt@gt_object_noHTML <- x9
      flag0 <- T
    },
    error = function(e) {
      message_egt(
        "Failed to create EnrichGT object! Please re-check your input."
      )
      flag0 <- F
    }
  )
  if (flag0 == F) {
    objegt <- NULL
  }
  return(objegt)
}

#' Filter Enrichment Results by Description Pattern
#'
#' Infix operator to filter enrichment results by matching against Description field.
#' For EnrichGT_obj objects, re-runs clustering analysis after filtering.
#'
#' @param x Either an EnrichGT_obj object or data.frame containing enrichment results
#' @param y Regular expression pattern to match against Description field
#'
#' @return For EnrichGT_obj input: A new EnrichGT_obj with filtered and re-clustered results.
#' For data.frame input: A filtered data.frame.
#'
#' @details This operator helps refine enrichment results by removing terms matching
#' the given pattern from the Description field. When applied to EnrichGT_obj, it
#' preserves all original parameters and re-runs the clustering analysis on the
#' filtered results.
#'
#' @examples
#' \dontrun{
#' # Filter out "ribosome" related terms
#' filtered_results <- reenrichment_obj %-delete->% "ribosome"
#'
#' # Filter data.frame directly
#' filtered_df <- df %-delete->% "metabolism"
#' }
#'
#' @rdname MaskTerms
#' @export
`%-delete->%` <- function(x, y) {
  if (class(x) == "EnrichGT_obj") {
    paramList <- x@param
    paramList[["ClusterNum"]] -> ClusterNum
    paramList[["P.adj"]] -> P.adj
    paramList[["force"]] -> force
    paramList[["objname"]] -> objname
    paramList[["nTop"]] -> nTop
    paramList[["method"]] -> method
    x <- x@raw_enriched_result
    x <- x |> dplyr::filter(!grepl(y, x$Description))
    res2 <- egt_recluster_analysis(
      x = x,
      ClusterNum = ClusterNum,
      P.adj = P.adj,
      force = force,
      nTop = nTop,
      method = method
    )
    return(res2)
  } else if (is.data.frame(x)) {
    x <- x |> dplyr::filter(!grepl(y, x$Description))
    return(x)
  } else {
    cli::cli_abort("Please provide enriched results. ")
  }
}

#' @importFrom cli cli_ul cli_end
summary_use_local <- function(x, name) {
  resTable <- x@enriched_result |> dplyr::filter(Cluster == name)
  maxPrint <- ifelse(nrow(resTable) > 5, 5, nrow(resTable))
  DescriptPrint <- paste(resTable$Description[1:maxPrint], collapse = ", ")
  if (length(resTable$Description) == 0) {
    cli::cli_abort(
      "Please check your input. For example, you can use `Object$C1`, `Object$c1` or  `Object$Cluster_1` to reach the cluster 1\'s result summary"
    )
  }
  CandidateGenes <- paste(x@gene_modules[[name]], collapse = ", ")
  cli::cli_h1(glue::glue("Enrichment Result of {name} (Local Summary)"))
  ul <- cli::cli_ul()
  cli::cli_li(glue::glue("This cluster contains {DescriptPrint} ..."))
  cli::cli_li(glue::glue(
    "Candidate genes includes {CandidateGenes} (We will print all genes. Please stroll to top to read)"
  ))
  cli::cli_end(ul)
}

summary_use_local_html <- function(x, name) {
  resTable <- x@enriched_result |> dplyr::filter(Cluster == name)
  maxPrint <- ifelse(nrow(resTable) > 5, 5, nrow(resTable))
  DescriptPrint <- paste(resTable$Description[1:maxPrint], collapse = ", ")
  if (length(resTable$Description) == 0) {
    return(div(style = "color: #dc3545; margin: 15px;",
               "Please check your input. For example, you can use Object$C1, Object$c1 or Object$Cluster_1 to reach the cluster 1's result summary"))
  }
  CandidateGenes <- paste(x@gene_modules[[name]], collapse = ", ")

  div(style = "font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; margin: 20px;",
    h2(style = "color: #495057; border-bottom: 2px solid #007bff; padding-bottom: 10px; margin-bottom: 20px;",
       htmltools::HTML(glue::glue("Enrichment Result of {name} (Local Summary)"))),
    div(style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 15px;",
      p(style = "color: #495057; line-height: 1.6; margin-bottom: 10px;",
        htmltools::HTML(glue::glue("<strong>This cluster contains:</strong> {DescriptPrint} ..."))),
      p(style = "color: #495057; line-height: 1.6;",
        htmltools::HTML(glue::glue("<strong>Candidate genes include:</strong> {CandidateGenes}")))
    ),
    div(style = "color: #6c757d; font-size: 0.9em; font-style: italic;",
      "We can't find LLM annotations. Please using egt_llm_* functions for further analysis. ")
  )
}


summary_use_llm <- function(x, name) {
  tryCatch(
    {
      cli::cli_h1(glue::glue("Enrichment Result of {name} (LLM Summary)"))
      cli::cli_h2(x@LLM_Annotation@genes_and_title$resultsTitle[[which(
        x@LLM_Annotation@genes_and_title$clustersName == name
      )]])
      ul <- cli::cli_ul()
      cli::cli_li(x@LLM_Annotation@pathways$results[[which(
        x@LLM_Annotation@pathways$cluster_names == name
      )]])
      cli::cli_li(x@LLM_Annotation@genes_and_title$results[[which(
        x@LLM_Annotation@genes_and_title$clustersName == name
      )]])
      cli::cli_end(ul)
    },
    error = function(e) {
      cli::cli_alert_danger(
        "You already get LLM summary but EnrichGT can't fetch final results. Please recheck your API or network. "
      )
      summary_use_local(x, name)
    }
  )
}

summary_use_llm_html <- function(x, name) {
  tryCatch({
    title_idx <- which(x@LLM_Annotation@genes_and_title$clustersName == name)
    pathway_idx <- which(x@LLM_Annotation@pathways$cluster_names == name)

    if (length(title_idx) == 0 || length(pathway_idx) == 0) {
      return(div(style = "color: #dc3545; margin: 15px;",
                 "You already get LLM summary but EnrichGT can't fetch final results. Please recheck your API or network."))
    }

    title <- x@LLM_Annotation@genes_and_title$resultsTitle[[title_idx]]
    pathway_summary <- x@LLM_Annotation@pathways$results[[pathway_idx]]
    gene_summary <- x@LLM_Annotation@genes_and_title$results[[title_idx]]

    div(style = "font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; margin: 20px;",
      h2(style = "color: #495057; border-bottom: 2px solid #28a745; padding-bottom: 10px; margin-bottom: 20px;",
         htmltools::HTML(glue::glue("Enrichment Result of {name} (LLM Summary)"))),
      div(style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;",
        h3(style = "color: #28a745; margin-bottom: 15px;", title),
        div(style = "margin-bottom: 15px;",
          h4(style = "color: #17a2b8; margin-bottom: 10px;", "Pathway Analysis:"),
          div(style = "color: #495057; line-height: 1.6; white-space: pre-wrap; background-color: white; padding: 15px; border-radius: 5px; border-left: 4px solid #17a2b8;",
              htmltools::HTML(gsub("\n", "<br>", pathway_summary)))
        ),
        div(style = "margin-bottom: 10px;",
          h4(style = "color: #6f42c1; margin-bottom: 10px;", "Gene Analysis:"),
          div(style = "color: #495057; line-height: 1.6; white-space: pre-wrap; background-color: white; padding: 15px; border-radius: 5px; border-left: 4px solid #6f42c1;",
              htmltools::HTML(gsub("\n", "<br>", gene_summary)))
        )
      )
    )
  }, error = function(e) {
    summary_use_local_html(x, name)
  })
}

is_numeric_string <- function(x) {
  grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", x)
}

summary_use_llm_comparison <- function(x, name) {
  tryCatch({
    if (length(x@LLM_Comparison@llm_results) == 0) {
      cli::cli_alert_warning("No LLM comparison results available. Use egt_llm_multi_summary() first.")
      return(summary_use_local(x, name))
    }

    cli::cli_h1(glue::glue("Multi-LLM Comparison Result of {name}"))

    # Get comparison data for this cluster
    comparison_data <- x@LLM_Comparison@comparison_summary[[name]]
    model_names <- x@LLM_Comparison@model_names

    if (is.null(comparison_data)) {
      cli::cli_alert_warning(paste0("No comparison data found for cluster ", name))
      return(summary_use_local(x, name))
    }

    # Display results from each model
    cli::cli_h2("Individual Model Results")
    for (model in model_names) {
      if (!is.null(comparison_data$model_results[[model]])) {
        cli::cli_h3(glue::glue("Model: {model}"))
        ul1 <- cli::cli_ul()
        cli::cli_li(glue::glue("Title: {comparison_data$model_results[[model]]$title}"))
        cli::cli_li(glue::glue("Pathway Analysis: {comparison_data$model_results[[model]]$pathway_summary}"))
        cli::cli_li(glue::glue("Gene Analysis: {comparison_data$model_results[[model]]$gene_summary}"))
        cli::cli_end(ul1)
      }
    }

  }, error = function(e) {
    cli::cli_alert_danger("Error displaying LLM comparison results. Falling back to local summary.")
    summary_use_local(x, name)
  })
}

summary_use_llm_comparison_html <- function(x, name) {
  tryCatch({
    if (length(x@LLM_Comparison@llm_results) == 0) {
      return(div(style = "color: #dc3545; margin: 15px;",
                 "No LLM comparison results available. Use egt_llm_multi_summary() first."))
    }

    # Get comparison data for this cluster
    comparison_data <- x@LLM_Comparison@comparison_summary[[name]]
    model_names <- x@LLM_Comparison@model_names

    if (is.null(comparison_data)) {
      return(div(style = "color: #dc3545; margin: 15px;",
                 paste0("No comparison data found for cluster ", name)))
    }

    # Create model cards
    model_cards <- lapply(model_names, function(model) {
      if (!is.null(comparison_data$model_results[[model]])) {
        div(style = "background-color: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); border: 1px solid #e9ecef; height: 100%; display: flex; flex-direction: column;",
          h4(style = "color: #e83e8c; margin-bottom: 15px; padding-bottom: 8px; border-bottom: 2px solid #e83e8c; font-size: 16px;", 
             htmltools::HTML(glue::glue("Model: {model}"))),
          div(style = "flex-grow: 1;",
            div(style = "margin-bottom: 12px;",
              p(style = "color: #495057; line-height: 1.5; margin-bottom: 8px; font-size: 14px;",
                htmltools::HTML(glue::glue("<strong>Title:</strong> {comparison_data$model_results[[model]]$title}")))
            ),
            div(style = "margin-bottom: 12px;",
              p(style = "color: #495057; line-height: 1.5; margin-bottom: 5px; font-size: 14px; font-weight: 600;", 
                "Pathway Analysis:"),
              p(style = "color: #495057; line-height: 1.5; font-size: 13px; background-color: #f8f9fa; padding: 8px; border-radius: 4px; margin: 0;",
                htmltools::HTML(glue::glue("{comparison_data$model_results[[model]]$pathway_summary}")))
            ),
            div(style = "margin-bottom: 0;",
              p(style = "color: #495057; line-height: 1.5; margin-bottom: 5px; font-size: 14px; font-weight: 600;", 
                "Gene Analysis:"),
              p(style = "color: #495057; line-height: 1.5; font-size: 13px; background-color: #f8f9fa; padding: 8px; border-radius: 4px; margin: 0;",
                htmltools::HTML(glue::glue("{comparison_data$model_results[[model]]$gene_summary}")))
            )
          )
        )
      }
    })

    div(style = "font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; margin: 20px;",
      h2(style = "color: #495057; border-bottom: 2px solid #fd7e14; padding-bottom: 10px; margin-bottom: 20px;",
         htmltools::HTML(glue::glue("Multi-LLM Comparison Result of {name}"))),

      # Individual Model Results
      div(style = "margin-bottom: 25px;",
        # h3(style = "color: #495057; margin-bottom: 15px;", "Individual Model Results"),
        div(style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px;",
          model_cards
        )
      )
    )
  }, error = function(e) {
    summary_use_local_html(x, name)
  })
}
summary_use_single_llm_comparison <- function(x, name) {
  tryCatch({
    if (length(x@LLM_Comparison@llm_results) == 0) {
      cli::cli_alert_warning("No LLM comparison results available.")
      return(summary_use_local(x, name))
    }

    # Get the single model name and results
    model_name <- x@LLM_Comparison@model_names[1]
    llm_result <- x@LLM_Comparison@llm_results[[1]]

    cli::cli_h1(glue::glue("Enrichment Result of {name} (LLM Summary - {model_name})"))

    # Find the cluster index
    cluster_idx <- which(llm_result@pathways$cluster_names == name)

    if (length(cluster_idx) == 0) {
      cli::cli_alert_warning(paste0("No LLM summary found for cluster ", name))
      return(summary_use_local(x, name))
    }

    # Display title
    title <- llm_result@genes_and_title$resultsTitle[[cluster_idx]]
    cli::cli_h2(title)

    # Display pathway analysis
    pathway_summary <- llm_result@pathways$results[[cluster_idx]]
    cli::cli_h3("Pathway Analysis")
    ul1 <- cli::cli_ul()
    cli::cli_li(pathway_summary)
    cli::cli_end(ul1)

    # Display gene analysis
    gene_summary <- llm_result@genes_and_title$results[[cluster_idx]]
    cli::cli_h3("Gene Analysis")
    ul2 <- cli::cli_ul()
    cli::cli_li(gene_summary)
    cli::cli_end(ul2)

    # Display model information
    cli::cli_h3("Model Information")
    ul3 <- cli::cli_ul()
    cli::cli_li(glue::glue("Model: {model_name}"))
    cli::cli_li(glue::glue("Analysis Time: {llm_result@llm_model_info}"))
    cli::cli_end(ul3)

  }, error = function(e) {
    cli::cli_alert_danger("Error displaying single LLM comparison results. Falling back to local summary.")
    summary_use_local(x, name)
  })
}

summary_use_single_llm_comparison_html <- function(x, name) {
  tryCatch({
    if (length(x@LLM_Comparison@llm_results) == 0) {
      return(div(style = "color: #dc3545; margin: 15px;",
                 "No LLM comparison results available."))
    }

    # Get the single model name and results
    model_name <- x@LLM_Comparison@model_names[1]
    llm_result <- x@LLM_Comparison@llm_results[[1]]

    # Find the cluster index
    cluster_idx <- which(llm_result@pathways$cluster_names == name)

    if (length(cluster_idx) == 0) {
      return(div(style = "color: #dc3545; margin: 15px;",
                 paste0("No LLM summary found for cluster ", name)))
    }

    # Get data
    title <- llm_result@genes_and_title$resultsTitle[[cluster_idx]]
    pathway_summary <- llm_result@pathways$results[[cluster_idx]]
    gene_summary <- llm_result@genes_and_title$results[[cluster_idx]]

    div(style = "font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; margin: 20px;",
      h2(style = "color: #495057; border-bottom: 2px solid #6f42c1; padding-bottom: 10px; margin-bottom: 20px;",
         htmltools::HTML(glue::glue("Enrichment Result of {name} (LLM Summary - {model_name})"))),

      div(style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;",
        h3(style = "color: #6f42c1; margin-bottom: 15px;", title),

        div(style = "background-color: white; padding: 15px; border-radius: 5px; margin-bottom: 15px; border-left: 4px solid #28a745;",
          h4(style = "color: #28a745; margin-bottom: 10px;", "Pathway Analysis"),
          p(style = "color: #495057; line-height: 1.6; white-space: pre-wrap;",
            htmltools::HTML(gsub("\n", "<br>", pathway_summary)))
        ),

        div(style = "background-color: white; padding: 15px; border-radius: 5px; margin-bottom: 15px; border-left: 4px solid #17a2b8;",
          h4(style = "color: #17a2b8; margin-bottom: 10px;", "Gene Analysis"),
          p(style = "color: #495057; line-height: 1.6; white-space: pre-wrap;",
            htmltools::HTML(gsub("\n", "<br>", gene_summary)))
        ),

        div(style = "background-color: white; padding: 15px; border-radius: 5px; border-left: 4px solid #6c757d;",
          h4(style = "color: #6c757d; margin-bottom: 10px;", "Model Information"),
          p(style = "color: #495057; line-height: 1.5; margin-bottom: 5px;",
            htmltools::HTML(glue::glue("<strong>Model:</strong> {model_name}"))),
          p(style = "color: #495057; line-height: 1.5;",
            htmltools::HTML(glue::glue("<strong>Analysis Time:</strong> {llm_result@llm_model_info}")))
        )
      )
    )
  }, error = function(e) {
    summary_use_local_html(x, name)
  })
}
