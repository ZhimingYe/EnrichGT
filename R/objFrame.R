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

is_numeric_string <- function(x) {
  grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", x)
}

#' Display Multi-LLM Comparison Summary
#' @keywords internal
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
    
    # Display consensus analysis
    cli::cli_h2("Consensus Analysis")
    ul2 <- cli::cli_ul()
    cli::cli_li(glue::glue("Common Themes: {comparison_data$consensus$common_themes}"))
    cli::cli_li(glue::glue("Agreement Level: {comparison_data$consensus$agreement_level}"))
    cli::cli_end(ul2)
    
    # Display differences
    cli::cli_h2("Model Differences")
    ul3 <- cli::cli_ul()
    cli::cli_li(glue::glue("Model-specific Insights: {comparison_data$differences$model_specific}"))
    cli::cli_li(glue::glue("Conflicting Views: {comparison_data$differences$conflicting_views}"))
    cli::cli_end(ul3)
    
    # Display confidence scores
    cli::cli_h2("Confidence Assessment")
    ul4 <- cli::cli_ul()
    cli::cli_li(glue::glue("Overall Confidence: {comparison_data$confidence_scores$overall}"))
    for (model in model_names) {
      if (!is.null(comparison_data$confidence_scores$per_model[[model]])) {
        cli::cli_li(glue::glue("{model} Confidence: {comparison_data$confidence_scores$per_model[[model]]}"))
      }
    }
    cli::cli_end(ul4)
    
  }, error = function(e) {
    cli::cli_alert_danger("Error displaying LLM comparison results. Falling back to local summary.")
    summary_use_local(x, name)
  })
}

#' Display Single LLM Comparison Summary
#' @keywords internal
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
