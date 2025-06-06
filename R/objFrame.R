setClass(
  "egt_llm",
  slots = list(
    pathways = "list",
    genes_and_title = "list"
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
    LLM_Annotation = "egt_llm"
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
  if ((x@LLM_Annotation@pathways |> length()) < 2) {
    summary_use_local(x, name)
  } else {
    summary_use_llm(x, name)
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
