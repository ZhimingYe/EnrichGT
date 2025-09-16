retry_function <- function(FUN, ntry = 5, delay = 3, ...) {
  for (attempt in 1:ntry) {
    result <- tryCatch(
      {
        FUN(...)
      },
      error = function(e) {
        if (attempt < ntry) {
          cli::cli_alert_info(sprintf(
            "Trying %d/%d Failed, Wait for %.1f secs for retry...",
            attempt,
            ntry,
            delay
          ))
          Sys.sleep(delay)
        }
        e
      }
    )
    if (!inherits(result, "error")) {
      return(result)
    }
  }
  cli::cli_abort(sprintf("Error after %d tries: %s", ntry, result$message))
}

check_is_ellmer_obj <- function(x) {
  if (sum(class(x) %in% c("Chat", "R6")) >= 2) T else wrong_llm()
}

chatter <- function(x, Prompt) {
  if (check_is_ellmer_obj(x)) {
    res <- retry_function(x$chat, ntry = 5, delay = 3, Prompt, echo = "none")
  } else {
    wrong_llm()
  }
  res
}
wrong_llm <- function() {
  cli::cli_alert_info(
    "The LLM support of `EnrichGT` is based on package `ellmer`\n For Example: "
  )
  cli::cli_code(
    "library(ellmer)\nchat <- chat_deepseek(api_key = YOUR_API, model = \"deepseek-chat\", system_prompt = \"\")\n# Please check https://ellmer.tidyverse.org/index.html for details"
  )
  cli::cli_abort(
    "LLM model registration Failed"
  )
}

summarize_clusters <- function(x, chat, prompt_type = "English") {
  promptList <- readRDS(system.file(
    "extdata",
    "egtllm.data",
    package = "EnrichGT"
  ))
  if (!inherits(x, "EnrichGT_obj")) {
    cli::cli_abort("Input must be an EnrichGT_obj object")
  }

  prompt <- switch(
    prompt_type,
    "English" = promptList$getModuleSummary_English,
    "Chinese" = promptList$getModuleSummary_Chinese,
    cli::cli_abort("Invalid prompt_type, must be 'English' or 'Chinese'")
  )

  cluster_names <- names(x@pathway_clusters)

  results <- lapply(
    cli::cli_progress_along(cluster_names, name = "Summarizing clusters"),
    function(cluster) {
      tryCatch(
        {
          cluster_prompt <- gsub(
            "<PATHWAYLIST>",
            paste(x@pathway_clusters[[cluster]], collapse = ", "),
            prompt
          )
          chatter(chat, cluster_prompt)
        },
        error = function(e) {
          cli::cli_alert_warning(
            glue::glue("Failed to process cluster {cluster}. ")
          )
          NULL
        }
      )
    }
  )
  return(list(results = results, cluster_names = cluster_names))
}


summarize_genes <- function(x, y, chat, prompt_type = "English") {
  promptList <- readRDS(system.file(
    "extdata",
    "egtllm.data",
    package = "EnrichGT"
  ))
  clustersName <- y[["cluster_names"]]
  y <- y[["results"]]
  have_exclude <- F
  excludeIndex <- c()
  if (!inherits(x, "EnrichGT_obj")) {
    cli::cli_abort("Input must be an EnrichGT_obj object")
  }
  if (sum(sapply(y, is.null)) > 0) {
    have_exclude <- T
    excludeIndex <- which(sapply(y, is.null))
    y <- y[-excludeIndex]
    clustersName <- clustersName[-excludeIndex]
    excludeIndex_warning <- paste(excludeIndex, collapse = ", ")
    cli::cli_alert_danger(glue::glue(
      "In cluster summarize, {excludeIndex_warning} were faild to fetch LLM summary, in next steps they will be excluded. "
    ))
  }
  prompt <- switch(
    prompt_type,
    "English" = promptList$getGeneSummary_English,
    "Chinese" = promptList$getGeneSummary_Chinese,
    cli::cli_abort("Invalid prompt_type, must be 'English' or 'Chinese'")
  )
  prompt_ <- switch(
    prompt_type,
    "English" = promptList$getTitle_English,
    "Chinese" = promptList$getTitle_Chinese,
    cli::cli_abort("Invalid prompt_type, must be 'English' or 'Chinese'")
  )
  need2summarize <- list()
  need2summarize_names <- list()
  index0 <- 1
  for (i in clustersName) {
    prompt2 <- gsub(
      "<GENELIST>",
      paste(x@gene_modules[[i]], collapse = ", "),
      prompt
    )
    prompt2 <- gsub("<BIOFUNS>", y[[which(clustersName == i)]], prompt2)
    prompt3 <- gsub(
      "<PRERES>",
      y[[which(clustersName == i)]],
      prompt_
    )
    need2summarize[index0] <- prompt2
    names(need2summarize[index0]) <- i
    need2summarize_names[index0] <- prompt3
    names(need2summarize_names[index0]) <- i
    index0 <- index0 + 1
  }
  results <- lapply(
    cli::cli_progress_along(clustersName, name = "Summarizing genes"),
    function(cluster) {
      tryCatch(
        {
          cluster_prompt <- need2summarize[[cluster]]
          chatter(chat, cluster_prompt)
        },
        error = function(e) {
          cli::cli_alert_warning(
            glue::glue("Failed to process cluster {cluster}. ")
          )
          NULL
        }
      )
    }
  )

  results2 <- lapply(
    cli::cli_progress_along(clustersName, name = "Summarizing Titles"),
    function(cluster) {
      tryCatch(
        {
          cluster_prompt <- need2summarize_names[[cluster]]
          chatter(chat, cluster_prompt)
        },
        error = function(e) {
          cli::cli_alert_warning(
            glue::glue("Failed to process cluster {cluster}. ")
          )
          NULL
        }
      )
    }
  )
  return(list(
    results = results,
    resultsTitle = results2,
    have_exclude = have_exclude,
    excludeIndex = excludeIndex,
    clustersName = clustersName
  ))
}

#' Summarize EnrichGT results using LLM
#'
#' This function uses a Large Language Model (LLM) to generate summaries for
#' pathway clusters and gene modules in an EnrichGT_obj object.
#'
#' @param x An EnrichGT_obj object created by \code{\link{egt_recluster_analysis}}.
#' @param chat An LLM chat object created by the \code{ellmer} package.
#' @param lang Language pass to LLM. Can be \code{English} or \code{Chinese}.
#'
#' @return Returns the input EnrichGT_obj object with added LLM annotations in
#' the \code{LLM_Annotation} slot. The annotations include:
#' \itemize{
#'   \item \code{pathways}: Summaries of pathway clusters
#'   \item \code{genes_and_title}: Summaries of gene modules and their titles
#' }
#'
#' @examples
#' \dontrun{
#' # Create LLM chat object
#' chat <- chat_deepseek(api_key = YOUR_API_KEY, model = "deepseek-chat", system_prompt = "")
#'
#' # Run enrichment analysis and get EnrichGT_obj
#' re_enrichment_results <- egt_recluster_analysis(...)
#'
#' # Get LLM summaries
#' re_enrichment_results <- egt_llm_summary(re_enrichment_results, chat)
#' }
#'
#' @references
#' For more information about creating chat objects, see the
#' \href{https://ellmer.tidyverse.org/index.html}{ellmer package documentation}.
#'
#' @note
#' It is recommended not to add system prompts when creating the chat object.
#' The function provides its own carefully crafted prompts for biological analysis.
#'
#' @seealso \code{\link{egt_recluster_analysis}} to create the input object.
#' @export
egt_llm_summary <- function(x, chat, lang = "English", model_name = NULL) {
  if (sum(lang %in% c("English", "Chinese")) != 1) {
    cli::cli_abort("Invalid prompt_type, must be 'English' or 'Chinese'")
  }
  if (class(x) != "EnrichGT_obj")
    cli::cli_abort("Please run `egt_recluster_analysis()` before summarizing. ")

  # Detect model name if not provided
  if (is.null(model_name)) {
    model_name <- tryCatch({
      if (exists("model", envir = chat)) {
        chat$model
      } else {
        "Unknown_Model"
      }
    }, error = function(e) "Unknown_Model")
  }

  a1 <- summarize_clusters(x, chat, lang)
  a2 <- summarize_genes(x, a1, chat, lang)
  obj_llm <- new("egt_llm")
  obj_llm@pathways <- a1
  obj_llm@genes_and_title <- a2
  obj_llm@llm_model_info <- paste0(model_name, "_", lang, "_", Sys.time())
  x@LLM_Annotation <- obj_llm
  return(x)
}

#' Compare Multiple LLM Summaries
#'
#' Generate summaries using multiple LLM models and create a comparison object.
#'
#' @param x An EnrichGT_obj object from egt_recluster_analysis()
#' @param chat_list A named list of LLM chat objects
#' @param lang Language for summaries ("English" or "Chinese")
#' @param comparison_prompt Custom prompt for generating comparison summary
#'
#' @return EnrichGT_obj with LLM_Comparison slot filled
#' @export
egt_llm_multi_summary <- function(x, chat_list, lang = "English", comparison_prompt = NULL) {
  if (class(x) != "EnrichGT_obj")
    cli::cli_abort("Please run `egt_recluster_analysis()` before summarizing.")

  if (!is.list(chat_list) || is.null(names(chat_list)))
    cli::cli_abort("chat_list must be a named list of LLM chat objects.")

  cli::cli_alert_info("Starting multi-LLM comparison...")

  # Generate summaries for each LLM
  llm_results <- list()
  model_names <- names(chat_list)

  for (i in seq_along(chat_list)) {
    model_name <- model_names[i]
    chat <- chat_list[[i]]

    cli::cli_alert_info(paste0("Generating summary with ", model_name, "..."))

    tryCatch({
      temp_result <- egt_llm_summary(x, chat, lang, model_name)
      llm_results[[model_name]] <- temp_result@LLM_Annotation
    }, error = function(e) {
      cli::cli_alert_warning(paste0("Failed to generate summary with ", model_name, ": ", e$message))
      llm_results[[model_name]] <- NULL
    })
  }

  # Remove failed results
  llm_results <- llm_results[!sapply(llm_results, is.null)]
  model_names <- names(llm_results)

  if (length(llm_results) == 0) {
    cli::cli_abort("No successful LLM summaries were generated.")
  }

  # Create comparison object
  comp_obj <- new("egt_llm_comparison")
  comp_obj@llm_results <- llm_results
  comp_obj@model_names <- model_names

  # For single LLM, store the result directly without comparison summary
  if (length(llm_results) == 1) {
    comp_obj@comparison_summary <- list()
    cli::cli_alert_success(paste0("Single LLM summary completed with ", model_names[1], "."))
  } else {
    # For multiple LLMs, generate comparison summary
    comparison_summary <- generate_comparison_summary(llm_results, model_names, lang, comparison_prompt)
    comp_obj@comparison_summary <- comparison_summary
    cli::cli_alert_success(paste0("Multi-LLM comparison completed with ", length(model_names), " models."))
  }

  x@LLM_Comparison <- comp_obj
  return(x)
}

generate_comparison_summary <- function(llm_results, model_names, lang, comparison_prompt = NULL) {
  cluster_names <- llm_results[[1]]@pathways$cluster_names
  comparison_summary <- list()

  for (i in seq_along(cluster_names)) {
    cluster <- cluster_names[i]
    cluster_comparison <- list()

    # Extract results for each model for this cluster
    for (model in model_names) {
      # Check if cluster exists in this model's results
      model_cluster_names <- llm_results[[model]]@genes_and_title$clustersName
      model_idx <- which(model_cluster_names == cluster)

      if (length(model_idx) > 0) {
        cluster_comparison[[model]] <- list(
          pathway_summary = llm_results[[model]]@pathways$results[[i]],
          gene_summary = llm_results[[model]]@genes_and_title$results[[model_idx]],
          title = llm_results[[model]]@genes_and_title$resultsTitle[[model_idx]]
        )
      } else {
        # If cluster doesn't exist in this model, use placeholder
        cluster_comparison[[model]] <- list(
          pathway_summary = "No pathway summary available for this cluster",
          gene_summary = "No gene summary available for this cluster",
          title = "No title available for this cluster"
        )
      }
    }

    # Store only model results without consensus analysis
    comparison_summary[[cluster]] <- list(
      model_results = cluster_comparison
    )
  }

  return(comparison_summary)
}
