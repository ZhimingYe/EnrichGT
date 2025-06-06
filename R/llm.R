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
  prompt3 <- switch(
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
      prompt3
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
egt_llm_summary <- function(x, chat, lang = "English") {
  if (sum(lang %in% c("English", "Chinese")) != 1) {
    cli::cli_abort("Invalid prompt_type, must be 'English' or 'Chinese'")
  }
  if (class(x) != "EnrichGT_obj")
    cli::cli_abort("Please run `egt_recluster_analysis()` before summarizing. ")
  a1 <- summarize_clusters(x, chat, lang)
  a2 <- summarize_genes(x, a1, chat, lang)
  obj_llm <- new("egt_llm")
  obj_llm@pathways <- a1
  obj_llm@genes_and_title <- a2
  x@LLM_Annotation <- obj_llm
  return(x)
}
