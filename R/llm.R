promptList <- list(
  getModuleSummary_Chinese = "这是一些生物学富集分析后重新再次整理后得到的生物学主题列表。你是一名专业生物学家,帮我用一些话概括这些通路的主题是什么,反映了什么,有什么生物学解释的角度。请简洁有力,必要时提供依据。请您还需要为这个模块拟定一个代表性名字。请使用学术写作风格。请控制字数在250字左右,避免长篇大论。请输出纯文本,不包括markdown格式。参与这个模块的通路包括:<PATHWAYLIST>",
  getGeneSummary_Chinese = "这是一系列Gene,可能具有类似生物学功能。前面已经总结他们的生物学功能可能包括\n```biofuns\n<BIOFUNS>\n```\n,帮我总结探索这些基因可能潜在的生物学角色,哪一些更加值得关注,提出您的见解。请不要分点,写成一段,请控制字数在250字左右。请输出纯文本,不包括markdown格式。这些Gene包括:<GENELIST>",
  getTitle_Chinese = "你是一名生物学家,我提供给你一段生物学文本\n``` text\n <PRERES> \n```\n,请帮我的生物学主题拟定一个title。请您只返回标题本身,只返回最优的唯一(注意是唯一)标题,请不要返回任何其他的内容.",
  getModuleSummary_English = "Here is a list of reorganized biological themes derived from enrichment analysis. As a professional biologist, please concisely summarize the overarching theme of these pathways, their biological implications, and possible interpretations. Provide evidence where necessary. Also, propose a representative name for this module. Use an academic writing style and limit the response to ~250 words. Output plain text only (no markdown). The included pathways are:<PATHWAYLIST>",
  getGeneSummary_English = "Here is a set of genes that may share similar biological functions. Previous analysis suggests their potential functions may include \n```biofuns\n<BIOFUNS>\n```\n. As a professional biologist, please help summarize their potential biological roles, highlight which ones may be most noteworthy, and provide your insights. Write in a single paragraph (~250 words) without bullet points. Output plain text only (no markdown). The genes include:<GENELIST>",
  getTitle_English = "You are a biologist. I will provide you with a piece of biological text. \n``` text\n <PRERES> \n```\n Please help me generate a title for the biological topic. Only return the title itself—choose the single best (and only one) title. Do not return any other content."
)

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
    res <- retry_function(x$chat, ntry = 3, delay = 3, Prompt, echo = "none")
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
  return(list(results = results,cluster_names = cluster_names))
}


summarize_genes <- function(x, y, chat, prompt_type = "English") {
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
egt_llm_summary <- function(x, chat){
  if (class(x) != "EnrichGT_obj") cli::cli_abort("Please run `egt_recluster_analysis()` before summarizing. ")
  a1 <- summarize_clusters(x, chat)
  a2 <- summarize_genes(x, a1, chat)
  obj_llm <- new("egt_llm")
  obj_llm@pathways <- a1
  obj_llm@genes_and_title <- a2
  x@LLM_Annotation <- obj_llm
  return(x)
}
