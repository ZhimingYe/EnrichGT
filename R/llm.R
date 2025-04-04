promptList <- list(
  getModuleSummary_Chinese = "这是一些生物学富集分析后重新再次整理后得到的生物学主题列表。你是一名专业生物学家,帮我用一些话概括这些通路的主题是什么,反映了什么,有什么生物学解释的角度。请简洁有力,必要时提供依据。请您还需要为这个模块拟定一个代表性名字。请使用学术写作风格。请控制字数在250字左右,避免长篇大论。请输出纯文本,不包括markdown格式。参与这个模块的通路包括:<PATHWAYLIST>",
  getGeneSummary_Chinese = "这是一系列Gene,可能具有类似生物学功能。前面已经总结他们的生物学功能可能包括```biofuns\n<BIOFUNS>\n```\n,帮我总结探索这些基因可能潜在的生物学角色,哪一些更加值得关注,提出您的见解。请不要分点,写成一段,请控制字数在250字左右。请输出纯文本,不包括markdown格式。这些Gene包括:<GENELIST>",
  getModuleSummary_English = "Here is a list of reorganized biological themes derived from enrichment analysis. As a professional biologist, please concisely summarize the overarching theme of these pathways, their biological implications, and possible interpretations. Provide evidence where necessary. Also, propose a representative name for this module. Use an academic writing style and limit the response to ~250 words. Output plain text only (no markdown). The included pathways are:<PATHWAYLIST>",
  getGeneSummary_English = "Here is a set of genes that may share similar biological functions. Previous analysis suggests their potential functions may include ```biofuns\n<BIOFUNS>\n```\n. As a professional biologist, please help summarize their potential biological roles, highlight which ones may be most noteworthy, and provide your insights. Write in a single paragraph (~250 words) without bullet points. Output plain text only (no markdown). The genes include:<GENELIST>"
)

retry_function <- function(FUN, ntry = 3, delay = 3, ...) {
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
    "The LLM support of `EnrichGT` is based on package `ellmer`"
  )
  cli::cli_abort(
    "Please check package `ellmer` installation or your model provider & API key. We can't get any results from this `chat` object. Please check https://ellmer.tidyverse.org/index.html for your LLM model registration and provide it properly in EnrichGT. "
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
  results
}


summarize_genes <- function(x, y, chat, prompt_type = "English") {
  excludeIndex <- (-1)
  if (!inherits(x, "EnrichGT_obj")) {
    cli::cli_abort("Input must be an EnrichGT_obj object")
  }
  if (sum(sapply(y, is.null)) > 0) {
    excludeIndex <- which(sapply(y, is.null))
    x@gene_modules <- x@gene_modules[-excludeIndex]
    y <- y[-excludeIndex]
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
  need2summarize <- list()
  for (i in 1:length(y)) {
    prompt2 <- gsub(
      "<GENELIST>",
      paste(x@gene_modules[[i]], collapse = ", "),
      prompt
    )
    prompt2 <- gsub("<BIOFUNS>", y[[i]], prompt2)
    need2summarize[i] <- prompt2
  }
  cluster_names <- names(x@gene_modules)
  names(need2summarize) <- cluster_names
  results <- lapply(
    cli::cli_progress_along(cluster_names, name = "Summarizing genes"),
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
  results
}
