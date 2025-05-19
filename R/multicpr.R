
#'
#' Framework for comparing and analyzing enrichment results
#' across multiple groups
#' 
#' @details
#' Type should be one of "ORA" or "GSEA". 
#' For functions inside reactor, please see below 'See also' (in the bottom of this doc)
#' 
#' @examples
#' \dontrun{
#' # ORA example
#' reactor <- egt_comparison_reactor("ORA")
#' reactor$append_enriched_result(ora_result1, "group_1_go")
#' reactor$append_enriched_result(ora_result3, "group_3_go")
#' reactor$prefilter_by_p_adj(0.05)
#' reactor$make_plans(group = c("group_1_go","group_3_go"), use_value = "p")
#' reactor$find_relationship(Num = 3)
#' wordcloudFigure <- reactor$fetch_biological_theme()
#' wordcloudFigure[[1]]
#' reactor$split_by_cluster()
#' reactor$do_recluster(ClusterNum = 10)
#' cls_res <- reactor$get_recluster_result()
#' relation_df <- reactor$fetch_relationship()
#' 
#' # GSEA can do more with: 
#' gsea_reactor <- egt_comparison_reactor("GSEA")
#' gsea_reactor$prefilter_by_NES(1.5)
#' }
#' 
#' 
#' @seealso 
#' \code{\link{comparison_reactor_base}} for the base class documentation
#' \code{\link{comparison_reactor_ora}} for ORA-specific functionality
#' \code{\link{comparison_reactor_gsea}} for GSEA-specific functionality
#' @importFrom cowplot plot_grid
#' @importFrom R6 R6Class
#' @export

egt_comparison_reactor <- function(Type = NULL) {
  error00 <- function() {
    cli::cli_alert_danger(
      "When create a comparison analysis object, please pinpoint the type of this reactor: "
    )
    cli::cli_code(
      "comparison_object <- egt_comparison_reactor('ORA')\ncomparison_object <- egt_comparison_reactor('GSEA')"
    )
    cli::cli_abort("Undefined reactor type")
  }
  if (is.null(Type) | sum(Type %in% c("ORA", "ora", "GSEA", "gsea")) != 1) {
    error00()
  } else {
    if (Type == "ORA" | Type == "ora") {
      kk <- comparison_reactor_ora$new()
    } else if (Type == "GSEA" | Type == "gsea") {
      kk <- comparison_reactor_gsea$new()
    } else {
      error00()
    }
  }
  kk
}

#' Comparison Reactor Base
#' 
#' @description
#' Typically created via \code{\link{egt_comparison_reactor}()}
#' A framework for comparing and analyzing enrichment results across multiple groups.
#' This is the base class that provides core functionality for both ORA and GSEA analysis.
#'
#' @details
#' The comparison reactor allows:
#' - Appending multiple enrichment results from different groups
#' - Filtering results by p-value or NES score
#' - Comparing results between groups
#' - Identifying relationships between enriched terms
#' - Visualizing biological themes
#' - Sub-clustering results for deeper analysis
comparison_reactor_base <- R6::R6Class(
  "comparison_reactor_base",
  public = list(
    #' @description
    #' Create a new comparison reactor base object
    #' @param Type The analysis type ("ORA" or "GSEA")
    initialize = function(Type = NULL) {
      cli::cli_h1("EnrichGT comparison reactor")
      private$group_name <- tibble::tibble(
        Group = "Initial",
        Index = "G_0",
        Included = F
      )
    },
    #' @description
    #' Add enrichment results to the reactor
    #' @param x Data.frame of enrichment results (must contain p-value/adjusted p-value columns)
    #' @param group Character string naming the group for these results
    #' @return The reactor object (invisible) for method chaining
    #' @examples
    #' \dontrun{
    #' reactor <- egt_comparison_reactor("ORA")
    #' reactor$append_enriched_result(ora_result1, "group1")
    #' }
    append_enriched_result = function(x, group) {
      if (length(group) > 1) {
        cli::cli_abort("Please only include one group per time")
      }
      if (group %in% c("Initial", "auto")) {
        cli::cli_abort("Please change the group name.")
      }
      now_length <- nrow(private$group_name)
      new_info <- tibble::tibble(
        Group = group,
        Index = glue::glue("G_{now_length}"),
        Included = F
      )
      private$group_name -> aa
      if (group %in% aa$Group) {
        cli::cli_abort("Please do not use duplicated group name. ")
      }
      aa <- rbind(aa, new_info)

      private$raw_enriched_result -> bb
      tdf <- cpres_internal_getter(x)
      tdf <- .egt_mcp_helper_wash_p_val(tdf, type = "default")
      if (nrow(tdf) < 10) {
        cli::cli_abort("Invalid input. ")
      }
      tdf$CPRID <- tdf$Description
      
      private$latest_added_item_names -> latestTerm
      if(length(latestTerm) > 0){
        overlapNum <- sum(tdf$CPRID %in% latestTerm) 
        latestTotal <- length(latestTerm)
        pctOverlap <- round(((overlapNum/latestTotal) * 100), 2)
        cli::cli_alert_info(glue::glue("Overlap rate of new added data and the latest data:{pctOverlap}. Please ensure there are overlaps among appended data. "))
        if(pctOverlap < 25){
          cli::cli_alert_danger("Less than 25% overlap terms. Please ensure you are using same source and same param in these series of data")
        }
        private$latest_added_item_names <- latestTerm
      } else {
        private$latest_added_item_names <- latestTerm
      }
      

      bb[[new_info$Index[1]]] <- tdf

      private$group_name <- aa
      private$raw_enriched_result <- bb

      cli::cli_alert_success(glue::glue("Appended data into group {group}."))
      invisible(self)
    },
    #' @description 
    #' Print summary of groups in reactor
    summarize = function() {
      tbl0 <- private$group_name
      tbl0 <- tbl0 |> dplyr::filter(Group != "Initial")
      print(tbl0)
      cli::cli_alert_info(glue::glue(
        "Will make comparasion of multiple {private$type} results. "
      ))
    },
    #' @description
    #' Filter enrichment results by adjusted p-value cutoff
    #' @param x Numeric cutoff for adjusted p-values (default 0.05)
    #' @return The reactor object (invisible) for method chaining
    #' @examples
    #' \dontrun{
    #' reactor$prefilter_by_p_adj(0.01) # Use stricter cutoff
    #' }
    prefilter_by_p_adj = function(x = 0.05) {
      private$check_appended_data()
      private$raw_enriched_result -> list0
      list1 <- lapply(list0, function(q) q |> dplyr::filter(p.adjust < x))
      names(list1) <- names(list0)
      private$raw_enriched_result <- list1
      cli::cli_alert_success(glue::glue("Filter according to p adjust < {x}"))
      invisible(self)
    },
    #' @description
    #' Identify relationships between enriched terms across groups
    #' @param Num Number of top terms to consider from each group
    #' @param dist_method Distance method for clustering (default "euclidean")
    #' @param hclust_method Hierarchical clustering method (default "ward.D2")
    #' @param ... Additional parameters passed to heatmap function
    #' @return The reactor object (invisible) for method chaining
    #' @examples
    #' \dontrun{
    #' reactor$find_relationship(Num = 5, dist_method = "manhattan")
    #' }
    find_relationship = function(
      Num,
      dist_method = "euclidean",
      hclust_method = "ward.D2",
      ...
    ) {
      # Get data from agg_df
      mat <- as.matrix(private$agg_df[, -1]) # Remove CPRID column
      rownames(mat) <- private$agg_df$CPRID
      mat[is.na(mat)] <- 0 # Replace NA with 0

      # Perform hierarchical clustering
      d <- proxy::dist(mat, method = dist_method)
      hc <- stats::hclust(d, method = hclust_method)

      # Save clustering results
      clusters <- stats::cutree(hc, k = Num)
      private$cluster_results <- list(
        clusters = clusters,
        dist_method = dist_method,
        hclust_method = hclust_method,
        num_clusters = Num
      )
      private$dendrogram <- stats::as.dendrogram(hc)

      # Print cluster info
      cli::cli_alert_info(glue::glue(
        "Clustering with {dist_method} distance and {hclust_method} linkage, cut into {Num} clusters"
      ))
      print(table(clusters))

      base_pheatmap(
        mat,
        cluster_rows = hc,
        cutree_rows = Num,
        show_rownames = F,
        ...
      )
      .egt_mcp_helper_best_cluster(clusters)
      invisible(self)
    },
    #' @description
    #' Retrieve relationship data between terms
    #' @return Data.frame containing term relationships and cluster assignments
    #' @examples
    #' \dontrun{
    #' relation_df <- reactor$fetch_relationship()
    #' head(relation_df)
    #' }
    fetch_relationship = function() {
      if (is.null(private$cluster_results)) {
        cli::cli_abort("Please run find_relationship() first")
      }

      # Create data frame with cluster info
      cluster_df <- data.frame(
        CPRID = private$agg_df$CPRID,
        Cluster = private$cluster_results$clusters,
        stringsAsFactors = FALSE
      )
      # Merge with original data
      result_df <- base::merge.data.frame(
        cluster_df,
        private$agg_df,
        by = "CPRID",
        all.x = TRUE
      )
      result_df <- tibble::as_tibble(result_df)

      return(result_df)
    },
    #' @description
    #' Generate wordcloud visualization of biological themes
    #' @param ... Additional parameters passed to wordcloud generator
    #' @return List containing ggplot2 wordcloud objects
    #' @examples
    #' \dontrun{
    #' wordclouds <- reactor$fetch_biological_theme()
    #' wordclouds[[1]] # View first wordcloud
    #' }
    fetch_biological_theme = function(...) {
      df_export <- self$fetch_relationship()
      keep_CLS <- .egt_mcp_helper_best_cluster(df_export$Cluster)

      df_export2 <- df_export |> dplyr::filter(Cluster %in% keep_CLS)
      df_export2_list <- split(df_export2$CPRID, df_export2$Cluster)
      lapply(df_export2_list, function(x) {
        if (length(x) <= 3) {
          cli::cli_abort("Too small population. You can view data directly")
        }
      })
      figlist <- lapply(df_export2_list, wordcloud_generator2)
      figlist2 <- lapply(figlist, function(x) ggplotGrob(x))
      figout <- cowplot::plot_grid(
        plotlist = figlist2,
        ncol = 1,
        align = "h"
      )
      print(figout)
      invisible(figlist)
    },
    #' @description
    #' Split results by identified clusters for further analysis
    #' @param ... Additional parameters
    #' @return The reactor object (invisible) for method chaining
    split_by_cluster = function(...) {
      df_export <- self$fetch_relationship()
      keep_CLS <- .egt_mcp_helper_best_cluster(df_export$Cluster)
      df_export <- df_export |> dplyr::filter(Cluster %in% keep_CLS)
      aalist <- private$raw_enriched_result
      finalList <- list()
      for (i in names(table(df_export$Cluster))) {
        sub_cls_df <- df_export |> dplyr::filter(Cluster %in% i)
        if (nrow(sub_cls_df) <= 3) {
          next
        }
        for (j in names(aalist)) {
          bb0 <- aalist[[j]] |> dplyr::filter(CPRID %in% sub_cls_df$CPRID)
          if (nrow(bb0) <= 5) {
            next
          }
          cc0 <- private$group_name
          j1 <- cc0[cc0$Index == j, "Group", drop = T]
          finalList[[glue::glue("Data#{j1}|Cluster#{i}")]] <- bb0
        }
      }
      private$final_result <- finalList
      invisible(self)
    },
    #' @description
    #' Get list of results split by cluster
    #' @return List of data frames split by cluster
    get_splited_list = function() {
      if (length(private$result) < 0) {
        cli::cli_abort(
          "Please run previous analysis. This should run after `split_by_cluster()` already executed. "
        )
      }
      private$final_result
    },
    #' @description
    #' Perform sub-clustering within existing clusters
    #' @param ClusterNum Number of sub-clusters to generate
    #' @param P.adj Adjusted p-value cutoff (default 0.05)
    #' @param force Whether to force reclustering (default FALSE)
    #' @param nTop Number of top terms to consider (default 10)
    #' @param method Clustering method (default "ward.D2")
    #' @param ... Additional parameters passed to clustering functions
    #' @return The reactor object (invisible) for method chaining
    #' @examples
    #' \dontrun{
    #' reactor$do_recluster(ClusterNum = 5, method = "complete")
    #' }
    do_recluster = function(
      ClusterNum = 10,
      P.adj = 0.05,
      force = F,
      nTop = 10,
      method = "ward.D2",
      ...
    ) {
      if (length(private$final_result) < 0) {
        cli::cli_abort(
          "Please run previous analysis. This should run after `split_by_cluster()` already executed. "
        )
      }
      resList <- private$final_result
      resList2 <- lapply(
        resList,
        function(x)
          tryCatch(
            {
              egt_recluster_analysis(
                x = x,
                ClusterNum = ClusterNum,
                P.adj = P.adj,
                force = force,
                nTop = nTop,
                method = method,
                ...
              )
            },
            error = function(e) list()
          )
      )
      names(resList2) <- names(resList)
      private$recluster_result <- resList2
      invisible(self)
    },
    #' @description
    #' Retrieve sub-clustering results
    #' @return List containing reclustering results
    #' @examples
    #' \dontrun{
    #' recluster_results <- reactor$get_recluster_result()
    #' names(recluster_results)
    #' }
    get_recluster_result = function() {
      if (length(private$recluster_result) < 0) {
        cli::cli_abort(
          "Please run previous analysis. This should run after `do_recluster()` already executed. "
        )
      }
      private$recluster_result
    },
    print = function(...) {
      self$summarize()
    }
  ),
  private = list(
    group_name = tibble::tibble(),
    plan = c(),
    raw_enriched_result = list(),
    agg_df = data.frame(),
    cluster_results = list(),
    summary_stats = data.frame(),
    type = "Unknown",
    dendrogram = list(),
    final_result = list(),
    recluster_result = list(),
    latest_added_item_names = c(),
    check_appended_data = function() {
      if (sum(private$group_name$Group != "Initial") < 2) {
        cli::cli_abort(
          "Please add data into reactor by using YOUR_reactor_object$append_enriched_result(x = YOUR RESULT, group = FROM WHICH GROUP)"
        )
      }
    },
    make_plans_internal = function(x, use_value) {
      private$check_appended_data()
      if (x[1] == "auto" & length(x) == 1) {
        private$plan <- private$group_name$Group
        private$plan <- private$plan[private$plan != "Initial"]
      } else if (sum(x %in% private$group_name$Group) == length(x)) {
        private$plan <- x
      } else {
        cli::cli_abort("Please re-check your plan. Can't find enough groups. ")
      }
      x <- private$group_name[
        private$group_name$Group %in% private$plan,
        "Index",
        drop = T
      ]
      private$raw_enriched_result -> list0
      use_what0 <- .egt_mcp_helper_get_col(
        use_value = use_value,
        type = private$type
      )
      a1 <- lapply(x, function(q) {
        if (sum(colnames(list0[[q]]) %in% use_what0) < 1)
          cli::cli_abort(glue::glue(
            "Please check the input, can't find column: {use_what0}"
          ))
        d1 <- list0[[q]] |> dplyr::select(CPRID, !!dplyr::sym(use_what0))
        colnames(d1)[2] <- paste0(q, "|", use_what0)
        d1
      })
      # print(a1)
      a2 <- Reduce(
        function(a, b) {
          a <- as.data.frame(a)
          b <- as.data.frame(b)
          dplyr::full_join(a, b, by = "CPRID")
        },
        a1
      )
      a2[is.na(a2)] <- 0 # may have bugs
      private$agg_df <- a2
      # invisible(self)
    }
  )
)

#' ORA Comparison Reactor
#' 
#' @description
#' Typically created via \code{\link{egt_comparison_reactor}("ORA")}
#' 
#' @seealso \code{\link{comparison_reactor_base}} for inherited methods
#' A specialized reactor for comparing Over-Representation Analysis (ORA) results.
#' Inherits from comparison_reactor_base and provides ORA-specific functionality.
#'
#' @details
#' This reactor is optimized for comparing ORA results across multiple groups,
#' with methods tailored for p-value based comparisons.
comparison_reactor_ora <- R6::R6Class(
  "comparison_reactor_ora",
  inherit = comparison_reactor_base,
  public = list(
    #' @description
    #' Create a new ORA comparison reactor
    initialize = function() {
      super$initialize()
      private$type <- "ORA"
    },
    #' @description
    #' Create comparison plans between specified ORA groups
    #' @param group Character vector of group names to compare or "auto" for all groups
    #' @param use_value Which value to use for comparison ("p" for p-value or "padj" for adjusted p-value)
    #' @return The reactor object (invisible) for method chaining
    #' @examples
    #' \dontrun{
    #' ora_reactor$make_plans(group = c("group1", "group2"), use_value = "padj")
    #' }
    make_plans = function(group = "auto", use_value = "p") {
      group -> x
      private$make_plans_internal(x = x, use_value = use_value)
      invisible(self)
    }
  ),
  private = list()
)

#' GSEA Comparison Reactor
#' 
#' @description
#' Typically created via \code{\link{egt_comparison_reactor}("GSEA")}
#' 
#' @seealso \code{\link{comparison_reactor_base}} for inherited methods
#' A specialized reactor for comparing Gene Set Enrichment Analysis (GSEA) results.
#' Inherits from comparison_reactor_base and provides GSEA-specific functionality.
#'
#' @details
#' This reactor is optimized for comparing GSEA results across multiple groups,
#' with methods tailored for NES (Normalized Enrichment Score) based comparisons.
comparison_reactor_gsea <- R6::R6Class(
  "comparison_reactor_gsea",
  inherit = comparison_reactor_base,
  public = list(
    #' @description
    #' Create a new GSEA comparison reactor
    initialize = function() {
      super$initialize()
      private$type <- "GSEA"
    },
    #' @description
    #' Create comparison plans between specified GSEA groups
    #' @param group Character vector of group names to compare or "auto" for all groups
    #' @param use_value Which value to use for comparison ("NES" for normalized enrichment score)
    #' @return The reactor object (invisible) for method chaining
    #' @examples
    #' \dontrun{
    #' gsea_reactor$make_plans(group = c("group1", "group2"))
    #' }
    make_plans = function(group = "auto", use_value = "NES") {
      group -> x
      private$make_plans_internal(x = x, use_value = use_value)
      invisible(self)
    },
    #' @description
    #' Filter GSEA results by normalized enrichment score cutoff
    #' @param x Numeric cutoff for absolute NES values (default 1)
    #' @return The reactor object (invisible) for method chaining
    #' @examples
    #' \dontrun{
    #' gsea_reactor$prefilter_by_NES(1.5) # Filter for stronger effects
    #' }
    prefilter_by_NES = function(x = 1) {
      private$check_appended_data()
      if (private$type != "GSEA") cli::cli_abort("Not GSEA! ")
      private$raw_enriched_result -> list0
      list1 <- lapply(list0, function(q) q |> dplyr::filter(abs(NES) > x))
      names(list1) <- names(list0)
      private$raw_enriched_result <- list1
      cli::cli_alert_success(glue::glue("Filter according to NES > {x}"))
      invisible(self)
    }
  ),
  private = list()
)

.egt_mcp_helper_wash_p_val <- function(x, type) {
  if (type == "default") {
    x <- x |>
      dplyr::filter(!duplicated(Description)) |>
      dplyr::mutate(p_Adj = -log10(p.adjust))
  }
  x
  # Can add more filtering way
}

.egt_mcp_helper_get_col <- function(use_value, type) {
  if (use_value == "NES" & type == "GSEA") {
    use_what0 <- "NES" # Already contains
  } else if (sum(use_value %in% c("p", "pval", "P", "Padj", "pAdj")) > 0) {
    use_what0 <- "p_Adj" # Created in .egt_mcp_helper_wash_p_val fun
  }
  # This should return the agg_df 's true contained row names
  use_what0
}

.egt_mcp_helper_best_cluster <- function(clusters) {
  cluster_count <- table(clusters)
  keep_CLS <- names(cluster_count)[cluster_count >= 5]
  cli::cli_alert_success(glue::glue(
    "Suggest include: {paste(keep_CLS,collapse = ',')}"
  ))
  excluded_CLS <- names(cluster_count)[cluster_count < 5]
  cli::cli_alert_danger(glue::glue(
    "Suggest exclude: {paste(excluded_CLS,collapse = ',')}"
  ))
  invisible(keep_CLS)
}


base_pheatmap <- function(
  mat,
  cluster_rows = NULL,
  cluster_cols = FALSE,
  cutree_rows = NULL,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = RColorBrewer::brewer.pal(11, "Spectral"),
  legend = TRUE,
  legend_pos = "topright",
  ...
) {
  if (!is.null(cluster_rows) && inherits(cluster_rows, "hclust")) {
    rdend <- stats::as.dendrogram(cluster_rows)
  } else {
    rdend <- NA
  }
  if (!is.null(cutree_rows) && inherits(cluster_rows, "hclust")) {
    row_groups <- stats::cutree(cluster_rows, k = cutree_rows)
    side_colors <- cocol(cutree_rows, 3)[row_groups]
  } else {
    side_colors <- NULL
  }
  labRow <- if (show_rownames) rownames(mat) else NA
  labCol <- if (show_colnames) colnames(mat) else NA
  cols <- grDevices::colorRampPalette(color)(100) |> rev()
  colsA <- cols[1]
  colsB <- cols[length(cols)]
  # cols <- c(rep(colsA, 1), cols, rep(colsB, 1)) |> rev()
  brks <- seq(from = min(mat), to = max(mat), length.out = length(cols) + 1)
  stats::heatmap(
    mat,
    Rowv = rdend,
    Colv = if (cluster_cols) NULL else NA,
    labRow = labRow,
    labCol = labCol,
    col = cols,
    breaks = brks,
    RowSideColors = side_colors,
    scale = "none",
    margins = c(5, 5),
    cexRow = 0.6,
    cexCol = 0.8,
    ...
  )
  if (legend && !is.null(cutree_rows) && inherits(cluster_rows, "hclust")) {
    graphics::legend(
      legend_pos,
      legend = paste0("Cluster ", seq_len(cutree_rows)),
      fill = cocol(cutree_rows, 3),
      border = NA,
      bty = "n",
      xpd = T,
      cex = 0.5
    )
  }
}


#' @importFrom ggplot2 ggplot theme_minimal element_rect scale_color_gradient
#' @importFrom ggwordcloud geom_text_wordcloud_area
wordcloud_generator2 <- function(
  sentences,
  max_words = 100,
  remove_stop = TRUE,
  stopwords = NULL,
  min_freq = 1,
  bg_color = "white",
  ...
) {
  txt <- tolower(paste(sentences, collapse = " "))
  txt <- gsub("[[:punct:]]|[[:digit:]]", " ", txt)
  words <- unlist(strsplit(txt, "\\s+"))
  if (remove_stop) {
    if (is.null(stopwords)) {
      stopwords <- c(
        stp_wds(),
        "system",
        "positive",
        "negative",
        "cell",
        "pathway",
        "regulation",
        "type",
        "c",
        "p",
        "pp",
        "k",
        "ier",
        "l",
        "t",
        "s",
        "g",
        "ii",
        "via",
        "cell",
        "b",
        "nf",
        "cd",
        "like"
      )
    }
    words <- words[!words %in% stopwords]
  }
  words <- words[nzchar(words)]
  freq <- sort(table(words), decreasing = TRUE)
  freq <- freq[freq >= min_freq]
  if (length(freq) > max_words) freq <- head(freq, max_words)
  df <- data.frame(word = names(freq), freq = as.numeric(freq))

  ggplot(df, aes(label = word, size = freq, colour = freq)) +
    ggwordcloud::geom_text_wordcloud_area(
      rm_outside = TRUE,
      shape = "circle",
      eccentricity = 0.6
    ) +
    scale_size_area(max_size = 12) +
    scale_color_gradient(low = "steelblue", high = "darkred") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = bg_color, colour = NA),
      legend.position = "none",
      ...
    )
}


stp_wds <- function() {
  data(stopWords_list, package = "EnrichGT", envir = db_getter_env)
  db_getter_env$stopWords_list
}
