egt_comparison_reactor <- function(Type = NULL) {
  error00 <- function(){
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
    if (Type == "ORA" | Type == "ora"){
      kk <- comparison_reactor_ora$new()
    } else if(Type == "GSEA" | Type == "gsea"){
      kk <- comparison_reactor_gsea$new()
    } else {
      error00()
    }
  }
  kk
}


#' @importFrom R6 R6Class
comparison_reactor_base <- R6::R6Class(
  "comparison_reactor_base",
  public = list(
    initialize = function(Type = NULL) {
      cli::cli_h1("EnrichGT comparison reactor")
      private$group_name <- tibble::tibble(
        Group = "Initial",
        Index = "G_0",
        Included = F
      )
    },
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
      bb[[new_info$Index[1]]] <- tdf

      private$group_name <- aa
      private$raw_enriched_result <- bb

      cli::cli_alert_success(glue::glue("Appended data into group {group}."))
      invisible(self)
    },
    summarize = function() {
      tbl0 <- private$group_name
      tbl0 <- tbl0 |> dplyr::filter(Group != "Initial")
      print(tbl0)
      cli::cli_alert_info(glue::glue(
        "Will make comparasion of multiple {private$type} results. "
      ))
    },
    prefilter_by_p_adj = function(x = 0.05) {
      private$check_appended_data()
      private$raw_enriched_result -> list0
      list1 <- lapply(list0, function(q) q |> dplyr::filter(p.adjust < x))
      names(list1) <- names(list0)
      private$raw_enriched_result <- list1
      cli::cli_alert_success(glue::glue("Filter according to p adjust < {x}"))
      invisible(self)
    },
    find_relationship = function(
      Num,
      dist_method = "euclidean",
      hclust_method = "complete",
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

      invisible(self)
    },
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
    do_recluster = function(...) {
      # todo
      NULL
    },
    fetch_recluster_results = function() {
      # todo
      NULL
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

comparison_reactor_ora <- R6::R6Class(
  "comparison_reactor_ora",
  inherit = comparison_reactor_base,
  public = list(
    initialize = function() {
      super$initialize()
      private$type <- "ORA"
    },
    make_plans = function(x = "auto", use_value = "p") {
      private$make_plans_internal(x = x, use_value = use_value)
      invisible(self)
    }
  ),
  private = list()
)

comparison_reactor_gsea <- R6::R6Class(
  "comparison_reactor_gsea",
  inherit = comparison_reactor_base,
  public = list(
    initialize = function() {
      super$initialize()
      private$type <- "GSEA"
    },
    make_plans = function(x = "auto", use_value = "NES") {
      private$make_plans_internal(x = x, use_value = use_value)
      invisible(self)
    },
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


base_pheatmap <- function(
  mat,
  cluster_rows    = NULL, 
  cluster_cols    = FALSE,
  cutree_rows     = NULL, 
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  color           = heat.colors(100),
  legend          = TRUE,
  legend_pos      = "topright",
  ...
){
  if (!is.null(cluster_rows) && inherits(cluster_rows, "hclust")) {
    rdend <- stats::as.dendrogram(cluster_rows)
  } else {
    rdend <- NA
  }
  if (!is.null(cutree_rows) && inherits(cluster_rows, "hclust")) {
    row_groups  <- stats::cutree(cluster_rows, k = cutree_rows)
    side_colors <- grDevices::rainbow(cutree_rows)[row_groups]
  } else {
    side_colors <- NULL
  }
  labRow <- if (show_rownames) rownames(mat) else NA
  labCol <- if (show_colnames) colnames(mat) else NA
  stats::heatmap(
    mat,
    Rowv          = rdend,
    Colv          = if (cluster_cols) NULL else NA,
    labRow        = labRow,
    labCol        = labCol,
    col           = color,
    RowSideColors = side_colors,
    ...
  )
  if (legend && !is.null(cutree_rows) && inherits(cluster_rows, "hclust")) {
    legend(
      legend_pos,
      legend = paste0("Cluster ", seq_len(cutree_rows)),
      fill   = rainbow(cutree_rows),
      border = NA,
      bty    = "n"
    )
  }
}

