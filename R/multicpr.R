#' @importFrom R6 R6Class
egt_mcp <- R6::R6Class(
  "comparison_reactor",
  public = list(
    initialize = function(Type = NULL) {
      if (Type %in% c("ORA", "ora")) {
        private$type <- "ORA"
        private$group_name <- tibble::tibble(
          Group = "Initial",
          Index = "G_0",
          Included = F
        )
      } else if (Type %in% c("GSEA", "gsea")) {
        private$type <- "GSEA"
        private$group_name <- tibble::tibble(
          Group = "Initial",
          Index = "G_0",
          Included = F
        )
      } else {
        cli::cli_abort("Please define correct analysis type: `ORA` or `GSEA`. ")
      }
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
      tdf <- cpres_internal_getter(x) |>
        dplyr::filter(!duplicated(Description)) |>
        dplyr::mutate(p_Adj = -log10(p.adjust))
      if (nrow(tdf) < 10) {
        cli::cli_abort("Invalid input. ")
      }
      tdf$CPRID <- tdf$Description
      bb[[new_info$Index[1]]] <- tdf

      private$group_name <- aa
      private$raw_enriched_result <- bb
      invisible(self)
    },
    make_plans = function(x = "auto", use_value = "p") {
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
      if (use_value == "NES" & private$type == "GSEA") {
        use_what0 <- "NES"
      } else {
        use_what0 <- "p_Adj"
      }
      a1 <- lapply(x, function(q) {
        d1 <- list0[[q]] |> dplyr::select(CPRID, !!sym(use_what0))
        colnames(d1)[2] <- paste0(q, "|", use_value)
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
    prefilter_by_NES = function(x = 1) {
      private$check_appended_data()
      if (private$type != "GSEA") cli::cli_abort("Not GSEA! ")
      private$raw_enriched_result -> list0
      list1 <- lapply(list0, function(q) q |> dplyr::filter(abs(NES) > x))
      names(list1) <- names(list0)
      private$raw_enriched_result <- list1
      cli::cli_alert_success(glue::glue("Filter according to NES > {x}"))
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
      pheatmap::pheatmap(
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
    }
  )
)
egt_comparison_reactor <- function(Type = NULL) {
  if (is.null(Type) | sum(Type %in% c("ORA", "ora", "GSEA", "gsea")) != 1) {
    cli::cli_alert_danger(
      "When create a comparison analysis object, please pinpoint the type of this reactor: "
    )
    cli::cli_code(
      "comparison_object <- egt_comparison_reactor('ORA')\ncomparison_object <- egt_comparison_reactor('GSEA')"
    )
    cli::cli_abort("Undefined reactor type")
  } else {
    kk <- egt_mcp$new(Type)
  }
  kk
}

# kk <- egt_comparison_reactor("ORA")
# kk$append_enriched_result(ora_result1, "group_1_go")
# kk$append_enriched_result(ora_result3, "group_2_go")
# kk$append_enriched_result(ora_result3, "group_3_go")
# kk$prefilter_by_p_adj(0.05)
# kk$make_plans()
# kk$make_plans(c("group_1_go","group_3_go"))
# kk$find_relationship(3)
# kk$fetch_relationship()

# a20 <- a2
