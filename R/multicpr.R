# egt_mcp <- R6Class(
#   "egt_multi_compare",
#   public = list(
#     initialize = function(Type = NULL) {
#       if (Type %in% c("ORA", "ora")) {
#         private$type <- "ORA"
#         private$group_name <- tibble::tibble(
#           Group = "Initial",
#           Index = "G_0",
#           Included = F
#         )
#       } else if (Type %in% c("GSEA", "gsea")) {
#         private$type <- "GSEA"
#         private$group_name <- tibble::tibble(
#           Group = "Initial",
#           Index = "G_0",
#           Included = F
#         )
#       } else {
#         cli::cli_abort("Please define correct analysis type: `ORA` or `GSEA`. ")
#       }
#     },
#     append_enriched_result = function(x, group) {
#       if (length(group) > 1) {
#         cli::cli_abort("Please only include one group per time")
#       }
#       if (group %in% c("Initial", "auto")) {
#         cli::cli_abort("Please change the group name.")
#       }
#       now_length <- nrow(private$group_name) - 1
#       new_info <- tibble::tibble(
#         Group = "Initial",
#         Index = glue::glue("G_{now_length}"),
#         Included = F
#       )
#       private$group_name -> aa
#       aa <- rbind(aa, new_info)

#       private$raw_enriched_result -> bb
#       tdf <- cpres_internal_getter(x) |>
#         dplyr::filter(!duplicated(Description)) |>
#         dplyr::mutate(p_Adj = -log10(p.adjust))
#       if (nrow(tdf) < 10) {
#         cli::cli_abort("Invalid input. ")
#       }
#       tdf$CPRID <- tdf$Description
#       bb[[new_info$Index[1]]] <- tdf

#       private$group_name <- aa
#       private$raw_enriched_result <- bb
#       invisible(self)
#     },
#     make_plans = function(x = "auto", use_value = "p") {
#       private$check_appended_data()
#       if (x == "auto") {
#         private$plan <- private$group_name$Group
#       } else if (sum(x %in% private$group_name$Group) == length(x)) {
#         private$plan <- x
#       } else {
#         cli::cli_abort("Please re-check your plan. Can't find enough groups. ")
#       }
#       private$raw_enriched_result -> list0
#       if (use_value == "NES" & private$type == "GSEA") {
#         use_what0 <- "NES"
#       } else {
#         use_what0 <- "p_Adj"
#       }
#       a1 <- lapply(x, function(q) {
#         d1 <- list0[[q]] |> dplyr::select(CPRID, !!sym(use_what0))
#         colnames(d1)[2] <- paste0(c(x), "|", use_value)
#         d1
#       })
#       a2 <- Reduce(
#         function(a, b) {
#           a <- as.data.frame(a)
#           b <- as.data.frame(b)
#           dplyr::full_join(a, b, by = "CPRID")
#         },
#         a1
#       )
#       private$agg_df <- a2
#       invisible(self)
#     },
#     summarize = function() {
#       tbl0 <- private$group_name
#       tbl0 <- tbl0 |> dplyr::filter(Group != "Initial")
#       print(tbl0)
#       cli::cli_alert_info(glue::glue(
#         "Will make comparasion of multiple {private$type} results. "
#       ))
#     },
#     prefilter_by_p_adj = function(x) {
#       private$check_appended_data()
#       private$raw_enriched_result -> list0
#       list1 <- lapply(list0, function(x) x |> dplyr::filter(p.adjust < x))
#       names(list1) <- names(list0)
#       private$raw_enriched_result <- list1
#       cli::cli_alert_success(glue::glue("Filter according to p adjust < {x}"))
#       invisible(self)
#     },
#     prefilter_by_NES = function(x) {
#       private$check_appended_data()
#       if (private$type != "GSEA") cli::cli_abort("Not GSEA! ")
#       private$raw_enriched_result -> list0
#       list1 <- lapply(list0, function(x) x |> dplyr::filter(abs(NES) > x))
#       names(list1) <- names(list0)
#       private$raw_enriched_result <- list1
#       cli::cli_alert_success(glue::glue("Filter according to NES > {x}"))
#       invisible(self)
#     },
#     find_relationship = function(Num, dist_method = "euclidean", hclust_method = "complete", ...) {
#       # Get data from agg_df
#       mat <- as.matrix(private$agg_df[, -1])  # Remove CPRID column
#       rownames(mat) <- private$agg_df$CPRID
#       mat[is.na(mat)] <- 0  # Replace NA with 0
      
#       # Perform hierarchical clustering
#       d <- dist(mat, method = dist_method)
#       hc <- hclust(d, method = hclust_method)
      
#       # Save clustering results
#       clusters <- cutree(hc, k = Num)
#       private$cluster_results <- list(
#         clusters = clusters,
#         dist_method = dist_method,
#         hclust_method = hclust_method,
#         num_clusters = Num
#       )
#       private$dendrogram <- as.dendrogram(hc)
      
#       # Print cluster info
#       cli::cli_alert_info(glue::glue(
#         "Clustering with {dist_method} distance and {hclust_method} linkage, cut into {Num} clusters"
#       ))
#       print(table(clusters))
      
#       # Try pheatmap first
#       tryCatch({
#         pheatmap::pheatmap(
#           mat,
#           cluster_rows = hc,
#           cutree_rows = Num,
#           ...
#         )
#       }, error = function(e) {
#         # Fallback to ComplexHeatmap if pheatmap fails
#         if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
#           stop("Please install ComplexHeatmap package: BiocManager::install('ComplexHeatmap')")
#         }
#         ComplexHeatmap::Heatmap(
#           mat,
#           cluster_rows = hc,
#           row_split = Num,
#           ...
#         )
#       })
      
#       invisible(self)
#     },
#     fetch_relationship = function() {
#       if (is.null(private$cluster_results)) {
#         cli::cli_abort("Please run find_relationship() first")
#       }
      
#       # Create data frame with cluster info
#       cluster_df <- data.frame(
#         CPRID = private$agg_df$CPRID,
#         Cluster = private$cluster_results$clusters,
#         Distance_Method = private$cluster_results$dist_method,
#         Linkage_Method = private$cluster_results$hclust_method,
#         Num_Clusters = private$cluster_results$num_clusters,
#         stringsAsFactors = FALSE
#       )
      
#       # Merge with original data
#       result_df <- merge(
#         cluster_df,
#         private$agg_df,
#         by = "CPRID",
#         all.x = TRUE
#       )
      
#       return(result_df)
#     },
#     do_recluster = function(...) {
#     },
#     fetch_recluster_results = function() {
#     },
#     print = function(...) {
#       self$summarize()
#     }
#   ),
#   private = list(
#     group_name = tibble::tibble(),
#     plan = c(),
#     raw_enriched_result = list(),
#     agg_df = data.frame(),
#     cluster_results = list(),
#     summary_stats = data.frame(),
#     type = "Unknown",
#     check_appended_data = function() {
#       if (sum(private$group_name$Group != "Initial") < 2) {
#         cli::cli_abort(
#           "Please add data into reactor by using YOUR_reactor_object$append_enriched_result(x = YOUR RESULT, group = FROM WHICH GROUP)"
#         )
#       }
#     }
#   )
# )
