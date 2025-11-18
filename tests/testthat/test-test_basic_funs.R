library(dplyr)
library(tibble)
library(ggplot2)
if(!require(org.Hs.eg.db)){
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)
library(EnrichGT)
library(readr)
library(testthat)
test_that("Data filtering and basic analysis works", {
  data("DEGexample", package = "EnrichGT")
  DEGexample2 <- DEGexample |> dplyr::filter(pvalue < 0.05)
  expect_gt(nrow(DEGexample2), 0)
  DEGexample_UpReg <- DEGexample |> dplyr::filter(pvalue < 0.05, log2FoldChange > 0.7)
  expect_gt(nrow(DEGexample_UpReg), 0)
})

test_that("Enrichment analysis functions work", {
  data("DEGexample", package = "EnrichGT")
  DEGexample_UpReg <- DEGexample |> dplyr::filter(pvalue < 0.05, log2FoldChange > 0.7)
  DEGs <- DEGexample_UpReg$...1
  ora_result <- egt_enrichment_analysis(
    genes = DEGs,
    database = database_GO_BP(OrgDB = org.Hs.eg.db)
  )
  expect_gt(nrow(ora_result), 2)
  keggRes <- database_KEGG()
  expect_gt(nrow(keggRes), 2)
})

test_that("Recluster analysis works", {
  data("DEGexample", package = "EnrichGT")
  DEGexample_UpReg <- DEGexample |> dplyr::filter(pvalue < 0.05, log2FoldChange > 0.7)
  DEGs <- DEGexample_UpReg$...1
  ora_result <- egt_enrichment_analysis(
    genes = DEGs,
    database = database_GO_BP(OrgDB = org.Hs.eg.db)
  )
  re_enrichment_results <- egt_recluster_analysis(ora_result)
  expect_s4_class(re_enrichment_results, "EnrichGT_obj")
  expect_true("gt_tbl" %in% class(re_enrichment_results@gt_object))
})


test_that("Fusing analysis works", {
  data("DEGexample", package = "EnrichGT")
  DEGexample_UpReg <- DEGexample |> dplyr::filter(pvalue < 0.05, log2FoldChange > 0.7)
  DEGexample_DownReg <- DEGexample |> dplyr::filter(pvalue < 0.05, log2FoldChange < (-0.7))
  DEGs <- DEGexample_UpReg$...1
  ora_result1 <- egt_enrichment_analysis(
    genes = DEGs,
    database = database_GO_BP(OrgDB = org.Hs.eg.db)
  )
  ora_result2 <- egt_enrichment_analysis(
    genes = DEGs,
    database = database_Reactome(OrgDB = org.Hs.eg.db)
  )
  ora_result3 <- egt_enrichment_analysis(
    genes = DEGexample_DownReg$...1,
    database = database_Reactome(OrgDB = org.Hs.eg.db)
  )
  re_enrichment_results <- egt_recluster_analysis(list(ora_result1,ora_result2))
  re_enrichment_results
  re_enrichment_results2 <- egt_compare_groups(ora_result3,ora_result2)
  re_enrichment_results2
  expect_s4_class(re_enrichment_results, "EnrichGT_obj")
  expect_true("gt_tbl" %in% class(re_enrichment_results@gt_object))
  expect_s4_class(re_enrichment_results2[[1]], "EnrichGT_obj")
})


test_that("Visualization functions work", {
  data("DEGexample", package = "EnrichGT")
  DEGexample_UpReg <- DEGexample |> dplyr::filter(pvalue < 0.05, log2FoldChange > 0.7)
  DEGs <- DEGexample_UpReg$...1
  ora_result <- egt_enrichment_analysis(
    genes = DEGs,
    database = database_GO_BP(OrgDB = org.Hs.eg.db)
  )
  plt1 <- egt_plot_results(ora_result, showIDs = TRUE, ntop = 20)
  expect_s3_class(plt1, "gg")
  re_enrichment_results <- egt_recluster_analysis(ora_result)
  plt2 <- egt_plot_results(re_enrichment_results)
  expect_s3_class(plt2, "gg")
  plt3 <- egt_fetch_biological_theme(re_enrichment_results, 3)
  expect_s3_class(plt3, "gg")
})

test_that("GSEA analysis works", {
  data("DEGexample", package = "EnrichGT")
  DEGexample2 <- DEGexample |> dplyr::filter(pvalue < 0.05)
  genes_with_weights <- genes_with_weights(DEGexample2$...1, DEGexample2$log2FoldChange)
  GSEAexample <- egt_gsea_analysis(
    genes = genes_with_weights,
    database = database_Reactome(OrgDB = org.Hs.eg.db)
  )
  expect_gt(nrow(GSEAexample), 2)
  ll <- egt_plot_gsea(
    GSEAexample$Description[99],
    genes = genes_with_weights,
    database = database_Reactome(OrgDB = org.Hs.eg.db)
  )
  expect_s3_class(ll, "gg")
  ll2 <- egt_plot_gsea(
    GSEAexample[1:3, ],
    genes = genes_with_weights,
    database = database_Reactome(OrgDB = org.Hs.eg.db)
  )
  expect_s3_class(ll2, "gg")
  saveRDS(GSEAexample, file = "GSEAexample.rds")
  re_enrichment_results_gsea <- egt_recluster_analysis(GSEAexample, force = T)
  expect_s4_class(re_enrichment_results_gsea, "EnrichGT_obj")
})



data("MultiDEGExample", package = "EnrichGT")
ora_result_g1 <- egt_enrichment_analysis(
  genes = MultiDEGExample$liver |> dplyr::filter(P.Value < 0.05, logFC > 0.5) |> dplyr::pull(...1),
  database = database_GO_BP(OrgDB = org.Hs.eg.db)
)
ora_result_g2 <- egt_enrichment_analysis(
  genes = MultiDEGExample$kidney |> dplyr::filter(P.Value < 0.05, logFC > 0.5) |> dplyr::pull(...1),
  database = database_GO_BP(OrgDB = org.Hs.eg.db)
)
ora_result_g3 <- egt_enrichment_analysis(
  genes = MultiDEGExample$muscle |> dplyr::filter(P.Value < 0.05, logFC > 0.5) |> dplyr::pull(...1),
  database = database_GO_BP(OrgDB = org.Hs.eg.db)
)
ora_result_g4 <- egt_enrichment_analysis(
  genes = MultiDEGExample$pancreas |> dplyr::filter(P.Value < 0.05, logFC > 0.5) |> dplyr::pull(...1),
  database = database_GO_BP(OrgDB = org.Hs.eg.db)
)
ora_result_g5 <- egt_enrichment_analysis(
  genes = MultiDEGExample$spleen |> dplyr::filter(P.Value < 0.05, logFC > 0.5) |> dplyr::pull(...1),
  database = database_GO_BP(OrgDB = org.Hs.eg.db)
)


ora_result_gsea1 <- egt_gsea_analysis(
  genes = genes_with_weights(MultiDEGExample$liver$...1, MultiDEGExample$liver$logFC),
  database = database_GO_BP(OrgDB = org.Hs.eg.db)
)
ora_result_gsea2 <- egt_gsea_analysis(
  genes = genes_with_weights(MultiDEGExample$kidney$...1, MultiDEGExample$kidney$logFC),
  database = database_GO_BP(OrgDB = org.Hs.eg.db)
)
ora_result_gsea3 <- egt_gsea_analysis(
  genes = genes_with_weights(MultiDEGExample$muscle$...1, MultiDEGExample$muscle$logFC),
  database = database_GO_BP(OrgDB = org.Hs.eg.db)
)

ora_result_r1 <- egt_enrichment_analysis(
  genes = MultiDEGExample$liver |> dplyr::filter(P.Value < 0.05, logFC > 0.5) |> dplyr::pull(...1),
  database = database_Reactome(OrgDB = org.Hs.eg.db)
)
ora_result_r2 <- egt_enrichment_analysis(
  genes = MultiDEGExample$kidney |> dplyr::filter(P.Value < 0.05, logFC > 0.5) |> dplyr::pull(...1),
  database = database_Reactome(OrgDB = org.Hs.eg.db)
)
ora_result_r3 <- egt_enrichment_analysis(
  genes = MultiDEGExample$muscle |> dplyr::filter(P.Value < 0.05, logFC > 0.5) |> dplyr::pull(...1),
  database = database_Reactome(OrgDB = org.Hs.eg.db)
)
ora_result_r4 <- egt_enrichment_analysis(
  genes = MultiDEGExample$pancreas |> dplyr::filter(P.Value < 0.05, logFC > 0.5) |> dplyr::pull(...1),
  database = database_Reactome(OrgDB = org.Hs.eg.db)
)
ora_result_r5 <- egt_enrichment_analysis(
  genes = MultiDEGExample$spleen |> dplyr::filter(P.Value < 0.05, logFC > 0.5) |> dplyr::pull(...1),
  database = database_Reactome(OrgDB = org.Hs.eg.db)
)

ora_fuse1 <- egt_recluster_analysis(list(ora_result_g1, ora_result_r1))
ora_fuse2 <- egt_recluster_analysis(list(ora_result_g2, ora_result_r2))
ora_fuse3 <- egt_recluster_analysis(list(ora_result_g3, ora_result_r3))
ora_fuse4 <- egt_recluster_analysis(list(ora_result_g4, ora_result_r4))

test_that("Delete Function works", {
  ora_result_g1A <- ora_result_g1 %-delete->% "test"
  ora_result_g1A <- ora_fuse3 %-delete->% "test"
  expect_s4_class(ora_result_g1A, "EnrichGT_obj")
})

test_that("ORA not fused comparing", {
  reactor1 <- egt_comparison_reactor("ora")
  reactor1$append_enriched_result(ora_result_g1, "liver_GO")
  reactor1$append_enriched_result(ora_result_g2, "kidney_GO")
  reactor1$append_enriched_result(ora_result_g3, "muscle_GO")
  reactor1$append_enriched_result(ora_result_g4, "pancreas_GO")
  reactor1$append_enriched_result(ora_result_g5, "spleen_GO")
  expect_error(reactor1$append_enriched_result(ora_result_g1, "liver_GO"))
  reactor1$prefilter_by_p_adj(0.05)
  reactor1$make_plans()
  expect_equal((reactor1$.__enclos_env__$private$agg_df |> ncol()), 6)
  reactor1$make_plans(c("liver_GO","kidney_GO","spleen_GO"))
  expect_equal((reactor1$.__enclos_env__$private$agg_df |> ncol()), 4)
  reactor1$make_plans()
  reactor1$find_relationship(5)
  expect_equal((reactor1$fetch_relationship()$Cluster |> table() |> names() |> length()), 5)
  figlist0 <- reactor1$fetch_biological_theme()
  expect_s3_class(figlist0[[1]], "gg")
  reactor1$split_by_cluster()
  # expect_equal(length(names(reactor1$get_splited_list())), 5*5)
  reactor1$do_recluster()
  res <- reactor1$get_recluster_result()
  expect_s4_class(res[[1]], "EnrichGT_obj")

})




test_that("ORA have fused comparing", {
  reactor1 <- egt_comparison_reactor("ora")
  reactor1$append_enriched_result(ora_result_g1, "liver_F")
  reactor1$append_enriched_result(ora_result_g2, "kidney_F")
  reactor1$append_enriched_result(ora_result_g3, "muscle_F")
  expect_error(reactor1$append_enriched_result(ora_result_g1, "liver_F"))
  reactor1$make_plans()
  expect_equal((reactor1$.__enclos_env__$private$agg_df |> ncol()), 4)
  reactor1$prefilter_by_p_adj(0.05)
  reactor1$make_plans(c("liver_F","kidney_F","muscle_F"))
  expect_equal((reactor1$.__enclos_env__$private$agg_df |> ncol()), 4)
  reactor1$make_plans()
  reactor1$find_relationship(5)
  expect_equal((reactor1$fetch_relationship()$Cluster |> table() |> names() |> length()), 5)
  figlist0 <- reactor1$fetch_biological_theme()
  expect_s3_class(figlist0[[1]], "gg")
  reactor1$split_by_cluster()
  # expect_equal(length(names(reactor1$get_splited_list())), 3*5)
  reactor1$do_recluster()
  res <- reactor1$get_recluster_result()
  expect_s4_class(res[[1]], "EnrichGT_obj")

})



test_that("GSEA comparing", {
  reactor1 <- egt_comparison_reactor("gsea")
  reactor1$append_enriched_result(ora_result_gsea1, "liver_F")
  reactor1$append_enriched_result(ora_result_gsea2, "kidney_F")
  reactor1$append_enriched_result(ora_result_gsea3, "muscle_F")
  expect_error(reactor1$append_enriched_result(ora_result_g1, "liver_F"))
  reactor1$make_plans()
  expect_equal((reactor1$.__enclos_env__$private$agg_df |> ncol()), 4)
  reactor1$prefilter_by_p_adj(0.05)
  reactor1$prefilter_by_NES(1)
  reactor1$make_plans(c("liver_F","kidney_F","muscle_F"))
  expect_equal((reactor1$.__enclos_env__$private$agg_df |> ncol()), 4)
  reactor1$make_plans()
  reactor1$find_relationship(5)
  expect_equal((reactor1$fetch_relationship()$Cluster |> table() |> names() |> length()), 5)
  figlist0 <- reactor1$fetch_biological_theme()
  expect_s3_class(figlist0[[1]], "gg")
  reactor1$split_by_cluster()
  # expect_equal(length(names(reactor1$get_splited_list())), 3*5)
  # Because have no overlap, skip this check
  reactor1$do_recluster()
  res <- reactor1$get_recluster_result()
  expect_s4_class(res[[1]], "EnrichGT_obj")

})

