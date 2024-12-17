library(dplyr)
library(tibble)
library(org.Hs.eg.db)
library(gt)
library(testthat)
library(withr)
library(EnrichGT)
library(readr)
DEGexample <- read_csv("tests/testthat/DEGexample.csv")
DEGexample_UpReg <- DEGexample |> dplyr::filter(pvalue<0.05,log2FoldChange>0.7)
dbgobp <- database_GO_BP(org.Hs.eg.db)
ora_result <- egt_enrichment_analysis(genes = DEGexample_UpReg$...1,database = dbgobp)
re_enrich <- egt_recluster_analysis(ora_result)
re_enrich
egt_plot_results(re_enrich)
egt_plot_umap(re_enrich)
gsea_result <- egt_gsea_analysis(genes_with_weights(DEGexample$...1,DEGexample$log2FoldChange),database = dbgobp)
re_enrich2 <- egt_recluster_analysis(gsea_result,ClusterNum = 30)
re_enrich2
egt_plot_results(re_enrich)
egt_plot_umap(re_enrich)


