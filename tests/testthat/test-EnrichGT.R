library(clusterProfiler)
library(org.Hs.eg.db)
library(gt)
library(testthat)
library(withr)

test_that("EnrichGT creates four HTML files", {
  tmp_dir <- withr::local_tempdir()
  file1 <- file.path(tmp_dir, "test1.html")
  file2 <- file.path(tmp_dir, "test2.html")
  file3 <- file.path(tmp_dir, "test3.html")
  file4 <- file.path(tmp_dir, "test4.html")
  file5 <- file.path(tmp_dir, "test5.html")
  data(geneList, package="DOSE")
  gene <- names(geneList)[abs(geneList) > 2]
  gene1 <- names(geneList)[abs(geneList) > 1]
  gene1a<-gene1[sample(1000,400)]
  gene1b<-gene1[sample(1000,400)]
  ego <- enrichGO(gene          = gene,
                  universe      = names(geneList),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  egoa <- enrichGO(gene          = gene1a,
                  universe      = names(geneList),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  egob <- enrichGO(gene          = gene1b,
                  universe      = names(geneList),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  kk <- enrichKEGG(gene         = gene,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
  obj1<-EnrichGT(ego, ClusterNum = 15, P.adj = 1,method = "average")
  obj2<-EnrichGT(ego, ClusterNum = 15, P.adj = 1)
  obj1@gt_object |> gt::gtsave(file1)
  expect_true(file.exists(file1), info = "test1.html should be created")
  expect_true(obj1@clustering_tree[["method"]]=="average",info="method works!")
  expect_true(obj2@clustering_tree[["method"]]=="cosine",info="method works_2!")
  compareGT(egoa,egob)->dd
  expect_true(dd[["Overlap_Control"]]@clustering_tree[["dist.method"]]=="cosine",info = "cprgt works!")
  EnrichGT(kk, ClusterNum = 100, P.adj = 1)@gt_object |> gt::gtsave(file2)
  expect_true(file.exists(file2), info = "test2.html should be created")
  expect_true(nrow(EnrichGT(ego, ClusterNum = 10, P.adj = 1,nTop = 5)@enriched_result)==10,info="nTop.ora Works")
  ego3 <- gseGO(geneList     = geneList,
                OrgDb        = org.Hs.eg.db,
                ont          = "CC",
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
  EnrichGT(ego3, ClusterNum = 100, P.adj = 1)@gt_object |> gt::gtsave(file3)
  expect_true(nrow(EnrichGT(ego3, ClusterNum = 3, P.adj = 1,nTop = 2)@enriched_result)==6,info="nTop.gsea Works")
  expect_true(file.exists(file3), info = "test3.html should be created")
  data(gcSample)
  ck <- compareCluster(geneCluster = gcSample, fun = enrichGO, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.5, qvalueCutoff = 0.5,ont="BP")
  ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  sd <- EnrichGT(ck, ClusterNum = 3, P.adj = 1)
  das<-ck@compareClusterResult$Cluster |> table()

  dza<-sd[[names(das[das==max(das)])]]
  dza@gt_object |> gt::gtsave(file4)
  expect_true(file.exists(file4), info = "test4.html should be created")

  EnrichGT(list(ego,kk), ClusterNum = 15, P.adj = 1)@gt_object |> gt::gtsave(file5)
  expect_true(file.exists(file5), info = "test5.html should be created")
})
