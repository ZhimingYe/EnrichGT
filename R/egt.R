setGeneric("doEnrichGT", function(x, ...) standardGeneric("doEnrichGT"))
setMethod("doEnrichGT", signature(x = "enrichResult"), function(x, ...) {
  if (sum(grepl("^GO", names(x@geneSets))) > 5) {
    nsimp()
  }
  y <- .genGT(x = x@result, ...)
  return(y)
})
setMethod(
  "doEnrichGT",
  signature(x = "compareClusterResult"),
  function(x, ...) {
    x <- x@compareClusterResult
    y <- .cprres(x, ...)
    if (sum(grepl("^GO", names(y$ID))) > 5) {
      nsimp()
    }
    return(y)
  }
)
setMethod("doEnrichGT", signature(x = "gseaResult"), function(x, ...) {
  y <- .genGSEAGT(x = x@result, ...)
  return(y)
})
setMethod("doEnrichGT", signature(x = "data.frame"), function(x, ...) {
  if ("NES" %in% colnames(x)) {
    y <- .genGSEAGT(x, ...)
    return(y)
  } else if ("Cluster" %in% colnames(x)) {
    y <- .cprres(x, ...)
    if (sum(grepl("^GO", names(y$ID))) > 5) {
      nsimp()
    }
    return(y)
  } else {
    if (sum(grepl("^GO", names(x$ID))) > 5) {
      nsimp()
    }
    y <- .genGT(x = x, ...)
    return(y)
  }
})
setMethod("doEnrichGT", signature(x = "list"), function(x, ...) {
  y <- lapply(x, function(z) {
    q <- cpres_internal_getter(z)
    return(q)
  })
  y2 <- do.call(rbind, y)
  y2 <- y2 |> as.data.frame()
  y2 <- y2[!duplicated(y2$ID), ]
  y2 <- y2[!duplicated(y2$Description), ]
  y2 <- as.data.frame(y2)
  y3 <- universalGT(y2, ...)
  return(y3)
})

.genGT <- function(
  x,
  ClusterNum,
  P.adj = 0.05,
  force = F,
  objname,
  nTop,
  method,
  ...
) {
  InnerDF <- x
  InnerDF -> InnerDF_000
  InnerDF <- InnerDF |> dplyr::filter(pvalue < 0.05)
  if (dim(x)[1] <= 5) {
    cli::cli_abort(
      "less than 5 enrichment result satisfied that the origin P<0.05(not adjusted)"
    )
  }
  .checkRowNames(x, "ORA")
  ClusterNum0 <- ClusterNum
  ClusterNum <- .genClusterNum(
    x = x,
    ClusterNum = ClusterNum0,
    force = force
  ) |>
    round()
  clsObj <- .enrichpws(InnerDF$ID, InnerDF$geneID, ClusterNum, method)
  suppressMessages({
    InnerDF <- InnerDF |> dplyr::left_join(clsObj[[1]])
  }) # Merge according to "ID"
  if (sum(colnames(InnerDF) == "Up_Vs_Down") > 0) {
    InnerDF <- InnerDF |>
      dplyr::filter(pvalue < 0.05, Count >= 5, p.adjust < P.adj) |>
      dplyr::mutate(
        geneID = paste0("@Up_Genes@: /", UpGenes, "/@Down_Genes@: /", DownGenes)
      ) |>
      dplyr::select(
        ID,
        Description,
        GeneRatio,
        `p.adjust`,
        geneID,
        Cluster,
        Count,
        Up_Vs_Down
      )
  } else {
    InnerDF <- InnerDF |>
      dplyr::filter(pvalue < 0.05, Count >= 5, p.adjust < P.adj) |>
      dplyr::select(
        ID,
        Description,
        GeneRatio,
        `p.adjust`,
        geneID,
        Cluster,
        Count
      ) # Need Fix
  }
  if (nrow(InnerDF) < 2) {
    cli::cli_abort(
      "No useable enriched result, because each term candidated genes is too low. E.g. only 3 genes take part in a pathway."
    )
  }
  InnerDF_raw <- InnerDF
  .checkNrows(InnerDF, force = force)
  obj <- InnerDF |>
    dplyr::mutate(
      PCT = sapply(InnerDF$GeneRatio, function(x) eval(parse(text = x))) * 100
    ) |>
    dplyr::mutate(Padj = signif(p.adjust, 2), PCT = signif(PCT, 2))
  if (sum(colnames(InnerDF) == "Up_Vs_Down") > 0) {
    obj <- obj |>
      dplyr::select(
        Description,
        ID,
        Count,
        Cluster,
        PCT,
        Padj,
        geneID,
        Up_Vs_Down
      ) |>
      dplyr::mutate(geneID = gsub("/", ", ", geneID)) |>
      dplyr::rename(`up/dn` = Up_Vs_Down)
  } else {
    obj <- obj |>
      dplyr::select(Description, ID, Count, Cluster, PCT, Padj, geneID) |>
      dplyr::mutate(geneID = gsub("/", ", ", geneID))
  }
  obj <- obj |>
    dplyr::group_by(Cluster) |>
    dplyr::arrange(Padj) |>
    dplyr::slice_head(n = nTop) |>
    dplyr::ungroup()
  obj0 <- obj |> gt_ora(ClusterNum = ClusterNum, objname = objname, ...)
  obj2 <- obj
  obj3 <- obj2 |> genMetaGM(type = "ORA")
  obj3_1 <- obj3[[1]] |> remove_more_updownInfo()
  obj3_2 <- obj3[[2]]
  objA <- new.egt(
    obj2,
    obj0,
    obj3_1,
    obj3_2,
    clsObj[[3]],
    clsObj[[2]],
    InnerDF_000,
    list(ClusterNum = ClusterNum, 
      P.adj = P.adj,
      force = force,
      objname = "CleanedObj",
      nTop = nTop,
      method = method)
  )
  return(objA)
}

remove_more_updownInfo <- function(x) {
  x <- lapply(x, function(x) {
    x <- x[!grepl("Up_Genes", x)]
    x <- x[!grepl("Down_Genes", x)]
    x <- x[x != ""]
    return(x)
  })
}

.genGSEAGT <- function(
  x,
  ClusterNum,
  P.adj = 0.05,
  force = F,
  objname,
  nTop,
  method,
  ...
) {
  InnerDF <- x
  InnerDF -> InnerDF_000
  InnerDF <- InnerDF |> dplyr::filter(pvalue < 0.05)
  if (dim(x)[1] <= 5) {
    cli::cli_abort(
      "less than 5 enrichment result satisfied that the origin P<0.05(not adjusted)"
    )
  }
  .checkRowNames(x, "GSEA")
  ClusterNum0 <- ClusterNum
  ClusterNum <- .genClusterNum(
    x = x,
    ClusterNum = ClusterNum0,
    force = force
  ) |>
    round()
  clsObj <- .enrichpws(InnerDF$ID, InnerDF$core_enrichment, ClusterNum, method)
  suppressMessages({
    InnerDF <- InnerDF |> dplyr::left_join(clsObj[[1]])
  }) # Merge according to "ID"
  InnerDF <- InnerDF |>
    dplyr::filter(pvalue < 0.05, abs(NES) >= 0.9, p.adjust < P.adj) |>
    dplyr::select(ID, Description, NES, `p.adjust`, Cluster, core_enrichment) # Need Fix
  if (nrow(InnerDF) < 2) {
    cli::cli_abort(
      "No useable enriched result, because each term candidated genes is too low. E.g. only 3 genes take part in a pathway."
    )
  }
  InnerDF_raw <- InnerDF
  .checkNrows(InnerDF, force = force)
  obj <- InnerDF |>
    dplyr::mutate(absNES = abs(NES)) |>
    dplyr::mutate(Reg = ifelse(NES > 0, "red", "forestgreen")) |>
    dplyr::mutate(Padj = signif(p.adjust, 2), absNES = signif(absNES, 4)) |>
    dplyr::select(
      Description,
      ID,
      Reg,
      absNES,
      Cluster,
      Padj,
      core_enrichment
    ) |>
    dplyr::mutate(core_enrichment = gsub("/", ", ", core_enrichment)) #|> need fix colors
  obj <- obj |>
    dplyr::group_by(Cluster) |>
    dplyr::arrange(Padj) |>
    dplyr::slice_head(n = nTop) |>
    dplyr::ungroup()
  obj0 <- obj |> gt_gsea(ClusterNum = ClusterNum, objname = objname, ...)
  obj2 <- obj |> dplyr::mutate(Reg = ifelse(Reg == "red", "UpReg", "DownReg"))
  obj3 <- obj2 |> genMetaGM(type = "GSEA")
  obj3_1 <- obj3[[1]]
  obj3_2 <- obj3[[2]]
  objA <- new.egt(
    obj2,
    obj0,
    obj3_1,
    obj3_2,
    clsObj[[3]],
    clsObj[[2]],
    InnerDF_000,
    list(ClusterNum = ClusterNum, 
      P.adj = P.adj,
      force = force,
      objname = "CleanedObj",
      nTop = nTop,
      method = method)
  )
  return(objA)
}
# attachment::att_amend_desc()