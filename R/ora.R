parse_gs <- function(gset) {
  colnames(gset) <- c("term", "gene")
  gset2 <- split(gset, gset$term)
  gset2 <- lapply(gset2, function(x) x$gene)
  gset2 -> gene_sets
  return(list(background_genes = 0, gene_sets = gene_sets))
}

#' @importFrom qvalue qvalue
#' @importFrom stats p.adjust
#' @importFrom utils stack
#' @importFrom dplyr filter arrange left_join select mutate case_when
#' @useDynLib EnrichGT
#' @importFrom Rcpp sourceCpp
doEnrich_Internal <- function(
  genes,
  database,
  p_adj_methods,
  p_val_cut_off,
  background_genes,
  min_geneset_size,
  max_geneset_size
) {
  has_direction <- F
  if (!is.null(names(genes))) {
    genes2 <- genes
    genes <- names(genes)
    genes_up <- names(genes2[genes2 > 0])
    genes_down <- names(genes2[genes2 < 0])
    has_direction <- T
  }

  dblist <- prepare_database(database, "TERMs")
  dblist$db0 -> db0
  dblist$database -> database

  termCount <- table(database$term) |> as.data.frame()
  colnames(termCount)[1] <- "TERMs"
  termCount <- termCount |>
    dplyr::filter(Freq >= min_geneset_size, Freq <= max_geneset_size)
  if (nrow(termCount) == 0) {
    cli::cli_abort("No gene set match the asked gene set size ")
  }
  termCount$TERMs <- as.character(termCount$TERMs)
  background_genes_ <- database$gene |> table() |> names()
  database0 <- database |> dplyr::filter(term %in% termCount$TERMs)

  tgtGs <- parse_gs(database0)
  if (!is.null(background_genes)) {
    tgtGs[[1]] <- background_genes
  } else {
    tgtGs[[1]] <- background_genes_
  }
  genes <- genes[genes %in% tgtGs[[1]]]

  HittedCheck <- data.frame(gene = genes, hitted = "Hitted")
  database2 <- database0 |> dplyr::left_join(HittedCheck, by = "gene")
  database2 <- database2[!is.na(database2$hitted), ]
  database2 <- database2[database2$hitted == "Hitted", ]
  if (!has_direction) {
    hitted_result <- database2 |>
      dplyr::group_by(term) |>
      dplyr::summarise(gene_list = paste(gene, collapse = "/")) |>
      dplyr::ungroup()
  } else {
    hitted_result <- database2 |>
      dplyr::group_by(term) |>
      dplyr::summarise(
        gene_list = paste(gene, collapse = "/"),
        UpGenes = paste(gene[gene %in% genes_up], collapse = "/"),
        DownGenes = paste(gene[gene %in% genes_down], collapse = "/"),
        UpSummary = sum(gene %in% genes_up),
        DownSummary = sum(gene %in% genes_down)
      ) |>
      dplyr::ungroup()
  }

  colnames(hitted_result)[1] <- "TERMs"
  hitted_result$TERMs <- as.character(hitted_result$TERMs)

  t1 <- Sys.time()
  if (length(genes) == 0) {
    cli::cli_abort("No gene overlap! Please recheck. ")
  }
  result <- ora_analysis(tgtGs[[2]], genes, tgtGs[[1]])
  result$TERMs <- tgtGs[[2]] |> names()
  df00 <- hitted_result |>
    dplyr::left_join(result, by = "TERMs") |>
    dplyr::left_join(termCount, by = "TERMs") |>
    dplyr::left_join(db0, by = "TERMs")
  df00$padj <- p.adjust(df00$PValue, method = p_adj_methods)
  qobj <- tryCatch(
    qvalue(p = df00$padj, lambda = 0.05, pi0.method = "bootstrap"),
    error = function(e) NULL
  )
  if (inherits(qobj, "qvalue")) {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }
  df00$qvalues <- qvalues
  df00$LengthOfInput <- length(genes)
  df00$LengthOfBG <- length(tgtGs[[1]])
  if (!has_direction) {
    res <- data.frame(
      ID = df00$ID,
      Description = df00$TERMs,
      GeneRatio = paste0(df00$Overlap, "/", df00$LengthOfInput),
      BgRatio = paste0(df00$Freq, "/", df00$LengthOfBG),
      pvalue = df00$PValue,
      p.adjust = df00$padj,
      qvalue = df00$qvalues,
      geneID = df00$gene_list,
      Count = df00$Overlap
    )
  } else {
    res <- data.frame(
      ID = df00$ID,
      Description = df00$TERMs,
      GeneRatio = paste0(df00$Overlap, "/", df00$LengthOfInput),
      BgRatio = paste0(df00$Freq, "/", df00$LengthOfBG),
      pvalue = df00$PValue,
      p.adjust = df00$padj,
      qvalue = df00$qvalues,
      geneID = df00$gene_list,
      Count = df00$Overlap,
      UpGenes = df00$UpGenes,
      DownGenes = df00$DownGenes,
      UpSummary = df00$UpSummary,
      DownSummary = df00$DownSummary,
      Up_Vs_Down = paste0(df00$UpSummary, "/", df00$DownSummary)
    )
  }
  t2 <- Sys.time()
  timeLast <- t2 - t1
  cli::cli_alert_success(paste0("Done ORA in ", timeLast, " sec."))
  attr(res, "Package") <- "EnrichGT"
  attr(res, "Input") <- genes
  attr(res, "Database") <- database
  attr(res, "Other_Params") <- list(p_val_cut_off,
                                    background_genes,
                                    min_geneset_size,
                                    max_geneset_size,
                                    p_adj_methods)
  attr(res, "Time") <- Sys.time()
  return(res)
}
