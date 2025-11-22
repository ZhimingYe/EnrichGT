#' Visualize enrichment results using simple plot
#'
#' @description
#' This plot is the most widely like `enrichplot::dotplot()`used method to visualize enriched terms. It shows the enrichment scores (e.g. p values) and gene ratio or NES as dot size and color / or bar height. Users can specify the number of terms using `ntop` or selected terms to color via the `low.col` and `hi.col`.
#' @param x a data frame form enriched result like `egt_enrichment_analysis()` or `egt_gsea_analysis()`, or an re-clustered `EnrichGT` object
#' @param ntop Show top N in each cluster. In default, for origin enriched result, showing top 15; for re-clustered object, showing top 5 in each cluster.
#' @param showIDs bool, show pathway IDs or not. Default is FALSE
#' @param low.col the color for the lowest
#' @param hi.col the color for the highest
#' @param max_len_descript the label format length, default as 40.
#' @param keepAll Do filtering to avoid overlap of same genes or not
#' @param maskNoise (Only works with re-enriched object) Cut-off value to mask rare population in cluster tree. Less than its value in a specifc child tree will be ignore (because it may hit only by coincidence). Set maskNoise = 0 to ignore this.
#' @param P.adj (Only works with origin data.frame) If pass an origin data.frame from original enriched result, you can specify the P-adjust value cut off. If is null, default is 0.05. When passing `EnrichGT_obj`, this filter is previously done by `egt_recluster_analysis`.
#' @param ... Other param
#'
#' @returns a ggplot2 object
#' @export
#' @importFrom ggplot2 fortify
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_gradient
#' @importFrom ggplot2 scale_color_continuous
#' @importFrom ggplot2 scale_fill_continuous
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 guide_colorbar
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_size
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 geom_col
#' @importFrom methods is
#' @importFrom forcats fct_reorder
#' @author Zhiming Ye
egt_plot_results <- function(
  x,
  ntop = NULL,
  showIDs = F,
  max_len_descript = 40,
  keepAll = F,
  maskNoise = 3,
  ...,
  P.adj = NULL
) {
  if (is.null(ntop)) {
    ntop <- 15
  }
  if (class(x)[1] == "EnrichGT_obj") {
    if (is.null(ntop)) {
      ntop <- 5
    }
    if (x@fused) {
      cli::cli_alert_info(
        "Found a fused EnrichGT_obj, you can set showIDs=T to show the source database the pathways from. "
      )
    }
    tryCatch(
      {
        if (maskNoise > ntop) {
          stop("`maskNoise` should not be greater than `ntop`")
        }
        if (maskNoise > 0) {
          ttt0 <- table(x@enriched_result$Cluster)
          x@enriched_result <- x@enriched_result |>
            dplyr::filter(Cluster %in% names(ttt0)[ttt0 > maskNoise])
        }
      },
      error = function(e) {
        cli::cli_abort("Not an EnrichGT_obj. or `maskNoise` param is too big")
      }
    )

    if (sum(colnames(x@enriched_result) == "absNES") == 0) {
      figure0 <- ORA2dp(
        x,
        ntop = ntop,
        showIDs = showIDs,
        max_len_descript = max_len_descript,
        ...
      )
    } else {
      figure0 <- GSEA2dp(
        x,
        ntop = ntop,
        showIDs = showIDs,
        max_len_descript = max_len_descript,
        ...
      )
    }
  } else if (is.data.frame(x)) {
    plotingTemp <- new.env()
    tryCatch(
      {
        assign("PadjVal", P.adj, envir = plotingTemp)
        if (is.null(P.adj)) {
          stop()
        }
      },
      error = function(e) {
        assign("PadjVal", 0.05, envir = plotingTemp)
        cli::cli_alert_info(
          "Use Default P-adjust cut-off 0.05. You can pass `P.adj=xxx` arugument to filter. "
        )
      }
    )
    if (sum(colnames(x) == "NES") != 0) {
      InnerDF <- x |>
        dplyr::filter(p.adjust < (plotingTemp$PadjVal)) |>
        dplyr::select(ID, Description, NES, `p.adjust`, core_enrichment) # Need Fix
      obj <- InnerDF |>
        dplyr::mutate(absNES = abs(NES)) |>
        dplyr::mutate(Reg = ifelse(NES > 0, "red", "forestgreen")) |>
        dplyr::mutate(Padj = signif(p.adjust, 2), absNES = signif(absNES, 4)) |>
        dplyr::select(Description, ID, Reg, absNES, Padj, core_enrichment) |>
        dplyr::mutate(core_enrichment = gsub("/", ", ", core_enrichment)) |>
        dplyr::mutate(Cluster = "GSEA Results")
      if (!keepAll) {
        obj <- obj |>
          dplyr::group_by(core_enrichment) |>
          dplyr::arrange(desc(absNES), .by_group = TRUE) |>
          dplyr::slice_head(n = 2) |>
          dplyr::ungroup()
      }
      obj <- obj |>
        dplyr::group_by(Reg) |>
        dplyr::arrange(desc(absNES), .by_group = TRUE) |>
        dplyr::slice_head(n = round(round(ntop / 2) + 1 + round(ntop / 20))) |> # For balance
        dplyr::ungroup() # Balance up and down
      figure0 <- GSEA2dp(
        obj,
        ntop = ntop,
        showIDs = showIDs,
        max_len_descript = max_len_descript,
        ...
      )
    } else if (sum(colnames(x) %in% c("Up_Vs_Down", "up_dn")) == 0) {
      InnerDF <- x |>
        dplyr::filter(p.adjust < (plotingTemp$PadjVal)) |>
        dplyr::select(ID, Description, GeneRatio, `p.adjust`, geneID, Count) # Need Fix
      obj <- InnerDF |>
        dplyr::mutate(
          PCT = sapply(InnerDF$GeneRatio, function(x) eval(parse(text = x))) *
            100
        ) |>
        dplyr::mutate(Padj = signif(p.adjust, 2), PCT = signif(PCT, 2)) |>
        dplyr::select(Description, ID, Count, PCT, Padj, geneID) |>
        dplyr::mutate(geneID = gsub("/", ", ", geneID)) |>
        dplyr::mutate(Cluster = "ORA Results")
      if (!keepAll) {
        obj <- obj |>
          dplyr::group_by(geneID) |>
          dplyr::arrange(Padj, .by_group = TRUE) |>
          dplyr::slice_head(n = 2) |>
          dplyr::ungroup()
      }
      figure0 <- ORA2dp(
        obj,
        ntop = ntop,
        showIDs = showIDs,
        max_len_descript = max_len_descript,
        ...
      )
    } else {
      if (sum(colnames(x) == "up_dn") > 0) {
        x <- x |> dplyr::rename(Up_Vs_Down = `up_dn`)
      }
      InnerDF <- x |>
        dplyr::filter(p.adjust < (plotingTemp$PadjVal)) |>
        dplyr::select(
          ID,
          Description,
          GeneRatio,
          `p.adjust`,
          geneID,
          Count,
          Up_Vs_Down
        ) # Need Fix
      obj <- InnerDF |>
        dplyr::mutate(
          PCT = sapply(InnerDF$GeneRatio, function(x) eval(parse(text = x))) *
            100
        ) |>
        dplyr::mutate(Padj = signif(p.adjust, 2), PCT = signif(PCT, 2)) |>
        dplyr::select(Description, ID, Count, PCT, Padj, geneID, Up_Vs_Down) |>
        dplyr::mutate(geneID = gsub("/", ", ", geneID)) |>
        dplyr::mutate(Cluster = "ORA Results")
      if (!keepAll) {
        obj <- obj |>
          dplyr::group_by(geneID) |>
          dplyr::arrange(Padj, .by_group = TRUE) |>
          dplyr::slice_head(n = 2) |>
          dplyr::ungroup()
      }
      figure0 <- ORA2dp(
        obj,
        ntop = ntop,
        showIDs = showIDs,
        max_len_descript = max_len_descript,
        ...
      )
    }
  }
  return(figure0)
}

shorten_labels_words <- function(label, max_length = 40) {
  shorten_one <- function(l) {
    words <- base::strsplit(l, " ")[[1]]
    cum <- base::cumsum(nchar(words) + 1)
    if (max(cum) <= max_length) {
      return(l)
    }
    cut <- max(which(cum <= max_length))
    paste(paste(words[1:cut], collapse = " "), "...")
  }
  out <- vapply(label, shorten_one, character(1))
  idx <- stats::ave(seq_along(out), out, FUN = seq_along)
  out[idx > 1] <- paste0(out[idx > 1], "(", idx[idx > 1], ")")
  out
}

ORA2dp <- function(
  x,
  ntop = 7,
  showIDs = F,
  low.col = "#ff6f81",
  hi.col = "#78cfe5",
  max_len_descript = 40,
  ...
) {
  if (is.list(x) & !is.data.frame(x)) {
    cli::cli_abort(
      "For a list object, please run plotting for every object inside list, instead of the whole list."
    )
  }
  TempPlotingEnv <- new.env()
  tryCatch(
    {
      if (
        dim(x@enriched_result)[1] < 2 |
          sum(colnames(x@enriched_result) == "Count") == 0
      ) {
        cli::cli_abort("ERROR! ")
      } else {
        if (sum(colnames(x@enriched_result) == "up_dn") > 0) {
          kk <- x@enriched_result |> dplyr::rename(Up_Vs_Down = `up_dn`)
        } else {
          kk <- x@enriched_result
        }
        assign("df0", kk, envir = TempPlotingEnv)
      }
    },
    error = function(e) {
      assign("df0", x, envir = TempPlotingEnv)
      cli::cli_alert_warning(
        "You are drawing origin results, for better result you can re-cluster it by egt_recluster_analysis()"
      )
    }
  )
  if (
    showIDs &
      (sum(TempPlotingEnv$df0$ID == TempPlotingEnv$df0$Description) <= 5)
  ) {
    TempPlotingEnv$df0$ID <- substr(TempPlotingEnv$df0$ID, 1, 15)
    TempPlotingEnv$df0$Description <- paste0(
      TempPlotingEnv$df0$ID,
      ":",
      TempPlotingEnv$df0$Description
    )
    max_len_descript <- max_len_descript + 16
  }
  tryCatch(
    {
      df <- TempPlotingEnv$df0 |>
        dplyr::group_by(Cluster) |>
        dplyr::slice_min(order_by = Padj, n = ntop, with_ties = FALSE) |>
        dplyr::ungroup()
    },
    error = function(e) {
      cli::cli_alert_warning("Subset ERROR! ")
      df <- TempPlotingEnv$df0
    }
  )
  df$Description <- shorten_labels_words(
    df$Description,
    max_length = max_len_descript
  )
  addUpDnRatio <- F
  if (sum(colnames(df) == "Up_Vs_Down") > 0) {
    df$Description <- paste0(df$Description, " ", df$Up_Vs_Down)
    addUpDnRatio <- T
  }
  px <- ggplot(
    df,
    aes(x = PCT, y = fct_reorder(Description, PCT), size = Count, color = Padj)
  ) +
    geom_point() +
    scale_color_continuous(
      low = low.col,
      high = hi.col,
      name = "adjustedP",
      guide = guide_colorbar(reverse = F)
    ) +
    scale_size(range = c(2, 8)) +
    xlab("Gene Ratio(%)") +
    ylab("Gene Sets") +
    facet_grid(Cluster ~ ., scales = "free", space = "free_y") +
    theme_bw()
  if (addUpDnRatio) {
    px <- px + ylab("Gene Sets [UP-reg/Down-reg Genes Nums]")
  }
  return(px)
}

GSEA2dp <- function(
  x,
  ntop = 7,
  showIDs = F,
  low.col = "#ff6f81",
  hi.col = "#78cfe5",
  max_len_descript = 40,
  ...
) {
  if (is.list(x) & !is.data.frame(x)) {
    cli::cli_abort(
      "For a list object, please run plotting for every object inside list, instead of the whole list."
    )
  }
  TempPlotingEnv <- new.env()
  tryCatch(
    {
      if (
        dim(x@enriched_result)[1] < 2 |
          sum(colnames(x@enriched_result) == "absNES") == 0
      ) {
        cli::cli_abort("ERROR! ")
      } else {
        assign("df0", x@enriched_result, envir = TempPlotingEnv)
      }
    },
    error = function(e) {
      assign("df0", x, envir = TempPlotingEnv)
      cli::cli_alert_warning(
        "You are drawing origin results, for better result you can re-cluster it by egt_recluster_analysis()"
      )
    }
  )
  if (
    showIDs &
      (sum(TempPlotingEnv$df0$ID == TempPlotingEnv$df0$Description) >= 5)
  ) {
    TempPlotingEnv$df0$ID <- substr(TempPlotingEnv$df0$ID, 1, 15)
    TempPlotingEnv$df0$Description <- paste0(
      TempPlotingEnv$df0$ID,
      ":",
      TempPlotingEnv$df0$Description
    )
    max_len_descript <- max_len_descript + 16
  }
  tryCatch(
    {
      df <- TempPlotingEnv$df0 |>
        dplyr::group_by(Cluster) |>
        dplyr::slice_max(order_by = absNES, n = ntop, with_ties = FALSE) |>
        dplyr::ungroup()
    },
    error = function(e) {
      cli::cli_alert_warning("Subset ERROR! ")
      df <- TempPlotingEnv$df0
    }
  )
  df$NES <- ifelse(
    (df$Reg == "UpReg" | df$Reg == "red"),
    df$absNES * (1),
    df$absNES * (-1)
  ) # For the different input
  df$Description <- shorten_labels_words(
    df$Description,
    max_length = max_len_descript
  )
  px <- ggplot(
    df,
    aes(x = NES, y = fct_reorder(Description, absNES), fill = Padj)
  ) +
    geom_col() +
    scale_fill_continuous(
      low = low.col,
      high = hi.col,
      name = "adjustedP",
      guide = guide_colorbar(reverse = F)
    ) +
    scale_size(range = c(2, 8)) +
    xlab("Normalize Enrichment Score(NES)") +
    ylab("Gene Sets") +
    facet_grid(Cluster ~ ., scales = "free", space = "free_y") +
    theme_bw()
  return(px)
}


cocol <- function(n, favor = 1, returnColor = F) {
  if (favor == 3) {
    colorSpace <- c(
      '#E41A1C',
      '#377EB8',
      '#4DAF4A',
      '#984EA3',
      '#F29403',
      '#F781BF',
      '#BC9DCC',
      '#A65628',
      '#54B0E4',
      '#222F75',
      '#1B9E77',
      '#B2DF8A',
      '#E3BE00',
      '#FB9A99',
      '#E7298A',
      '#910241',
      '#00CDD1',
      '#A6CEE3',
      '#CE1261',
      '#5E4FA2',
      '#8CA77B',
      '#00441B',
      '#DEDC00',
      '#DCF0B9',
      '#8DD3C7',
      '#999999'
    )
  } else if (favor == 2) {
    colorSpace <- c(
      "#7ca7ae",
      "#a3b3c9",
      "#788ab2",
      "#edbacd",
      "#687050",
      "#b8c0a8",
      "#908088",
      "#e1b19e",
      "#7fc4da",
      "#e8dff4",
      "#b7988f",
      "#c59d17",
      "#92a761",
      "#75aa7a",
      "#efdfbb",
      "#fabb6e",
      "#fc8002",
      "#addb88",
      "#369f2d",
      "#fac7b3",
      "#ee4431",
      "#b9181a",
      "#cedfef",
      "#92c2dd",
      "#4995c6",
      "#1663a9",
      "#bab4d5",
      "#614099",
      "#45aab4",
      "#038db2",
      "#f9637c",
      "#fe7966",
      "#fff4de",
      "#81d0bb",
      "#a5b8f3",
      "#feaac2",
      "#66C2A5",
      "#8DA0CB",
      "#E78AC3",
      "#A6D854",
      "#FFD92F",
      "#E5C494",
      "#B3B3B3"
    )
  } else if (favor == 1) {
    colorSpace <- c(
      '#0ca9ce',
      '#78cfe5',
      '#c6ecf1',
      '#ff6f81',
      '#ff9c8f',
      '#ffc2c0',
      '#d386bf',
      '#cdb1d2',
      '#fae6f0',
      '#eb6fa6',
      '#ff88b5',
      '#00b1a5',
      "#ffa68f",
      "#ffca75",
      "#b8d8c9",
      "#97bc83",
      "#009f93",
      "#448c99",
      "#db888e",
      "#e397a4",
      "#ead0c7",
      "#8f9898",
      "#bfcfcb"
    )
  }
  if (!returnColor) {
    if (n <= length(colorSpace)) {
      colors <- colorSpace[1:n]
    } else {
      colors <- grDevices::colorRampPalette(colorSpace)(n)
    }
  } else {
    colors <- colorSpace
  }
  return(colors)
}
