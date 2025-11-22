gseaScores3 <- function(geneList, geneSet, exponent = 1) {
  geneSet <- intersect(geneSet, names(geneList))
  if (length(geneSet) == 0L) {
    cli::cli_abort(
      "Can't perform leading edge analysis because the length of gene set is zero. "
    )
  }
  idx <- match(geneSet, names(geneList))
  idx <- idx[!is.na(idx)]
  out <- gseaScores2(
    stats = unname(geneList),
    hitIndices = as.integer(idx),
    exponent = exponent
  )
  df <- out$runningES
  df$gene <- names(geneList)
  list(ES = out$ES, runningES = df)
}

leading_edge <- function(observed_info) {
  n_sets <- length(observed_info)
  if (n_sets == 0L) {
    return(list(
      rank = numeric(0),
      tags = character(0),
      list = character(0),
      signal = character(0),
      leading_edge = character(0),
      core_enrichment = list()
    ))
  }

  N <- nrow(observed_info[[1]]$runningES)

  rank_pos <- numeric(n_sets)
  tags_frac <- numeric(n_sets)
  list_frac <- numeric(n_sets)
  set_size <- numeric(n_sets)
  core_enrich_gen <- vector("list", n_sets)

  for (k in seq_len(n_sets)) {
    info <- observed_info[[k]]
    es <- info$ES
    traj <- info$runningES

    if (!("runningScore" %in% names(traj)) || !("position" %in% names(traj))) {
      stop("runningES must contain columns 'runningScore' and 'position'.")
    }

    scores <- traj$runningScore
    hit_mask <- traj$position == 1L
    n_all <- length(scores)

    set_size[k] <- sum(hit_mask)

    if (es >= 0) {
      idx_ext_full <- which.max(scores)
    } else {
      idx_ext_full <- which.min(scores)
    }

    if (es >= 0) {
      rank_pos[k] <- idx_ext_full
    } else {
      rank_pos[k] <- n_all - idx_ext_full + 1L
    }

    scores_hits <- scores[hit_mask]
    genes_hits <- traj$gene[hit_mask]
    n_hits <- length(scores_hits)

    if (n_hits == 0L) {
      tags_frac[k] <- 0
      list_frac[k] <- 0
      core_enrich_gen[[k]] <- character(0)
      next
    }

    if (es >= 0) {
      idx_ext_hits <- which.max(scores_hits)
      tags_frac[k] <- idx_ext_hits / n_hits
      core_enrich_gen[[k]] <- genes_hits[seq_len(idx_ext_hits)]
    } else {
      idx_ext_hits <- which.min(scores_hits)
      tags_frac[k] <- (n_hits - idx_ext_hits + 1L) / n_hits
      core_enrich_gen[[k]] <- genes_hits[seq(from = idx_ext_hits, to = n_hits)]
    }

    if (es >= 0) {
      list_frac[k] <- idx_ext_full / n_all
    } else {
      list_frac[k] <- (n_all - idx_ext_full + 1L) / n_all
    }
  }

  # signal = tags * (1 - list) * N / (N - setSize)
  denom <- (N - set_size)
  signal_raw <- tags_frac * (1 - list_frac)
  factor_term <- ifelse(denom > 0, N / denom, NA_real_)
  signal_raw <- signal_raw * factor_term

  tags_chr <- paste0(round(tags_frac * 100), "%")
  list_chr <- paste0(round(list_frac * 100), "%")
  signal_chr <- paste0(round(signal_raw * 100), "%")

  leading_desc <- paste0(
    "tags=",
    tags_chr,
    ", list=",
    list_chr,
    ", signal=",
    signal_chr
  )

  out <- list(
    rank = rank_pos,
    tags = tags_chr,
    list = list_chr,
    signal = signal_chr,
    leading_edge = leading_desc,
    core_enrichment = core_enrich_gen
  )

  return(out)
}


egt_leading_edge <- function(geneList, geneSets, exponent = 1) {
  observed_info <- lapply(geneSets, function(gs) {
    gseaScores3(
      geneList = geneList,
      geneSet = gs,
      exponent = exponent
    )
  })
  leading_edge(observed_info)
}
