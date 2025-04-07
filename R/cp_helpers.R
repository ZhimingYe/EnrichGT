#' Create a ranked gene list for GSEA analysis
#'
#' Takes gene identifiers and corresponding weights (like log2 fold changes) and returns
#' a ranked vector suitable for Gene Set Enrichment Analysis (GSEA).
#'
#' @importFrom fgsea fgsea
#' @importFrom utils stack
#'
#' @param genes Character vector of gene identifiers (e.g., gene symbols or ENTREZ IDs)
#' @param weights Numeric vector of weights for each gene (typically log2 fold changes)
#'
#' @return A named numeric vector sorted in descending order by weight, where:
#'   - Names are gene identifiers
#'   - Values are the corresponding weights
#'
#' @export
#' @author Zhiming Ye
#'
#' @examples
#' \dontrun{
#' # Example using differential expression results
#' genes <- c("TP53", "BRCA1", "EGFR")
#' log2fc <- c(1.5, -2.1, 0.8)
#' ranked_genes <- genes_with_weights(genes, log2fc)
#' }
genes_with_weights <- function(genes, weights) {
  genes -> Gene
  weights -> Weight
  log2FC <- Weight
  Genetable <- data.frame(Gene = Gene, log2FC = log2FC)
  Genetable <- Genetable |> dplyr::arrange(desc(log2FC))
  GSElist <- as.numeric(Genetable$log2FC)
  names(GSElist) <- Genetable$Gene
  GSElist = sort(GSElist, decreasing = TRUE)
  return(GSElist)
}


#' Parse GMT format gene set files
#'
#' Reads gene set files in GMT format (e.g., from MSigDB or WikiPathways) and converts
#' them to a data frame suitable for enrichment analysis. Can optionally convert ENTREZ IDs
#' to gene symbols.
#'
#' @param gmtfile Path to GMT format file
#' @param OrgDB Annotation database for ID conversion (e.g., org.Hs.eg.db for human).
#'   Required if convert_2_symbols=TRUE.
#' @param convert_2_symbols Logical indicating whether to convert ENTREZ IDs to gene symbols.
#'   Default is TRUE.
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item term: Gene set name
#'     \item gene: Gene identifiers (symbols or ENTREZ IDs)
#'   }
#'   If input has 3 columns, includes an additional ID column.
#'
#' @importFrom utils stack
#' @importFrom tibble as_tibble
#' @export
#' @author Original GMT parser by Guangchuang Yu (https://github.com/YuLab-SMU/gson).
#'   Cache system and enhancements by Zhiming Ye.
#'
#' @examples
#' \dontrun{
#' # Read MSigDB hallmark gene sets
#' gmt_file <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "EnrichGT")
#' gene_sets <- database_from_gmt(gmt_file)
#'
#' # Read WikiPathways with ENTREZ to symbol conversion
#' gmt_file <- "wikipathways-20220310-gmt-Homo_sapiens.gmt"
#' gene_sets <- database_from_gmt(gmt_file, OrgDB = org.Hs.eg.db)
#' }

database_from_gmt <- function(gmtfile, OrgDB = NULL, convert_2_symbols = T) {
  x <- UniversalInternalDBFetcher("SelfBuild", NULL, gmtfile)
  if (convert_2_symbols) {
    if (is.null(OrgDB)) {
      cli::cli_abort("No valid OrgDB, please load it and input")
    }
    x <- x |> as_tibble()
    if (ncol(x) > 3) {
      cli::cli_abort("Not valid gmt file! ")
    }
    if (nrow(x) < 20) {
      Target <- nrow(x)
    } else {
      Target <- 20
    }
    if (
      ((is_numeric_string(x[1:Target, ncol(x), drop = T]) |> sum()) / Target) >
        0.75
    ) {
      cli::cli_alert_info(
        "Detectd numeric gene ids, try to convert to symols.\n You can set convert_2_symbols=F to disable this."
      )
      colnames(x)[ncol(x)] <- "ENTREZID"
      tryCatch(
        {
          qq <- convert_annotations_genes(
            x$ENTREZID,
            from_what = "ENTREZID",
            to_what = "SYMBOL",
            OrgDB = OrgDB
          )
          if (nrow(qq |> na.omit()) < 5) {
            cli::cli_abort(
              "Error in convert. Please set convert_2_symbols=F to disable this"
            )
          }

          x <- x |> left_join(qq, by = "ENTREZID")
          x <- x |> dplyr::select(-ENTREZID)
        },
        error = function(e) {
          cli::cli_abort(
            "Error in convert. Please set convert_2_symbols=F to disable this"
          )
        }
      )
    }
  }

  return(x)
}


#' Perform Over-Representation Analysis (ORA)
#'
#' Identifies enriched biological pathways or gene sets in a gene list using
#' high-performance C++ implementation with parallel processing support.
#'
#' @description
#' ORA compares the proportion of genes in your target list that belong to specific
#' categories (pathways, GO terms etc.) against the expected proportion in a background
#' set. This implementation uses hash tables for efficient gene counting and supports
#' parallel processing for analyzing multiple gene lists.
#'
#' @param genes Input genes, either:
#'   \itemize{
#'     \item Character vector of gene IDs (e.g., `c("TP53","BRCA1")`)
#'     \item Named numeric vector from `genes_with_weights()` (will split by expression direction)
#'     \item List of gene vectors for multiple comparisons (e.g., by cell type)
#'   }
#' @param database Gene set database, either:
#'   \itemize{
#'     \item Built-in database from `database_GO_BP()`, `database_KEGG()` etc.
#'     \item Custom data frame with columns: (ID, Pathway_Name, Genes) or (Pathway_Name, Genes)
#'     \item GMT file loaded via `database_from_gmt()`
#'   }
#' @param p_adj_method Multiple testing correction method (default "BH" for Benjamini-Hochberg)
#' @param p_val_cut_off Adjusted p-value cutoff (default 0.5)
#' @param background_genes Custom background genes (default: all genes in database)
#' @param min_geneset_size Minimum genes per set (default 10)
#' @param max_geneset_size Maximum genes per set (default 500)
#' @param multi_cores (Please don't use this since it has several known bugs) Number of cores for parallel processing (default 0 = serial)
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item ID: Gene set identifier
#'     \item Description: Gene set name
#'     \item GeneRatio: Enriched genes / input genes
#'     \item BgRatio: Set genes / background genes
#'     \item pvalue: Raw p-value
#'     \item p.adjust: Adjusted p-value
#'     \item geneID: Enriched genes
#'     \item Count: Number of enriched genes
#'   }
#'   For weighted input, additional columns show up/down regulated genes.
#'
#' @export
#' @author Zhiming Ye
#'
#' @examples
#' \dontrun{
#' # Basic ORA with GO Biological Processes
#' genes <- c("TP53", "BRCA1", "EGFR", "CDK2")
#' res <- egt_enrichment_analysis(genes, database_GO_BP())
#'
#' # ORA with DEG results (split by direction)
#' deg_genes <- genes_with_weights(DEG$gene, DEG$log2FC)
#' res <- egt_enrichment_analysis(deg_genes, database_KEGG())
#'
#' # Multi-group ORA with parallel processing
#' gene_lists <- list(
#'   Macrophages = c("CD68", "CD163", "CD169"),
#'   Fibroblasts = c("COL1A1", "COL1A2", "ACTA2")
#' )
#' res <- egt_enrichment_analysis(gene_lists, database_Reactome(), multi_cores=0)
#' }
egt_enrichment_analysis <- function(
  genes,
  database,
  p_adj_methods = "BH",
  p_val_cut_off = 0.5,
  background_genes = NULL,
  min_geneset_size = 10,
  max_geneset_size = 500,
  multi_cores = 0
) {
  if (is.character(genes) | is.numeric(genes)) {
    result <- doEnrich_Internal(
      genes,
      database,
      p_adj_methods,
      p_val_cut_off,
      background_genes,
      min_geneset_size,
      max_geneset_size
    )
  } else if (is.list(genes) & multi_cores <= 1) {
    result <- lapply(genes, function(x) {
      tryCatch(
        {
          res <- doEnrich_Internal(
            genes = x,
            database,
            p_adj_methods,
            p_val_cut_off,
            background_genes,
            min_geneset_size,
            max_geneset_size
          )
          return(res)
        },
        error = function(e) {
          return(data.frame(ERROR = "error..."))
        }
      )
    })
  } else if (is.list(genes) & multi_cores >= 2) {
    cli::cli_alert_warning(
      "Multi-cores is no longer support. Please set multi_cores to <= 1."
    )
  }
  tryCatch(
    {
      if (is.data.frame(result)) {
        result <- result |>
          dplyr::filter(pvalue < p_val_cut_off) |>
          dplyr::arrange(pvalue) |>
          tibble::as_tibble()
      } else if (is.list(result)) {
        result <- lapply(result, function(x) {
          if (is.data.frame(x)) {
            if (ncol(x) > 3 & nrow(x) > 3) {
              z <- x |>
                dplyr::filter(pvalue < p_val_cut_off) |>
                dplyr::arrange(pvalue) |>
                tibble::as_tibble()
              return(z)
            }
          }
        })
      }
    },
    error = function(e) {
      message_wrong_db()
      cli::cli_abort("No useable result!")
    }
  )
  return(result)
}

#' Perform Gene Set Enrichment Analysis (GSEA)
#'
#' Identifies enriched biological pathways in a ranked gene list using the fgsea algorithm.
#'
#' @description
#' GSEA analyzes whether predefined gene sets show statistically significant enrichment
#' at the top or bottom of a ranked gene list. This implementation uses the fast fgsea
#' algorithm from the fgsea package.
#'
#' @param genes A named numeric vector of ranked genes where:
#'   \itemize{
#'     \item Names are gene identifiers
#'     \item Values are ranking metric (e.g., log2 fold change, PCA loading)
#'   }
#'   Must be sorted in descending order (use `genes_with_weights()` to prepare)
#' @param database Gene set database, either:
#'   \itemize{
#'     \item Built-in database from `database_GO_BP()`, `database_KEGG()` etc.
#'     \item Custom data frame with columns: (ID, Pathway_Name, Genes) or (Pathway_Name, Genes)
#'     \item GMT file loaded via `database_from_gmt()`
#'   }
#' @param p_adj_method Multiple testing correction method (default "BH" for Benjamini-Hochberg)
#' @param p_val_cut_off Adjusted p-value cutoff (default 0.5)
#' @param min_geneset_size Minimum genes per set (default 10)
#' @param max_geneset_size Maximum genes per set (default 500)
#' @param gseaParam GSEA parameter controlling weight of ranking (default 1)
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item ID: Gene set identifier
#'     \item Description: Gene set name
#'     \item ES: Enrichment score
#'     \item NES: Normalized enrichment score
#'     \item pvalue: Raw p-value
#'     \item p.adjust: Adjusted p-value
#'     \item core_enrichment: Leading edge genes
#'   }
#'
#' @export
#' @author Based on fgsea package by Alexey Sergushichev et al.
#' @importFrom fgsea fgsea
#'
#' @examples
#' \dontrun{
#' # Using differential expression results
#' ranked_genes <- genes_with_weights(DEG$gene, DEG$log2FC)
#' res <- egt_gsea_analysis(ranked_genes, database_GO_BP())
#'
#' # Using PCA loadings
#' ranked_genes <- genes_with_weights(rownames(pca$rotation), pca$rotation[,1])
#' res <- egt_gsea_analysis(ranked_genes, database_KEGG())
#'
#' # Custom gene sets from GMT file
#' ranked_genes <- genes_with_weights(genes, weights)
#' res <- egt_gsea_analysis(ranked_genes, database_from_gmt("pathways.gmt"))
#' }
#'
egt_gsea_analysis <- function(
  genes,
  database,
  p_val_cut_off = 0.5,
  min_geneset_size = 10,
  max_geneset_size = 500,
  gseaParam = 1
) {
  if (!is.list(genes)) {
    res <- egt_gsea_analysis_internal(
      genes = genes,
      database = database,
      p_val_cut_off = p_val_cut_off,
      min_geneset_size = min_geneset_size,
      max_geneset_size = max_geneset_size,
      gseaParam = gseaParam
    )
  } else {
    res <- lapply(genes, function(x) {
      tryCatch(
        {
          output <- egt_gsea_analysis_internal(
            genes = x,
            database = database,
            p_val_cut_off = p_val_cut_off,
            min_geneset_size = min_geneset_size,
            max_geneset_size = max_geneset_size,
            gseaParam = gseaParam
          )
        },
        error = function(e) {
          output <- NULL
        }
      )
      return(output)
    })
  }
  res <- tibble::as_tibble(res) # better printing
  return(res)
}

egt_gsea_analysis_internal <- function(
  genes,
  database,
  p_val_cut_off = 0.5,
  min_geneset_size = 10,
  max_geneset_size = 500,
  gseaParam = 1,
  for_figures = FALSE
) {
  suppressPackageStartupMessages({
    requireNamespace("fgsea")
  })
  tryCatch(
    {
      if (ncol(database) != 2 & ncol(database) != 3) {
        message_wrong_db()
        cli::cli_abort("Not valid database")
      }
    },
    error = function(e) {
      message_wrong_db()
      cli::cli_abort("Not valid database")
    }
  )
  if (ncol(database) == 3) {
    colnames(database) <- c("ID", "term", "gene")
    db0 <- database[, c(1, 2)]
    database <- database[, c(2, 3)]
  } else {
    colnames(database)[1] <- "term"
    db0 <- data.frame(ID = database$term, term = database$term)
  }
  db0 <- db0 |>
    dplyr::mutate(CheckDup = paste0(ID, term)) |>
    dplyr::filter(!duplicated(CheckDup)) |>
    dplyr::select(-CheckDup) |>
    dplyr::rename(pathway = term) # Because of output is pathway
  colnames(database) <- c("term", "gene")
  database2 <- split(database$gene, database$term)
  t1 <- Sys.time()
  if (for_figures) {
    return(list(PWS = database2, RKS = genes))
  }
  fgseaRes <- fgsea::fgsea(
    pathways = database2,
    stats = genes,
    minSize = min_geneset_size,
    maxSize = max_geneset_size,
    gseaParam = gseaParam
  )
  fgseaRes <- fgseaRes |> dplyr::left_join(db0, by = "pathway")
  t2 <- Sys.time()
  cli::cli_alert_success(paste0(
    "Sucessful GSEA, time last ",
    (t2 - t1),
    " secs."
  ))
  fgseaRes2 <- data.frame(
    ID = fgseaRes$ID,
    Description = fgseaRes$pathway,
    ES = fgseaRes$ES,
    NES = fgseaRes$NES,
    pvalue = fgseaRes$pval,
    p.adjust = fgseaRes$padj,
    core_enrichment = sapply(
      fgseaRes$leadingEdge,
      function(x) paste(x, collapse = "/")
    )
  ) #,fgsea_leadingEdge = fgseaRes$leadingEdge
  tryCatch(
    {
      fgseaRes2 <- fgseaRes2 |>
        dplyr::filter(pvalue < p_val_cut_off) |>
        dplyr::arrange(desc(NES))
    },
    error = function(e) {
      message_wrong_db()
      cli::cli_abort("No useable result!")
    }
  )
  return(fgseaRes2)
}

message_wrong_db <- function() {
  cli::cli_alert_danger(
    "Please re-check the database OrgDB argument or input GMT file matchs the inputed genes. \nYou should input gene symbols when your database contains gene symbols, ENSEMBL IDs when ENSEMBL IDs. \nAnother possible reason is that you may forget to run: "
  )
  cli::cli_code(
    "library(org.Hs.eg.db) # for humans, if mouse use org.Mm.eg.db. Please re-check"
  )
  cli::cli_code(
    "See https://zhimingye.github.io/EnrichGT/database.html for further details"
  )
}
