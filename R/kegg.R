keggModuleList <- function(orgkegg) {
  k0 <- glue::glue("https://rest.kegg.jp/list/pathway/{orgkegg}")
  k1 <- glue::glue("https://rest.kegg.jp/link/pathway/{orgkegg}")
  k2 <- glue::glue("https://rest.kegg.jp/conv/{orgkegg}/ncbi-geneid")
  m0 <- glue::glue("https://rest.kegg.jp/list/module")
  m1 <- glue::glue("https://rest.kegg.jp/link/module/{orgkegg}")

  return(c(k0, k1, k2, m0, m1))
}
# y <- read.delim("https://rest.kegg.jp/link/module/mmu",quote = "\t",header = F)

#' Get KEGG database from KEGG website
#' @description
#' KEGG is a commercialized database. So EnrichGT can't pre-cache them locally. You can use this function to fetch KEGG database pathways and modules.
#' @usage database_KEGG(kegg_organism="hsa",OrgDB = org.Hs.eg.db,kegg_modules=F,local_cache=F)
#'
#' @param kegg_organism Determine which species data from KEGG will be fetch. For human, it would be `hsa`(in default); For mouse, it would be `mmu`. If you wants other species, see `database_kegg_show_organism()` for details.
#' @param OrgDB The AnnotationDbi database to convert KEGG gene ID to gene symbols. For human it would be `org.Hs.eg.db`, for mouse it would be `org.Mm.eg.db`. In AnnotationDbi there are many species, please search `AnnotationDbi` for other species annotation database.
#' @param kegg_modules If TRUE, returns KEGG module; If FALSE returns KEGG pathways. In default, this is setted to FALSE to get mouse commonly used KEGG pathways.
#' @param local_cache cache a copy in local working folder. It will be saved as a `.enrichgt_cache` file in working dictionary. The `.enrichgt_cache` is just a `.rds` file, feel free to read it using `readRDS()`.
#'
#' @returns data.frame contains KEGG annotations
#' @export
#'
#' @rdname KEGGhelp
database_KEGG <- function(
  kegg_organism = "hsa",
  OrgDB = org.Hs.eg.db,
  kegg_modules = F,
  local_cache = T
) {
  df <- database_KEGG_internal(
    kegg_organism = kegg_organism,
    OrgDB = OrgDB,
    kegg_modules = kegg_modules,
    local_cache = local_cache
  )
  return(df)
}
database_KEGG_internal <- function(
  kegg_organism = "hsa",
  OrgDB = org.Hs.eg.db,
  kegg_modules = F,
  local_cache = T
) {
  if (!kegg_modules) {
    fn <- paste0("KEGGPathway_", xfun::md5(keggModuleList(kegg_organism)))
  } else {
    fn <- paste0("KEGGModules_", xfun::md5(keggModuleList(kegg_organism)))
  }
  have_read <- F
  if (local_cache) {
    if (file.exists(paste0(fn, ".enrichgt_cache"))) {
      tryCatch(
        {
          cachedFile <- readRDS(paste0(fn, ".enrichgt_cache"))
          have_read <- T
        },
        error = function(e) {
          cli::cli_abort(
            "Load local cache error! Please re-check or set local_cache=F to using online files"
          )
        }
      )

      cli::cli_alert_success("Found on disk cache file and loaded. ")
      assign(fn, cachedFile, envir = db_getter_env)
    }
  }
  if (exists(fn, envir = db_getter_env)) {
    finalDF <- get(fn, envir = db_getter_env)
    cli::cli_alert_info(paste0("Use cached KEGG downloaded files: ", fn))
  } else {
    if (!kegg_modules) {
      targetsind <- 1:3
    } else {
      targetsind <- c(4, 5, 3)
    }
    keggs <- lapply(keggModuleList(kegg_organism)[targetsind], function(x) {
      tryCatch(
        {
          cli::cli_alert_info(glue::glue("downloading {x}..."))
          y <- retry_function(
            read.delim,
            ntry = 3,
            delay = 3,
            x,
            quote = "\t",
            header = F
          )
          Sys.sleep(2)
          return(y)
        },
        error = function(e) {
          cli::cli_abort(
            "Please re-check your kegg_organism input and OrgDB input. Or please check your network(Because of the limitation of KEGG license, you must download them instead of built in in EnrichGT). "
          )
        }
      )
    })
    colnames(keggs[[1]]) <- c("Pws", "Terms")
    colnames(keggs[[2]]) <- c("Kegggenes", "Pws")
    colnames(keggs[[3]]) <- c("Ncbigenes", "Kegggenes")
    if (!kegg_modules) {
      keggs[[2]]$Pws <- gsub("path:", "", keggs[[2]]$Pws)
    } else {
      keggs[[2]]$Pws <- s_(keggs[[2]]$Pws, "_", 2)
    }

    finalDF <- keggs[[2]] |>
      dplyr::left_join(keggs[[1]], by = "Pws", relationship = "many-to-many") |>
      dplyr::left_join(
        keggs[[3]],
        by = "Kegggenes",
        relationship = "many-to-many"
      )

    assign(fn, finalDF, envir = db_getter_env)
  }
  if (!local_cache | (local_cache & !have_read)) {
    finalDF$Ncbigenes <- gsub("ncbi-geneid:", "", finalDF$Ncbigenes)
    cvtDF <- convert_annotations_genes(
      finalDF$Ncbigenes,
      "ENTREZID",
      "SYMBOL",
      OrgDB = OrgDB
    )
    finalDF <- finalDF |> dplyr::rename(ENTREZID = Ncbigenes)
    finalDF <- finalDF |>
      dplyr::left_join(cvtDF, by = "ENTREZID", relationship = "many-to-many") |>
      dplyr::select(Pws, Terms, SYMBOL)
    finalDF$ALL <- paste0(finalDF$Pws, finalDF$Terms, finalDF$SYMBOL)
    finalDF <- finalDF[
      !duplicated(finalDF$ALL),
      -which(colnames(finalDF) == "ALL")
    ]

    finalDF$Terms <- s_(finalDF$Terms, " - ", 1)
  }

  res0 <- finalDF
  if (local_cache & !have_read) {
    saveRDS(res0, (paste0(fn, ".enrichgt_cache")))
    cli::cli_alert_success(paste0(
      "Wrote KEGG cache on disk: ",
      fn,
      ".enrichgt_cache"
    ))
  }
  return(res0)
}


#' @export
#'
#' @rdname KEGGhelp
database_KEGG_show_organism <- function() {
  x <- retry_function(
    read.delim,
    ntry = 3,
    delay = 3,
    "https://rest.kegg.jp/list/organism",
    quote = "\t",
    header = F
  )
  return(x)
}
