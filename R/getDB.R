db_getter_env <- new.env()


# The Database cache system
#' @importFrom xfun md5
#' @importFrom AnnotationDbi select Ontology keys mapIds
#' @importFrom GO.db GOTERM
#' @importFrom reactome.db reactomeEXTID2PATHID reactomePATHID2EXTID reactomePATHID2NAME
UniversalInternalDBFetcher <- function(
  Type,
  OrgDB = org.Hs.eg.db,
  ONTOLOGY = NULL,
  ...
) {
  Spec <- OrgDB$packageName
  if (Type == "GO") {
    FUNInternal <- database_GO
    key <- "GO_"
    if (is.null(ONTOLOGY)) {
      cli::cli_abort("please provide ont.")
    }
  } else if (Type == "Reactome") {
    FUNInternal <- database_RA
    key <- "RA_"
  } else if (Type == "SelfBuild") {
    FUNInternal <- getGMT
    DBreaded <- readLines(ONTOLOGY)
    hashval <- xfun::md5(DBreaded)
    key <- "Self_"
    ONTOLOGY <- hashval
    OrgDB <- DBreaded
  }
  ObjName <- paste0(key, ONTOLOGY, "_", Spec)
  if (exists(ObjName, envir = db_getter_env)) {
    cli::cli_alert_success(paste0("Use cached database: ", ObjName))
    x <- get(ObjName, envir = db_getter_env)
  } else {
    assign("GETFUN", FUNInternal, envir = db_getter_env)
    x <- db_getter_env$GETFUN(OrgDB = OrgDB, ONTOLOGY = ONTOLOGY)
    assign(ObjName, x, envir = db_getter_env)
  }
  return(x)
}


#' @importFrom stringr str_to_title
dbParser <- function(DB, species) {
  type <- dplyr::case_when(
    (DB == "progeny" & species == "human") ~ "1",
    (DB == "progeny" & species == "mouse") ~ "2",
    (DB == "collectri" & species == "human") ~ "3",
    (DB == "collectri" & species == "mouse") ~ "4",
    T ~ "Error"
  )
  if (type == "Error") {
    cli::cli_abort("Please check your DB and species Input! ")
  }
  if (type == "1") {
    data("pws_human")
    tdb0 <- pws_human |>
      dplyr::filter(p_value < 0.1) |>
      dplyr::mutate(Direction = ifelse(weight > 0, "Up", "Down")) |>
      dplyr::mutate(TERM = paste0(source, "|", Direction)) |>
      dplyr::select(TERM, target)
  } else if (type == "2") {
    data("pws_mouse")
    tdb0 <- pws_mouse |>
      dplyr::filter(p_value < 0.1) |>
      dplyr::mutate(Direction = ifelse(weight > 0, "Up", "Down")) |>
      dplyr::mutate(TERM = paste0(source, "|", Direction)) |>
      dplyr::select(TERM, target)
    tdb0$target <- stringr::str_to_title(tdb0$target)
  } else if (type == "3") {
    data("TF_human")
    tdb0 <- TF_human |>
      dplyr::mutate(Direction = ifelse(mor > 0, "Up", "Down")) |>
      dplyr::mutate(TERM = paste0(source, "|", Direction)) |>
      dplyr::select(TERM, target)
  } else {
    data("TF_mouse")
    tdb0 <- TF_mouse |>
      dplyr::mutate(Direction = ifelse(mor > 0, "Up", "Down")) |>
      dplyr::mutate(TERM = paste0(source, "|", Direction)) |>
      dplyr::select(TERM, target)
    tdb0$target <- stringr::str_to_title(tdb0$target)
  }
  cli::cli_alert_success(paste0("success loaded self-contained database"))
  return(tdb0)
}


cvgs <- function(genes, from_what, to_what, orgDB) {
  genes <- genes |> as.character()
  loadNamespace("AnnotationDbi")
  x <- AnnotationDbi::select(
    orgDB,
    keys = genes,
    keytype = from_what,
    columns = c(from_what, to_what)
  )
  return(x)
}

#' Convert gene annotations from any keys to any keys
#'
#' @param genes gene vector
#' @param from_what input type (like "SYMBOL","ENTREZID","ENSEMBL","GENENAME",...), keys should be supported by AnnotationDbi. Search for the help page of AnnotationDbi for further help.
#' @param to_what output type (like "SYMBOL","ENTREZID","ENSEMBL","GENENAME",...), keys should be supported by AnnotationDbi. Search for the help page of AnnotationDbi for further help. Can be multiple items E.g. `c("ENTREZID","ENSEMBL","GENENAME")`
#' @param OrgDB human = org.Hs.eg.db, mouse = org.Mm.eg.db, search BioConductor website for further help
#'
#' @returns a data.frame
#' @export
#'
convert_annotations_genes <- function(genes, from_what, to_what, OrgDB) {
  assign("cvgs", cvgs, envir = db_getter_env)
  x <- db_getter_env$cvgs(genes, from_what, to_what, orgDB = OrgDB) # The case should be noticed
  return(x)
}

database_GO <- function(OrgDB, ONTOLOGY, ...) {
  t1 <- Sys.time()
  loadNamespace("dplyr")
  loadNamespace("tibble")
  loadNamespace("AnnotationDbi")
  loadNamespace("GO.db")

  goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
  if (ONTOLOGY != "ALL") {
    goterms <- goterms[goterms == ONTOLOGY]
  }
  go2gene <- suppressMessages(
    AnnotationDbi::mapIds(
      OrgDB,
      keys = names(goterms),
      column = "SYMBOL",
      keytype = c("GOALL"),
      multiVals = 'list'
    )
  )
  goAnno <- stack(go2gene)
  goAnno <- unique(goAnno[!is.na(goAnno[, 1]), ])
  go_terms <- AnnotationDbi::Term(GO.db::GOTERM)
  go_terms <- as.data.frame(go_terms) |> tibble::rownames_to_column(var = "ind")
  goAnno <- goAnno |>
    dplyr::left_join(go_terms, by = "ind") |>
    dplyr::select(ind, go_terms, values) |>
    na.omit()
  t2 <- Sys.time()
  cli::cli_alert_success(paste0(
    "success loaded database, time used : ",
    (t2 - t1),
    " sec."
  ))
  return(goAnno)
}

database_RA <- function(OrgDB, ...) {
  start_time <- Sys.time()

  for (pkg in c("AnnotationDbi", "reactome.db")) {
    loadNamespace(pkg)
  }

  entrez_keys <- AnnotationDbi::keys(OrgDB, keytype = "ENTREZID")

  gene_mapping <- AnnotationDbi::select(
    OrgDB,
    keys = entrez_keys,
    keytype = "ENTREZID",
    columns = c("ENTREZID", "SYMBOL")
  )

  extid_to_path <- as.list(reactome.db::reactomeEXTID2PATHID)
  path_to_extid <- as.list(reactome.db::reactomePATHID2EXTID)
  path_to_name <- as.list(reactome.db::reactomePATHID2NAME)

  relevant_extids <- intersect(names(extid_to_path), entrez_keys)
  extid_to_path <- extid_to_path[relevant_extids]

  path_names_clean <- vapply(path_to_name, function(x) {
    name <- x[1L]
    sub("^[A-Za-z]+\\s+[A-Za-z]+:\\s*", "", name)
  }, character(1L))
  pathway_ids <- unique(unlist(extid_to_path, use.names = FALSE))

  valid_paths <- intersect(pathway_ids, names(path_to_extid))
  valid_paths <- intersect(valid_paths, names(path_names_clean))

  path_to_extid_filtered <- path_to_extid[valid_paths]
  path_names_filtered <- path_names_clean[valid_paths]

  pathway_lengths <- lengths(path_to_extid_filtered)
  total_rows <- sum(pathway_lengths)

  result_ids <- rep(names(path_to_extid_filtered), pathway_lengths)
  result_entrez <- unlist(path_to_extid_filtered, use.names = FALSE)

  result_df <- data.frame(
    ID = result_ids,
    ENTREZID = result_entrez,
    stringsAsFactors = FALSE
  )
  result_df <- result_df[result_df$ENTREZID %in% gene_mapping$ENTREZID, ]
  pathway_df <- data.frame(
    ID = names(path_names_filtered),
    PN = unname(path_names_filtered),
    stringsAsFactors = FALSE
  )

  result_df <- merge(result_df, pathway_df, by = "ID", all.x = TRUE)
  result_df <- merge(result_df, gene_mapping, by = "ENTREZID", all.x = TRUE)
  final_result <- result_df[, c("ID", "PN", "SYMBOL")]

  end_time <- Sys.time()
  time_diff <- round(as.numeric(end_time - start_time), 3)

  if (requireNamespace("cli", quietly = TRUE)) {
    cli::cli_alert_success(paste0(
      "success loaded database, time used : ",
      time_diff, " sec."
    ))
  } else {
    message(paste0("Success: database loaded in ", time_diff, " seconds"))
  }

  return(final_result)
}

getGMT <- function(OrgDB, ONTOLOGY) {
  parsed_data <- strsplit(OrgDB, "\t", fixed = TRUE)
  term_names <- character(length(parsed_data))
  for (i in seq_along(parsed_data)) {
    term_names[i] <- parsed_data[[i]][1L]
  }
  names(parsed_data) <- term_names

  gene_lists <- vector("list", length(parsed_data))
  names(gene_lists) <- term_names
  for (i in seq_along(parsed_data)) {
    current_row <- parsed_data[[i]]
    if (length(current_row) > 2L) {
      gene_lists[[i]] <- current_row[3L:length(current_row)]
    } else {
      gene_lists[[i]] <- character(0)
    }
  }
  term_lengths <- lengths(gene_lists)
  total_rows <- sum(term_lengths)
  # Pre-allocate result vectors
  data.frame(
    term = rep(names(gene_lists), term_lengths),
    gene = unlist(gene_lists, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}


#' @rdname get_database
#' @name database_...
#' @title Get database for enrichment analysis
#' @description
#' Get Gene Ontology (GO), Reactome, and other term-to-gene database, for enrichment analysis
#'
#' @param OrgDB The AnnotationDbi database to fetch pathway data and convert gene IDs to gene symbols. For human it would be `org.Hs.eg.db`, for mouse it would be `org.Mm.eg.db`. In AnnotationDbi there are many species, please search `AnnotationDbi` for other species annotation database. GO and Reactome should add this, progeny and collectri do not.
#' @returns a data.frame with ID, terms and genes
#' @author Zhiming Ye. Part of functions were inspired by `clusterProfiler` but with brand new implement.
#' @export
database_GO_BP <- function(OrgDB = org.Hs.eg.db) {
  x <- UniversalInternalDBFetcher("GO", OrgDB, ONTOLOGY = "BP")
  return(x)
}
#' @rdname get_database
#' @export
database_GO_CC <- function(OrgDB = org.Hs.eg.db) {
  x <- UniversalInternalDBFetcher("GO", OrgDB, ONTOLOGY = "CC")
  return(x)
}
#' @rdname get_database
#' @export
database_GO_MF <- function(OrgDB = org.Hs.eg.db) {
  x <- UniversalInternalDBFetcher("GO", OrgDB, ONTOLOGY = "MF")
  return(x)
}
#' @rdname get_database
#' @export
database_GO_ALL <- function(OrgDB = org.Hs.eg.db) {
  x <- UniversalInternalDBFetcher("GO", OrgDB, ONTOLOGY = "ALL")
  return(x)
}
#' @rdname get_database
#' @export
database_Reactome <- function(OrgDB = org.Hs.eg.db) {
  x <- UniversalInternalDBFetcher("Reactome", OrgDB)
  return(x)
}
#' @rdname get_database
#' @export
database_progeny_human <- function() {
  return(dbParser("progeny", "human"))
}
#' @rdname get_database
#' @export
database_progeny_mouse <- function() {
  return(dbParser("progeny", "mouse"))
}
#' @rdname get_database
#' @export
database_CollecTRI_human <- function() {
  return(dbParser("collectri", "human"))
}
#' @rdname get_database
#' @export
database_CollecTRI_mouse <- function() {
  return(dbParser("collectri", "mouse"))
}
