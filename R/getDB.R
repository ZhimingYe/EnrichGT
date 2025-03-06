db_getter_env <- new.env()


# The Database cache system
#' @importFrom xfun md5
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
#' @examples
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


# Function fetching Reactome is cited from https://github.com/YuLab-SMU/ReactomePA/blob/devel/R/gseAnalyzer.R
# with several modifies
# (c) Guangchuang Yu @ SMU ReactomePA
database_RA <- function(OrgDB, ...) {
  t1 <- Sys.time()
  loadNamespace("dplyr")
  loadNamespace("tibble")
  loadNamespace("AnnotationDbi")
  loadNamespace("reactome.db")
  eg <- AnnotationDbi::keys(OrgDB, keytype = c("ENTREZID"))
  eg2 <- AnnotationDbi::select(
    OrgDB,
    keys = eg,
    keytype = "ENTREZID",
    columns = c("ENTREZID", "SYMBOL")
  )
  EXTID2PATHID <- as.list(reactome.db::reactomeEXTID2PATHID)
  EXTID2PATHID <- EXTID2PATHID[names(EXTID2PATHID) %in% eg]
  PATHID2EXTID <- as.list(reactome.db::reactomePATHID2EXTID) ## also contains reactions
  PATHID2NAME <- as.list(reactome.db::reactomePATHID2NAME)
  PI <- names(PATHID2NAME)
  ## > PATHID2NAME[['68877']]
  ## [1] "Homo sapiens: Mitotic Prometaphase" "Homo sapiens: Mitotic Prometaphase"
  PATHID2NAME <- lapply(PATHID2NAME, function(x) x[1])
  names(PATHID2NAME) <- PI
  PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% names(PATHID2NAME)]
  PATHID2EXTID <- PATHID2EXTID[
    names(PATHID2EXTID) %in% unique(unlist(EXTID2PATHID))
  ]
  PATHID2NAME <- PATHID2NAME[names(PATHID2NAME) %in% names(PATHID2EXTID)]
  PATHID2NAME <- unlist(PATHID2NAME)
  PATHID2NAME <- gsub("^\\w+\\s\\w+:\\s+", "", PATHID2NAME) # remove leading spaces
  df <- data.frame(
    List_Name = rep(names(PATHID2EXTID), sapply(PATHID2EXTID, length)),
    Element = unlist(PATHID2EXTID)
  )
  df <- df[df$Element %in% eg2$ENTREZID, ]
  PATHID2NAME <- as.data.frame(PATHID2NAME) |>
    tibble::rownames_to_column(var = "ID")
  colnames(df) <- c("ID", "ENTREZID")
  df <- df |> dplyr::left_join(PATHID2NAME, by = "ID")
  df <- df |> dplyr::left_join(eg2, by = "ENTREZID")
  df <- df[, c(1, 3, 4)]
  t2 <- Sys.time()
  cli::cli_alert_success(paste0(
    "success loaded database, time used : ",
    (t2 - t1),
    " sec."
  ))
  return(df)
}


getGMT <- function(OrgDB, ONTOLOGY) {
  x <- OrgDB
  res <- strsplit(x, "\t")
  names(res) <- vapply(res, function(y) y[1], character(1))
  res <- lapply(res, "[", -c(1:2))
  ont2gene <- stack(res)
  ont2gene <- ont2gene[, c("ind", "values")]
  colnames(ont2gene) <- c("term", "gene")
  return(ont2gene)
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
