

db_getter_env<-new.env()

database_GO <- function(ONTOLOGY, OrgDB) {
  loadNamespace("dplyr")
  loadNamespace("tibble")
  loadNamespace("AnnotationDbi")
  loadNamespace("GO.db")


  go_annotations <- AnnotationDbi::select(
    OrgDB,
    keys = keys(OrgDB, keytype = "SYMBOL"),
    columns = c("GO", "ONTOLOGY"),
    keytype = "SYMBOL"
  )
  if (ONTOLOGY != "ALL"){
    bp_annotations <- go_annotations |> dplyr::filter(ONTOLOGY == ONTOLOGY)
  }

  bp_annotations <- na.omit(bp_annotations)

  go_terms <- AnnotationDbi::Term(GO.db::GOTERM)
  go_terms <- as.data.frame(go_terms) |>
    rownames_to_column(var = "GO") |>
    dplyr::left_join(bp_annotations, by = "GO") |>
    dplyr::select(GO, go_terms, SYMBOL) |>
    na.omit()

  return(go_terms)
}


# Function fetching Reactome is cited from https://github.com/YuLab-SMU/ReactomePA/blob/devel/R/gseAnalyzer.R
# with several modifies
# (c) Guangchuang Yu @ SMU ReactomePA
database_RA <- function(OrgDB) {
  loadNamespace("dplyr")
  loadNamespace("tibble")
  loadNamespace("AnnotationDbi")
  loadNamespace("reactome.db")
  eg <- AnnotationDbi::keys(org.Hs.eg.db, keytype=c("ENTREZID"))
  eg2 <- AnnotationDbi::select(org.Hs.eg.db, keys = eg,
                               keytype = "ENTREZID", columns = c("ENTREZID", "SYMBOL"))
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
  PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% unique(unlist(EXTID2PATHID))]
  PATHID2NAME <- PATHID2NAME[names(PATHID2NAME) %in% names(PATHID2EXTID)]
  PATHID2NAME <- unlist(PATHID2NAME)
  PATHID2NAME <- gsub("^\\w+\\s\\w+:\\s+", "", PATHID2NAME) # remove leading spaces
  df <- data.frame(
    List_Name = rep(names(PATHID2EXTID), sapply(PATHID2EXTID, length)),
    Element = unlist(PATHID2EXTID)
  )
  df <- df[df$Element%in%eg2$ENTREZID,]
  PATHID2NAME<-as.data.frame(PATHID2NAME) |> rownames_to_column(var="ID")
  colnames(df)<-c("ID","ENTREZID")
  df<-df |> left_join(PATHID2NAME,by="ID")
  df<-df |> left_join(eg2,by="ENTREZID")
  df<-df[,c(1,3,4)]
  return(df)
}



dbParser<-function(DB,species){
  type<-case_when((DB=="progeny"&species=="human")~"1",
                  (DB=="progeny"&species=="mouse")~"2",
                  (DB=="collectri"&species=="human")~"3",
                  (DB=="collectri"&species=="mouse")~"4",
                  T~"Error"
  )
  if(type=="Error"){
    cli::cli_abort("Please check your DB and species Input! ")
  }
  if(type=="1"){
    data("pws_human")
    tdb0<-pws_human |> dplyr::filter(p_value<0.1) |> dplyr::mutate(Direction=ifelse(weight>0,"Up","Down")) |> dplyr::mutate(TERM=paste0(source,"|",Direction)) |> dplyr::select(TERM,target)
  }
  else if(type=="2"){
    data("pws_mouse")
    tdb0<-pws_mouse |> dplyr::filter(p_value<0.1) |> dplyr::mutate(Direction=ifelse(weight>0,"Up","Down")) |> dplyr::mutate(TERM=paste0(source,"|",Direction)) |> dplyr::select(TERM,target)
    tdb0$target <- str_to_title(tdb0$target)
  }
  else if(type=="3"){
    data("TF_human")
    tdb0<-TF_human |> dplyr::mutate(Direction=ifelse(mor>0,"Up","Down")) |> dplyr::mutate(TERM=paste0(source,"|",Direction)) |> dplyr::select(TERM,target)
  }
  else{
    data("TF_mouse")
    tdb0<-TF_mouse |> dplyr::mutate(Direction=ifelse(mor>0,"Up","Down")) |> dplyr::mutate(TERM=paste0(source,"|",Direction)) |> dplyr::select(TERM,target)
    tdb0$target <- str_to_title(tdb0$target)
  }
  return(tdb0)
}





#' @rdname get_database
#' @title Get database form Gene ontology or Reactome Pathways
#' @param OrgDB org.DB form bioconductor, can be org.Hs.eg.db or org.Mm.eg.db,... GO and Reactome should add this, progeny and collectri do not.
#' @returns a data.frame with ID, terms and genes
#' @export
database_GO_BP <- function(OrgDB=org.Hs.eg.db){
  assign("database_GO",database_GO,envir = db_getter_env)
  eval({
    suppressWarnings(try({rm(x)}))
    x <- database_GO(ONTOLOGY="BP",OrgDB=OrgDB)
  },envir = db_getter_env)
  x <- db_getter_env$x
  return(x)
}
#' @rdname get_database
#' @export
database_GO_CC <- function(OrgDB=org.Hs.eg.db){
  assign("database_GO",database_GO,envir = db_getter_env)
  eval({
    suppressWarnings(try({rm(x)}))
    x <- database_GO(ONTOLOGY="CC",OrgDB=OrgDB)
  },envir = db_getter_env)
  x <- db_getter_env$x
  return(x)
}
#' @rdname get_database
#' @export
database_GO_MF <- function(OrgDB=org.Hs.eg.db){
  assign("database_GO",database_GO,envir = db_getter_env)
  eval({
    suppressWarnings(try({rm(x)}))
    x <- database_GO(ONTOLOGY="MF",OrgDB=OrgDB)
  },envir = db_getter_env)
  x <- db_getter_env$x
  return(x)
}
#' @rdname get_database
#' @export
database_GO_ALL <- function(OrgDB=org.Hs.eg.db){
  assign("database_GO",database_GO,envir = db_getter_env)
  eval({
    suppressWarnings(try({rm(x)}))
    x <- database_GO(ONTOLOGY="ALL",OrgDB=OrgDB)
  },envir = db_getter_env)
  x <- db_getter_env$x
  return(x)
}
#' @rdname get_database
#' @export
database_Reactome <- function(OrgDB=org.Hs.eg.db){
  assign("database_RA",database_RA,envir = db_getter_env)
  eval({
    suppressWarnings(try({rm(x)}))
    x <- database_RA(OrgDB=OrgDB)
  },envir = db_getter_env)
  x <- db_getter_env$x
  return(x)
}
#' @rdname get_database
#' @export
database_progeny_human <- function(){
  return(dbParser("progeny","human"))
}
#' @rdname get_database
#' @export
database_progeny_mouse <- function(){
  return(dbParser("progeny","mouse"))
}
#' @rdname get_database
#' @export
database_CollecTRI_human <- function(){
  return(dbParser("collectri","human"))
}
#' @rdname get_database
#' @export
database_CollecTRI_mouse <- function(){
  return(dbParser("collectri","mouse"))
}
