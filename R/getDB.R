

db_getter_env <- new.env()
db_getter_env$database_GO <- function(ONTOLOGY, OrgDB) {
  require(dplyr)
  require(tibble)
  require(AnnotationDbi)
  require(GO.db)


  go_annotations <- AnnotationDbi::select(
    OrgDB,
    keys = keys(OrgDB, keytype = "SYMBOL"),
    columns = c("GO", "ONTOLOGY"),
    keytype = "SYMBOL"
  )

  bp_annotations <- go_annotations |> dplyr::filter(ONTOLOGY == ONTOLOGY)
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
db_getter_env$database_RA <- function(OrgDB) {
  require(dplyr)
  require(tibble)
  require(AnnotationDbi)
  require(reactome.db)
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
  df<-df |> left_join(eg2)
  df<-df[,c(1,3,4)]
  return(df)
}


database_GO_BP <- function(OrgDB=org.Hs.eg.db){
  result <- db_getter_env$database_GO(ONTOLOGY="BP",OrgDB=OrgDB)
  return(result)
}

database_GO_CC <- function(OrgDB=org.Hs.eg.db){
  result <- db_getter_env$database_GO(ONTOLOGY="CC",OrgDB=OrgDB)
  return(result)
}

database_GO_MF <- function(OrgDB=org.Hs.eg.db){
  result <- db_getter_env$database_GO(ONTOLOGY="MF",OrgDB=OrgDB)
  return(result)
}

database_Reactome <- function(OrgDB=org.Hs.eg.db){
  result <- db_getter_env$database_RA(OrgDB=OrgDB)
  return(result)
}


