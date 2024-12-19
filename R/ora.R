parse_gs<-function(gset){
  colnames(gset) <- c("term", "gene")
  gset2<-split(gset,gset$term)
  gset2<-lapply(gset2,function(x)x$gene)
  gset2->gene_sets
  return(list(background_genes=0,gene_sets=gene_sets))
}

#' @importFrom qvalue qvalue
#' @importFrom stats p.adjust
#' @useDynLib EnrichGT
#' @importFrom Rcpp sourceCpp
doEnrich_Internal <- function(genes,database,p_adj_methods,p_val_cut_off,background_genes,min_geneset_size,max_geneset_size){

  tryCatch({
    if(ncol(database)!=2&ncol(database)!=3){
      cli::cli_abort("Not valid database")
    }
  },error=function(e){cli::cli_abort("Not valid database")})
  if(ncol(database)==3){
    colnames(database) <- c("ID","term","gene")
    db0 <- database[,c(1,2)]
    database <- database[,c(2,3)]
  }else{
    db0 <- data.frame(ID = database$term, term = database$term)
  }
  db0 <- db0 |> dplyr::mutate(CheckDup = paste0(ID,term)) |> dplyr::filter(!duplicated(CheckDup)) |> dplyr::select(-CheckDup) |> dplyr::rename(TERMs = term)
  colnames(database) <- c("term", "gene")

  termCount <- table(database$term) |> as.data.frame()
  colnames(termCount)[1] <- "TERMs"
  termCount <- termCount |> dplyr::filter(Freq >= min_geneset_size,Freq <= max_geneset_size)
  if(nrow(termCount)==0){
    cli::cli_abort("No gene set match the asked gene set size ")
  }
  termCount$TERMs <- as.character(termCount$TERMs)
  background_genes_<-database$gene |> table() |> names()
  database0 <-database |> dplyr::filter(term %in% termCount$TERMs)

  tgtGs<-parse_gs(database0)
  if(!is.null(background_genes)){
    tgtGs[[1]]<-background_genes
  }else{
    tgtGs[[1]]<-background_genes_
  }
  genes<-genes[genes%in%tgtGs[[1]]]

  HittedCheck<-data.frame(gene=genes,hitted="Hitted")
  database2 <- database0 |> left_join(HittedCheck,by="gene")
  database2 <- database2[!is.na(database2$hitted),]
  database2 <- database2[database2$hitted=="Hitted",]
  hitted_result <- database2 |>
    group_by(term) |>
    summarise(gene_list = paste(gene, collapse = "/")) |>
    ungroup()
  colnames(hitted_result)[1] <- "TERMs"
  hitted_result$TERMs <- as.character(hitted_result$TERMs)

  t1 <- Sys.time()
  if(length(genes)==0){
    cli::cli_abort("No gene overlap! Please recheck. ")
  }
  result <- ora_analysis(tgtGs[[2]], genes, tgtGs[[1]])
  result$TERMs <- tgtGs[[2]] |> names()
  df00 <- hitted_result |> left_join(result,by="TERMs") |> left_join(termCount,by="TERMs") |> left_join(db0,by="TERMs")
  df00$padj<-p.adjust(df00$PValue,method = p_adj_methods)
  qobj <- tryCatch(qvalue(p=df00$padj, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)
  if (inherits(qobj, "qvalue")) {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }
  df00$qvalues<-qvalues
  df00$LengthOfInput<-length(genes)
  df00$LengthOfBG<-length(tgtGs[[1]])
  res <- data.frame(ID=df00$ID,
                    Description=df00$TERMs,
                    GeneRatio=paste0(df00$Overlap,"/",df00$LengthOfInput),
                    BgRatio=paste0(df00$Freq,"/",df00$LengthOfBG),
                    pvalue=df00$PValue,
                    p.adjust=df00$padj,
                    qvalue=df00$qvalues,
                    geneID=df00$gene_list,
                    Count=df00$Overlap)
  t2 <- Sys.time()
  timeLast <- t2 - t1
  cli::cli_alert_success(paste0("Done ORA in ",timeLast," sec."))
  return(res)
}
