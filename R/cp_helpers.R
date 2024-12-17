#' @title Return ranked gene list which is use for "GSEA"
#'
#' @param Gene A vector containing genes
#' @param Weight A vector contain weight of genes, typically like log2FC from DEG analysis
#'
#' @return A ranked named numeric vector. Names of the numbers is the ENTREZID.
#' @export
#' @author Zhiming Ye
#'
Ranked.GS<-function(Gene,Weight){
  log2FC <- Weight
  Genetable <- data.frame(Gene=Gene,log2FC=log2FC)
  Genetable <- Genetable |> dplyr::arrange(desc(log2FC))
  GSElist <- as.numeric(Genetable$log2FC)
  names(GSElist)<-Genetable$Gene
  GSElist = sort(GSElist, decreasing = TRUE)
  return(GSElist)
}


# gmt file reader is cited from Yu lab's gson package. Author: Guangchuang Yu
# https://github.com/YuLab-SMU/gson/blob/main/R/GMT.R

##' parse gmt file to a data.frame
##'
##' @title getGMTFile
##' @description
##' read `.gmt` files
##'
##' @param gmtfile gmt file
##' @importFrom utils stack
##' @export
##' @return data.frame
##' @author cited from https://github.com/YuLab-SMU/gson/blob/main/R/GMT.R
getGMTFile <- function (gmtfile) {
  x <- readLines(gmtfile)
  res <- strsplit(x, "\t")
  names(res) <- vapply(res, function(y) y[1], character(1))
  res <- lapply(res, "[", -c(1:2))
  ont2gene <- stack(res)
  ont2gene <- ont2gene[, c("ind", "values")]
  colnames(ont2gene) <- c("term", "gene")
  return(ont2gene)
}


#' A C++ accelerated universal enrichment analyzer
#'
#' @param genes a vector of gene id
#' @param database a database
#' @param p_adj_methods one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param p_val_cut_off adjusted pvalue cutoff on enrichment tests to report
#' @param background_genes background genes. If missing, the all genes listed in the database
#' @param min_geneset_size minimal size of genes annotated for testing
#' @param max_geneset_size maximal size of genes annotated for testing
#' @param multi_cores multi_cores
#'
#' @returns
#' @export
#' @examples
doEnrich <- function(genes,database,p_adj_methods="BH",p_val_cut_off=0.5,background_genes=NULL,min_geneset_size=10,max_geneset_size=500,multi_cores=0){
  if(is.character(genes)){
    result <- doEnrich_Internal(genes,database,p_adj_methods,p_val_cut_off,background_genes,min_geneset_size,max_geneset_size)
  }else if(is.list(genes)&multi_cores<=1){
    result <- lapply(genes,function(x){
      tryCatch({
        res <- doEnrich_Internal(genes=x,database,p_adj_methods,p_val_cut_off,background_genes,min_geneset_size,max_geneset_size)
        return(res)
      },error=function(e){
        return(data.frame(ERROR=e))
      })
    })
  }else if(is.list(genes)&multi_cores>=2){
    require(parallel)
    result <- mclapply(genes,function(x){
      tryCatch({
        res <- doEnrich_Internal(genes=x,database,p_adj_methods,p_val_cut_off,background_genes,min_geneset_size,max_geneset_size)
        return(res)
      },error=function(e){
        return(data.frame(ERROR=e))
      })
    },mc.cores = multi_cores)
  }
  return(result)
}

#' GSEA analysis by EnrichGT
#' @description
#' A warpper of fgsea::fgsea
#'
#' @param genes a vector of gene id
#' @param database a database
#' @param p_adj_methods one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param p_val_cut_off adjusted pvalue cutoff on enrichment tests to report
#' @param min_geneset_size minimal size of genes annotated for testing
#' @param max_geneset_size maximal size of genes annotated for testing
#' @param gseaParam pass to fgsea
#'
#' @returns
#' @export
#' @examples
doGSEA <- function(genes,database,min_geneset_size=10,max_geneset_size=500,gseaParam=1){
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
  database2 <- split(database$gene,database$term)
  fgseaRes <- fgsea::fgsea(pathways = database2,
                    stats    = genes,
                    minSize  = min_geneset_size,
                    maxSize  = max_geneset_size,
                    gseaParam=gseaParam)
  return(fgseaRes)
}
