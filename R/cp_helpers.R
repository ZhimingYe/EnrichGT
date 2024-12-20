#' @title Return ranked gene list which is use for "GSEA" or other places
#'
#' @param genes A vector containing genes
#' @param weights A vector contain weight of genes, typically like log2FC from DEG analysis
#'
#' @return A ranked named numeric vector. Names of the numbers is the ENTREZID.
#' @export
#' @author Zhiming Ye
#'
genes_with_weights<-function(genes,weights){
  genes -> Gene
  weights -> Weight
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

# parse gmt form MsigDB file to a data.frame
#
#' @title parse GMT files file to a data.frame
#' @description
#' Read `.gmt` files. You can download them from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
#'
#' WikiPathway database also provides pre-built GMT files (https://data.wikipathways.org/current/gmt/). In default they are recorded as ENTREZ IDs, so you need to provide proper species database (e.g. org.Hs.eg.db for human), to database_from_gmt function and EnrichGT will automatically convert ENTREZ IDs to gene symbols for enrichment analysis.
#' @param gmtfile gmt file path
#' @param OrgDB Only need when converting genes, human = org.Hs.eg.db, mouse = org.Mm.eg.db, search BioConductor website for further help. Default is NULL.
#' @param convert_2_symbols Force to convert numeric gene ids (as ENTREZIDs) to gene symbols
#' @importFrom utils stack
#' @export
#' @return data.frame
#' @author cited from https://github.com/YuLab-SMU/gson/blob/main/R/GMT.R . The further Cache system is written by Zhiming Ye.

database_from_gmt <- function (gmtfile,OrgDB = NULL,convert_2_symbols=T) {
  x <- UniversalInternalDBFetcher("SelfBuild",NULL,gmtfile)
  if(convert_2_symbols){
    if(is.null(OrgDB)){
      cli::cli_abort("No valid OrgDB, please load it and input")
    }
    x <- x |> as_tibble()
    if(ncol(x)>3){
      cli::cli_abort("Not valid gmt file! ")
    }
    if(nrow(x)<20){
      Target <- nrow(x)
    }else{
      Target <-20
    }
    if(((is_numeric_string(x[1:Target,ncol(x),drop=T]) |> sum())/Target)>0.75){
      cli::cli_alert_info("Detectd numeric gene ids, try to convert to symols.\n You can set convert_2_symbols=F to disable this.")
      colnames(x)[ncol(x)]<-"ENTREZID"
      tryCatch({
        qq <- convert_annotations_genes(x$ENTREZID,from_what = "ENTREZID",to_what = "SYMBOL",OrgDB = OrgDB)
        if(nrow(qq |> na.omit())<5){
          cli::cli_abort("Error in convert. Please set convert_2_symbols=F to disable this")
        }

        x <- x |> left_join(qq,by="ENTREZID")
        x <- x |> dplyr::select(-ENTREZID)
      },error=function(e){
        cli::cli_abort("Error in convert. Please set convert_2_symbols=F to disable this")
      })

    }
  }

  return(x)
}


#' A C++ accelerated universal enrichment analyzer (Over-Representation Analysis (ORA))
#'
#' @description
#' ORA is a statistical method used to identify biological pathways or gene sets that are significantly enriched in a given list of genes (e.g., differentially expressed genes). The method compares the proportion of genes in the target list that belong to a specific category (e.g., pathways, GO terms) to the expected proportion in the background gene set.
#'
#' To accelerate the computation in ORA analysis, `EnrichGT` have implemented a function that leverages C++ for high-performance computation. The core algorithm utilizes hash tables for efficient lookup and counting of genes across categories. Also It provides multi-Core parallel calculation by package `parallel`.
#'
#' @usage res <- egt_enrichment_analysis(genes = DEGtable$Genes,
#' database = database_GO_BP())
#'
#' @usage res <- egt_enrichment_analysis(genes = c("TP53","CD169","CD68","CD163",...),
#' database = database_GO_ALL())
#'
#' @usage res <- egt_enrichment_analysis(genes = c("TP53","CD169","CD68","CD163",...),
#' database = database_from_gmt("MsigDB_Hallmark.gmt"))
#'
#' @usage res <- egt_enrichment_analysis(list(Macrophages=c("CD169","CD68","CD163"),
#' Fibroblast=c("COL1A2","COL1A3"),...),
#'  database = database_from_gmt("panglaoDB.gmt"))
#'
#' @param genes a vector of gene ids like `c("TP53","CD169","CD68","CD163"...)`.
#'
#' If you have genes from multiple source or experiment group, you can also pass a list with gene ids in it. For Example , `list(Macrophages=c("CD169","CD68","CD163"),Fibroblast=c("COL1A2","COL1A3))`.
#'
#' The genes should be match in the second param `database`'s `gene` column. For example, if database provides Ensembl IDs, you should input Ensembl IDs. But in default databases provided by `EnrichGT` is gene symbols.
#'
#' @param database a database data frame, can contain 3 columns (ID, Pathway_Name, Genes) or just 2 columns (Pathway_Name, Genes). You can read a data frame and pass it through this or run `database_GO_CC()` to get them, see example.
#'
#' You can run `database_GO_CC()` to see an example.
#'
#' The ID column is not necessary. EnrichGT contains several databases, functions about databases are named starts with `database_...`, like `database_GO_BP()` or `database_Reactome()`.
#'
#' The default gene in each database EnrichGT provided to input is `GENE SYMBOL` (like TP53, not 1256 or ENSG...), not `ENTREZ ID` or `Ensembl ID`.
#'
#' It will be more convince for new users. Avaliable databases includes `database_GO_BP()`, `database_GO_CC()`, `database_GO_MF()` and `database_Reactome()`.
#'
#' You can add more database by downloading MsigDB (https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)'s GMT files. It can be load by using `database_from_gmt(FILE_PATH)`.
#'
#' If you only have simple a table, you can also pass a data frame through this arguement. Of note, it should contains at least 2 coloumn (colnames(df) = c("Terms","Genes)), the first is term names and the second are the corresponding genes. If you have term ids, you can add a `ID` column at the first column, and `Terms` becomes the second column and `Genes` the third.
#'
#' @param p_adj_methods one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param p_val_cut_off adjusted pvalue cutoff on enrichment tests to report
#' @param background_genes background genes. If missing, the all genes listed in the database
#' @param min_geneset_size minimal size of genes annotated for testing
#' @param max_geneset_size maximal size of genes annotated for testing
#' @param multi_cores multi_cores (Experimental), only used when analysis a list of genes (multi-source or groups). Set to 0 or 1 to disable. May use lots of RAM.
#'
#' @returns a data frame with ORA results.
#' @export
#'
#' @author Zhiming Ye
egt_enrichment_analysis <- function(genes,database,p_adj_methods="BH",p_val_cut_off=0.5,background_genes=NULL,min_geneset_size=10,max_geneset_size=500,multi_cores=0){
  if(is.character(genes)){
    result <- doEnrich_Internal(genes,database,p_adj_methods,p_val_cut_off,background_genes,min_geneset_size,max_geneset_size)
  }else if(is.list(genes)&multi_cores<=1){
    result <- lapply(genes,function(x){
      tryCatch({
        res <- doEnrich_Internal(genes=x,database,p_adj_methods,p_val_cut_off,background_genes,min_geneset_size,max_geneset_size)
        return(res)
      },error=function(e){
        return(data.frame(ERROR="error..."))
      })
    })
  }else if(is.list(genes)&multi_cores>=2){
    require(parallel)
    result <- mclapply(genes,function(x){
      tryCatch({
        res <- doEnrich_Internal(genes=x,database,p_adj_methods,p_val_cut_off,background_genes,min_geneset_size,max_geneset_size)
        return(res)
      },error=function(e){
        return(data.frame(ERROR="error..."))
      })
    },mc.cores = multi_cores)
  }
  tryCatch({
    if(is.data.frame(result)){
      result <- result |> dplyr::filter(pvalue<p_val_cut_off) |> dplyr::arrange(pvalue)
    }else if(is.list(result)){
      result <- lapply(result,function(x){
        if(is.data.frame(x)){
          if(ncol(x)>3 & nrow(x)>3){
            z <- x |> dplyr::filter(pvalue<p_val_cut_off) |> dplyr::arrange(pvalue)
            return(z)
          }
        }

        })
    }

  },error=function(e){
    message_wrong_db()
    cli::cli_abort("No useable result!")
  })
  return(result)
}

#' Gene Set Enrichment Analysis (GSEA) by EnrichGT
#' @description
#' A warpper of `fgsea::fgsea()`.
#'
#' GSEA is a computational method used to determine whether predefined gene sets (e.g., pathways, GO terms) are statistically enriched in a ranked list of genes. Unlike ORA, GSEA considers the entire gene list and focuses on the cumulative distribution of gene ranks to identify coordinated changes.
#'
#' The fgsea (https://github.com/ctlab/fgsea) package is an R tool that implements an accelerated version of GSEA. It uses precomputed statistical methods and efficient algorithms to dramatically speed up enrichment analysis, especially for large datasets.
#'
#' @usage res <- egt_gsea_analysis(genes = genes_with_weights(genes = DEG$genes, weights = DEG$log2FoldChange),
#' database = database_GO_BP())
#'
#' @usage res <- egt_gsea_analysis(genes = genes_with_weights(genes = PCA_res$genes,weights =PCA_res$PC1_loading),
#' database = database_from_gmt("MsigDB_Hallmark.gmt"))
#' @param genes a named numeric vector, for example c(`TP53`=1.2,`KRT15`=1.1,`IL1B`=1.0,`PMP22` = 0.5,`FABP1` = -0.9, `GLUT1` = -2).
#'
#' The number is the weight of each gene, can use the logFC form DEG analysis results instead. Also NMF or PCA's loading can also be used.
#'
#' `EnrichGT` provides a `genes_with_weights(genes,weights)` function to build this numeric vector. Importantly, this vector should be !SORTED! for larger to smaller.
#'
#' @param database a database data frame, can contain 3 columns (ID, Pathway_Name, Genes) or just 2 columns (Pathway_Name, Genes). You can read a data frame and pass it through this or run `database_GO_CC()` to get them, see example.
#'
#' The ID column is not necessary. EnrichGT contains several databases, functions about databases are named starts with `database_...`, like `database_GO_BP()` or `database_Reactome()`.
#'
#' The default gene in each database EnrichGT provided to input is `GENE SYMBOL` (like TP53, not 1256 or ENSG...), not `ENTREZ ID` or `Ensembl ID`. It will be more convince for new users.
#'
#' Avaliable databases includes `database_GO_BP()`, `database_GO_CC()`, `database_GO_MF()` and `database_Reactome()`.
#'
#' You can add more database by downloading MsigDB (https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)'s GMT files. It can be load by using `database_from_gmt(FILE_PATH)`.
#'
#' If you only have simple a table, you can also pass a data frame through this arguement. Of note, it should contains at least 2 coloumn (colnames(df) = c("Terms","Genes)), the first is term names and the second are the corresponding genes. If you have term ids, you can add a `ID` column at the first column, and `Terms` becomes the second column and `Genes` the third.
#'
#' @param p_adj_methods one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param p_val_cut_off adjusted pvalue cutoff on enrichment tests to report
#' @param min_geneset_size minimal size of genes annotated for testing
#' @param max_geneset_size maximal size of genes annotated for testing
#' @param gseaParam other param passing to fgsea
#'
#' @returns a data frame
#' @export
#' @author warpped from fgsea package.
#'
egt_gsea_analysis <- function(genes,database,p_val_cut_off=0.5,min_geneset_size=10,max_geneset_size=500,gseaParam=1){
  if(!is.list(genes)){
    res <- egt_gsea_analysis_internal(genes=genes,database=database,p_val_cut_off=p_val_cut_off,min_geneset_size=min_geneset_size,max_geneset_size=max_geneset_size,gseaParam=gseaParam)
  }else{
    res <- lapply(genes,function(x){
      tryCatch({
        output <- egt_gsea_analysis_internal(genes=x,database=database,p_val_cut_off=p_val_cut_off,min_geneset_size=min_geneset_size,max_geneset_size=max_geneset_size,gseaParam=gseaParam)
      },error=function(e){
        output <- NULL
      })
      return(output)
    })
  }
  return(res)
}

egt_gsea_analysis_internal <- function(genes,database,p_val_cut_off=0.5,min_geneset_size=10,max_geneset_size=500,gseaParam=1){
  tryCatch({
    if(ncol(database)!=2&ncol(database)!=3){
      message_wrong_db()
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
  db0 <- db0 |> dplyr::mutate(CheckDup = paste0(ID,term)) |> dplyr::filter(!duplicated(CheckDup)) |> dplyr::select(-CheckDup) |> dplyr::rename(pathway = term) # Because of output is pathway
  colnames(database) <- c("term", "gene")
  database2 <- split(database$gene,database$term)
  t1 <- Sys.time()
  fgseaRes <- fgsea::fgsea(pathways = database2,
                    stats    = genes,
                    minSize  = min_geneset_size,
                    maxSize  = max_geneset_size,
                    gseaParam=gseaParam)
  fgseaRes <- fgseaRes |> left_join(db0,by="pathway")
  t2 <- Sys.time()
  cli::cli_alert_success(paste0("Sucessful GSEA, time last ", (t2-t1)," secs."))
  fgseaRes2 <- data.frame(ID = fgseaRes$ID,
                          Description = fgseaRes$pathway,
                          ES = fgseaRes$ES,
                          NES = fgseaRes$NES,
                          pvalue = fgseaRes$pval,
                          p.adjust = fgseaRes$padj,
                          core_enrichment = sapply(fgseaRes$leadingEdge,function(x)paste(x,collapse ="/")))
  tryCatch({
    fgseaRes2 <- fgseaRes2 |> dplyr::filter(pvalue<p_val_cut_off) |> dplyr::arrange(desc(NES))
  },error=function(e){
    message_wrong_db()
    cli::cli_abort("No useable result!")
  })
  return(fgseaRes2)
}

message_wrong_db<-function(){
  cli::cli_alert_danger("Please re-check the database OrgDB argument or input GMT file matchs the inputed genes. Some time you may forget to run: ")
  cli::cli_code("library(org.Hs.eg.db) # for humans, if mouse use org.Mm.eg.db. Please re-check")
  cli::cli_code("See https://zhimingye.github.io/EnrichGT/#built-in-database-form-annotationdbi for further details")
}
