#' Infering Pathway or Transcript Factors activity from EnrichGT meta-gene modules
#' @description
#' Only supports gene symbols. so you must use enrichedObj |> setReadable(OrgDb = xxx,keyType = "ENTREZID") |> EnrichGT() . Do Not supports ENTREZIDs!
#'
#' [PROGENy](https://saezlab.github.io/progeny/) is a comprehensive resource containing a curated collection of pathways and their target genes, with weights for each interaction.
#'
#'[CollecTRI](https://github.com/saezlab/CollecTRI) is a comprehensive resource containing a curated collection of TFs and their transcriptional targets compiled from 12 different resources. This collection provides an increased coverage of transcription factors and a superior performance in identifying perturbed TFs compared to our previous.
#' 
#' If when doing re-enrichment, you select a high number of clusters, that may cause low gene number in each meta-gene module, and then can't be infered sucessfully. So if result is empty, please increase the number of re-clustering when doing it. 
#'
#' @param x an EnrichGT_obj object.
#' @param DB can be "progeny" (the Pathway activity database), or "collectri" (TF activity database)
#' @param species can be "human" or "mouse"
#'
#' @return a `compareCluster` result from `clusterProfiler`
#' @export
#' @importFrom stringr str_to_title
#' @author Zhiming Ye, Saez-Rodriguez Lab (The decoupleR package, https://saezlab.github.io/decoupleR/)
#'
egt_infer_act <- function(x,DB="collectri",species="human"){
  genelist <- tryCatch(x@gene_modules,error=function(e){
    cli::cli_abort("NOT an enrichGT object. ")
  })
  cli::cli_alert_warning("If when doing re-enrichment, you select a high number of clusters, that may cause low gene number in each meta-gene module, and then can't be infered sucessfully. So if result is empty, please increase the number of re-clustering when doing it. ")
  InferedCpres <- egt_enrichment_analysis(genelist,p_val_cut_off = 1 ,min_geneset_size=2,max_geneset_size=1000,database=(dbParser(DB,species)))
  return(InferedCpres)
}
