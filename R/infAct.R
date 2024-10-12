#' Infering Pathway or Transcript Factors activity from EnrichGT meta-gene modules
#' @description
#' Only supports gene symbols. so you must use enrichedObj |> setReadable(OrgDb = xxx,keyType = "ENTREZID") |> EnrichGT() . Do Not supports ENTREZIDs!
#'
#' [PROGENy](https://saezlab.github.io/progeny/) is a comprehensive resource containing a curated collection of pathways and their target genes, with weights for each interaction.
#'
#'[CollecTRI](https://github.com/saezlab/CollecTRI) is a comprehensive resource containing a curated collection of TFs and their transcriptional targets compiled from 12 different resources. This collection provides an increased coverage of transcription factors and a superior performance in identifying perturbed TFs compared to our previous.
#'
#' @param x an EnrichGT_obj object.
#' @param DB can be "progeny" (the Pathway activity database), or "collectri" (TF activity database)
#' @param species can be "human" or "mouse"
#'
#' @return a `compareCluster` result from `clusterProfiler`
#' @export
#' @author Zhiming Ye, Saez-Rodriguez Lab (The decoupleR package, https://saezlab.github.io/decoupleR/)
#'
infering_regulator_act <- function(x,DB="collectri",species="human"){
  genelist <- tryCatch(x@gene_modules,error=function(e){
    cli::cli_abort("NOT an enrichGT object. ")
  })
  if(!require("clusterProfiler")){cli::cli_abort("Please install clusterProfiler form BioCondutor. ")}
  InferedCpres <- clusterProfiler::compareCluster(geneCluster = genelist, fun = clusterProfiler::enricher,pvalueCutoff = 1 ,qvalueCutoff = 1 ,minGSSize=2,TERM2GENE=(dbParser(DB,species)))
  return(InferedCpres)
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
    tdb0$target<-stringr::str_to_title(tdb0$target)
  }
  else if(type=="3"){
    data("TF_human")
    tdb0<-TF_human |> dplyr::mutate(Direction=ifelse(mor>0,"Up","Down")) |> dplyr::mutate(TERM=paste0(source,"|",Direction)) |> dplyr::select(TERM,target)
  }
  else{
    data("TF_mouse")
    tdb0<-TF_mouse |> dplyr::mutate(Direction=ifelse(mor>0,"Up","Down")) |> dplyr::mutate(TERM=paste0(source,"|",Direction)) |> dplyr::select(TERM,target)
    tdb0$target<-stringr::str_to_title(tdb0$target)
  }
  return(tdb0)
}




