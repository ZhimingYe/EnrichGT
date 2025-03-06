setClass("EnrichGT_obj",slots=list(enriched_result="data.frame",
                                         gt_object="gt_tbl",
                                         gene_modules="list",
                                         pathway_clusters="list",
                                         document_term_matrix="dgCMatrix",
                                         clustering_tree="hclust",
                                         raw_enriched_result="data.frame",
                                         fused="logical"))


setMethod("show", "EnrichGT_obj", function(object) {
  print(object@gt_object)
})

new.egt <- function(x1,x2,x3,x4,x5,x6,x7){
  require("gt")
  flag0<-F
  tryCatch({
    objegt<-new("EnrichGT_obj")
    objegt@enriched_result <- x1
    objegt@gt_object <- x2
    objegt@gene_modules <- x3
    objegt@pathway_clusters <- x4
    objegt@document_term_matrix <- x5
    objegt@clustering_tree <- x6
    objegt@raw_enriched_result <- x7
    objegt@fused <- F
    flag0<-T
  },error=function(e){
    message_egt("Failed to create EnrichGT object! Please re-check your input.")
    flag0<-F
  })
  if(flag0==F){
    objegt<-NULL
  }
  return(objegt)
}

`%del%` <- function(x,y){
  if(class(x)=="EnrichGT_obj"){
    if(nrow(x@enriched_result)>2){
      x@enriched_result <- x@enriched_result |> dplyr::filter(!grepl(y,x@enriched_result$Description))
      if(sum(colnames(x@enriched_result)=="NES")>0){
        xx <- x@enriched_result |> dplyr::mutate(Reg=ifelse(NES>0,"red","forestgreen"))
        x@gt_object <-xx |> gt_gsea(ClusterNum="previous",objname=glue::glue("remove:{x}"),...)
        obj3 <- x@enriched_result |> genMetaGM(type="GSEA")
        x@gene_modules <- obj3[[1]]
        x@pathway_clusters <- obj3[[2]]
      }
      return(x)
    }
  }
  else if(is.data.frame(x)){
      x <- x |> dplyr::filter(!grepl(y,x$Description))
      return(x)
  }
  else{
    cli::cli_abort("Please provide enriched results. ")
  }
  
}