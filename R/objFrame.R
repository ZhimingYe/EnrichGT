setClass("EnrichGT_obj",slots=list(enriched_result="data.frame",
                                         gt_object="list",
                                         gene_modules="list",
                                         pathway_clusters="list"))


setMethod("show", "EnrichGT_obj", function(object) {
  return(object@gt_object)
})

new.egt <- function(x1,x2,x3,x4){
  objegt<-new("EnrichGT_obj")
  objegt@enriched_result <- x1
  objegt@gt_object <- x2
  objegt@gene_modules <- x3
  objegt@pathway_clusters <- x4
  return(objegt)
}
