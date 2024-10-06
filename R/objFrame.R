setClass("EnrichGT_obj",slots=list(enriched_result="data.frame",
                                         gt_object="data.frame",
                                         gene_modules="list",
                                         pathway_clusters="list",
                                         hc_clustering_tree="list"))


setMethod("show", "EnrichGT_obj", function(object) {
  print(object@gt_object)
})

# gtobj<-new("EnrichGT_obj")
