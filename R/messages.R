
message_egt<-function(x,Type=0){
  cli_h1("EnrichGT message:")
  cli_alert_info(x)
}

nsimp<-function(){
  message_egt("\n=====[SUGGESTION]=====\nYou are passing an object from GO Enrichment.\nPlease ensure that `obj |> clusterProfiler::simplify()` is executed, to pre-simplify result,\nFor better enriched result.\n")
}
