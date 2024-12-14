egtDotPlot <- function(x,...){
  return(figure0)
}

ORA2dp<-function(x,clusterType="Empty",showStart=1,showEnd=20,orderBy,col_low="#1e6091",col_high="#99d98c",...){
  require(forcats)
  if(showEnd>nrow(x@Results)){
    if(nrow(x@Results)>30){
      showEnd<-30
    }
    else{
      showEnd<-nrow(x@Results)
    }
  }
  clusterType <- match.arg(clusterType, c("MatrixResult", "NetworkResult","Empty"))
  orderBy<-match.arg(orderBy,c("pValue","GeneRatio","Counts"))
  if(sum(c("GSCluster","NetworkCluster")%in%colnames(x@Results))==0|clusterType=="Empty"){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,Counts=x@Results$Count,GeneRatio=x@Results$GeneRatio,Cluster=NA,Regulation=x@Results$RegType)
  }
  if("GSCluster"%in%colnames(x@Results)&clusterType=="MatrixResult"){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,Counts=x@Results$Count,GeneRatio=x@Results$GeneRatio,Cluster=x@Results$GSCluster,Regulation=x@Results$RegType)
    FigDF$Cluster[is.na(FigDF$Cluster)]<-"UnKnown"
  }
  if("NetworkCluster"%in%colnames(x@Results)&clusterType=="NetworkResult"){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,Counts=x@Results$Count,GeneRatio=x@Results$GeneRatio,Cluster=x@Results$NetworkCluster,Regulation=x@Results$RegType)
    FigDF$Cluster[is.na(FigDF$Cluster)]<-"UnKnown"
  }
  if(nrow(FigDF)<2){
    stop("Enrichment result is too few!\n")
  }
  FigDF<-FigDF%>%dplyr::arrange(desc(!!sym(orderBy)))
  if(showEnd<showStart){
    stop("End is smaller than start!\n")
  }
  keep<-c(1:nrow(FigDF))[c(1:nrow(FigDF))%in%c(showStart:showEnd)]
  if(length(keep)==0){
    stop("NO overlap!\n")
  }
  p<-ggplot(FigDF[keep,],aes(x=parse_ratio_(GeneRatio), y=fct_reorder(Terms,!!sym(orderBy)), size=Counts, color=pValue,shape=Regulation))+geom_point()+scale_color_continuous(low=col_low, high=col_high, name = "adjustedP",guide=guide_colorbar(reverse=F))+scale_size(range=c(2, 8))+xlab("Gene Ratio")+ylab("Gene Set")+theme_bw()
  print(p)
  return(p)
}

GSEA2dp<-function(x,clusterType="Empty",showTerms,col_low="#1e6091",col_high="#99d98c",...){
  require(forcats)
  clusterType <- match.arg(clusterType, c("MatrixResult", "NetworkResult","Empty"))
  orderBy<-"NES"
  if(sum(c("GSCluster","NetworkCluster")%in%colnames(x@Results))==0|clusterType=="Empty"){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,NES=x@Results$NES,Cluster=NA)
  }
  if("GSCluster"%in%colnames(x@Results)&clusterType=="MatrixResult"){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,NES=x@Results$NES,Cluster=x@Results$GSCluster)
    FigDF$Cluster[is.na(FigDF$Cluster)]<-"UnKnown"
  }
  if("NetworkCluster"%in%colnames(x@Results)&clusterType=="NetworkResult"){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,NES=x@Results$NES,Cluster=x@Results$NetworkCluster)
    FigDF$Cluster[is.na(FigDF$Cluster)]<-"UnKnown"
  }
  if(nrow(FigDF)<2){
    stop("Enrichment result is too few!\n")
  }
  FigDF<-FigDF%>%dplyr::arrange(desc(!!sym(orderBy)))
  if(!(length(showTerms)>1)){
    stop("In GSEA, please pass a vector to showTerms to determine which items will show. \n")
  }
  keep<-c(1:nrow(FigDF))[c(1:nrow(FigDF))%in%showTerms]
  if(length(keep)==0){
    stop("NO overlap!\n")
  }
  p<-ggplot(FigDF[keep,],aes(x=NES, y=fct_reorder(Terms,NES),fill=pValue))+geom_col()+scale_fill_continuous(low=col_low, high=col_high, name = "adjustedP",guide=guide_colorbar(reverse=F))+xlab("NES")+ylab("Gene Set")+theme_bw()
  print(p)
  return(p)
}
