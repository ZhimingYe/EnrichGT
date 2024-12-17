#' Visualize results generated form `EnrichGT()` using simple plot
#'
#' @description
#' This plot is the most widely like `enrichplot::dotplot()`used method to visualize enriched terms. It shows the enrichment scores (e.g. p values) and gene ratio or NES as dot size and color / or bar height. Users can specify the number of terms using `ntop` or selected terms to color via the `low.col` and `hi.col`.
#' @param x an EnrichGT object
#' @param ntop Show top N in each cluster
#' @param low.col the color for the lowest
#' @param hi.col the color for the highest
#' @param max_len_descript the label format length, default as 40.
#' @param ... Other param
#'
#' @returns a ggplot2 object
#' @export
#' @importFrom ggplot2 fortify
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_gradient
#' @importFrom ggplot2 scale_color_continuous
#' @importFrom ggplot2 scale_fill_continuous
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 guide_colorbar
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_size
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 geom_col
#' @importFrom methods is
#' @importFrom forcats fct_reorder
#' @author Zhiming Ye
egt_plot_results <- function(x,...){
  if(sum(colnames(x@enriched_result)=="absNES")==0){
    figure0<-ORA2dp(x,...)
  }else{
    figure0<-GSEA2dp(x,...)
  }
  return(figure0)
}

shorten_labels_words <- function(label, max_length = 40) {
  sapply(label, function(l) {
    words <- unlist(strsplit(l, " "))
    cumulative_length <- cumsum(nchar(words) + 1)
    if (max(cumulative_length) <= max_length) {
      return(l)
    }
    cutoff <- max(which(cumulative_length <= max_length))
    paste(paste(words[1:cutoff], collapse = " "), "...")
  })
}

#' Visualize results generated form `EnrichGT()` using UMAP
#'
#' @description
#' A word frequency matrix represents the frequency of words or tokens across different documents or text samples. Each row corresponds to a document, and each column represents a word or token, with the cell values indicating the frequency of the respective word in that document.However, high-dimensional data like word frequency matrices can be challenging to interpret directly. To make such data easier to analyze, we can reduce its dimensionality and visualize the patterns or clusters in a 2D or 3D space. UMAP (Uniform Manifold Approximation and Projection) is a powerful, non-linear dimensionality reduction technique widely used for this purpose.
#'
#' @param x an EnrichGT object
#' @param ... Other param
#'
#' @returns a ggplot2 object
#' @export
#'
#' @author Zhiming Ye
egt_plot_umap <- function(x,...){
  px<-.egtUMAP(x,...)
  return(px)
}

ORA2dp<-function(x,ntop = 7,low.col="#78cfe5",hi.col="#ff6f81",max_len_descript=40,...){
  if(is.list(x)){
    cli::cli_abort("For a list object, please run plotting for every object inside list, instead of the whole list.")
  }
  tryCatch({
    if(dim(x@enriched_result)[1]<2 | sum(colnames(x@enriched_result)=="Count")==0){
      cli::cli_abort("ERROR! ")
    }
  },error=function(e){
    cli::cli_abort("Not EnrichGT object! Please run `EnrichGT()` first.")
  })
  tryCatch({
    df <- x@enriched_result |>
      group_by(Cluster) |>
      slice_min(order_by = Count, n = ntop, with_ties = FALSE) |>
      ungroup()
  },error=function(e){
    cli::cli_alert_warning("Subset ERROR! ")
    df <-x@enriched_result
  })
  df$Description<-shorten_labels_words(df$Description,max_length = max_len_descript)
  px<-ggplot(df,aes(x = PCT, y = fct_reorder(Description, PCT), size=Count, color=Padj))+geom_point()+scale_color_continuous(low=low.col, high=hi.col, name = "adjustedP",guide=guide_colorbar(reverse=F))+scale_size(range=c(2, 8))+xlab("Gene Ratio")+ylab("Gene Sets")+facet_grid(Cluster~.,scales="free",space="free_y")+theme_bw()
  return(px)
}

GSEA2dp<-function(x,ntop = 7,low.col="#78cfe5",hi.col="#ff6f81",max_len_descript=40,...){
  if(is.list(x)){
    cli::cli_abort("For a list object, please run plotting for every object inside list, instead of the whole list.")
  }
  tryCatch({
    if(dim(x@enriched_result)[1]<2 | sum(colnames(x@enriched_result)=="absNES")==0){
      cli::cli_abort("ERROR! ")
    }
  },error=function(e){
    cli::cli_abort("Not EnrichGT object! Please run `EnrichGT()` first.")
  })
  tryCatch({
    df <- x@enriched_result |>
      group_by(Cluster) |>
      slice_min(order_by = absNES, n = ntop, with_ties = FALSE) |>
      ungroup()
  },error=function(e){
    cli::cli_alert_warning("Subset ERROR! ")
    df <-x@enriched_result
  })
  df$NES<-ifelse(df$Reg=="UpReg",df$absNES*(1),df$absNES*(-1))
  df$Description<-shorten_labels_words(df$Description,max_length = max_len_descript)
  px<-ggplot(df,aes(x = NES, y = fct_reorder(Description, absNES), fill=Padj))+geom_col()+scale_fill_continuous(low=low.col, high=hi.col, name = "adjustedP",guide=guide_colorbar(reverse=F))+scale_size(range=c(2, 8))+xlab("Normalize Enrichment Score(NES)")+ylab("Gene Sets")+facet_grid(Cluster~.,scales="free",space="free_y")+theme_bw()
  return(px)
}
