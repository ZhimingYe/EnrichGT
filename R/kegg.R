#


database_kegg_show_organism <- function(){
  x <- read.delim("https://rest.kegg.jp/list/organism",quote = "\t",header = F)
  return(x)
}

keggModuleList <- function(orgkegg){
  k0<-glue::glue("https://rest.kegg.jp/list/pathway/{orgkegg}")
  k1<-glue::glue("https://rest.kegg.jp/link/pathway/{orgkegg}")
  k2<-glue::glue("https://rest.kegg.jp/conv/{orgkegg}/ncbi-geneid")
  m0<-glue::glue("https://rest.kegg.jp/list/module")
  m1<-glue::glue("https://rest.kegg.jp/link/module/{orgkegg}")

  return(c(k0,k1,k2,m0,m1))

}
# y <- read.delim("https://rest.kegg.jp/link/module/mmu",quote = "\t",header = F)
database_kegg <- function(kegg_organism="hsa",OrgDB = org.Hs.eg.db,kegg_modules=F,local_cache=F){
  if(!kegg_modules){
    fn <- paste0("KEGGPathway_",xfun::md5(keggModuleList(kegg_organism)))
  }else{
    fn <- paste0("KEGGModules_",xfun::md5(keggModuleList(kegg_organism)))
  }
  if(local_cache){
    if(file.exists(paste0(fn,".enrichgt_cache"))){
      tryCatch({
        cachedFile<-readRDS(paste0(fn,".enrichgt_cache"))
      },error=function(e){
        cli::cli_abort("Load local cache error! Please re-check or set local_cache=F to using online files")
      })

      cli::cli_alert_success("Found on disk cache file and loaded. ")
      assign(fn,cachedFile,envir = db_getter_env)
    }
  }
  if(exists(fn,envir = db_getter_env)){
    finalDF <- get(fn,envir = db_getter_env)
    cli::cli_alert_info(paste0("Use cached KEGG downloaded files: ",fn))
  }else{
    if(!kegg_modules){
      targetsind<-1:3
    }else{
      targetsind<-c(4,5,3)
    }
    keggs <- lapply(keggModuleList(kegg_organism)[targetsind],function(x){
      tryCatch({
        cli::cli_alert_info(glue::glue("downloading {x}..."))
        y <- read.delim(x,quote = "\t",header = F)
        Sys.sleep(2)
        return(y)
      },error=function(e){
        cli::cli_abort("Please re-check your kegg_organism input and OrgDB input. Or please check your network(Because of the limitation of KEGG license, you must download them instead of built in in EnrichGT). ")
      })

    })
    colnames(keggs[[1]])<-c("Pws","Terms")
    colnames(keggs[[2]])<-c("Kegggenes","Pws")
    colnames(keggs[[3]])<-c("Ncbigenes","Kegggenes")
    if(!kegg_modules){
      keggs[[2]]$Pws<-gsub("path:","",keggs[[2]]$Pws)
    }else{
      keggs[[2]]$Pws<-s_(keggs[[2]]$Pws,"_",2)
    }

    finalDF<-keggs[[2]] |> left_join(keggs[[1]],by="Pws",relationship = "many-to-many") |> left_join(keggs[[3]],by="Kegggenes",relationship = "many-to-many")

    assign(fn,finalDF,envir = db_getter_env)
  }
  if(!local_cache){
    finalDF$Ncbigenes<-gsub("ncbi-geneid:","",finalDF$Ncbigenes)
    cvtDF<-convert_annotations_genes(finalDF$Ncbigenes,"ENTREZID","SYMBOL",OrgDB = OrgDB)
    finalDF<-finalDF |> dplyr::rename(ENTREZID = Ncbigenes)
    finalDF<-finalDF |> left_join(cvtDF,by="ENTREZID",relationship = "many-to-many") |> dplyr::select(Pws,Terms,SYMBOL)
    finalDF$ALL <- paste0(finalDF$Pws,finalDF$Terms,finalDF$SYMBOL)
    finalDF<-finalDF[!duplicated(finalDF$ALL),-which(colnames(finalDF)=="ALL")]

    finalDF$Terms<-s_(finalDF$Terms," - ",1)
  }

  res0 <- finalDF
  if(local_cache){
    saveRDS(res0,(paste0(fn,".enrichgt_cache")))
    cli::cli_alert_success(paste0("Wrote KEGG cache on disk: ",fn,".enrichgt_cache"))
  }
  return(res0)

}

# library(org.Hs.eg.db)
# ddd<-database_kegg_pathway(kegg_modules = F,local_cache = T)
