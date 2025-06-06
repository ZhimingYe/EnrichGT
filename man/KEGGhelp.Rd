% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kegg.R
\name{database_KEGG}
\alias{database_KEGG}
\alias{database_KEGG_show_organism}
\title{Get KEGG database from KEGG website}
\usage{
database_KEGG(kegg_organism="hsa",OrgDB = org.Hs.eg.db,kegg_modules=F,local_cache=F)

database_KEGG_show_organism()
}
\arguments{
\item{kegg_organism}{Determine which species data from KEGG will be fetch. For human, it would be `hsa`(in default); For mouse, it would be `mmu`. If you wants other species, see `database_kegg_show_organism()` for details.}

\item{OrgDB}{The AnnotationDbi database to convert KEGG gene ID to gene symbols. For human it would be `org.Hs.eg.db`, for mouse it would be `org.Mm.eg.db`. In AnnotationDbi there are many species, please search `AnnotationDbi` for other species annotation database.}

\item{kegg_modules}{If TRUE, returns KEGG module; If FALSE returns KEGG pathways. In default, this is setted to FALSE to get mouse commonly used KEGG pathways.}

\item{local_cache}{cache a copy in local working folder. It will be saved as a `.enrichgt_cache` file in working dictionary. The `.enrichgt_cache` is just a `.rds` file, feel free to read it using `readRDS()`.}
}
\value{
data.frame contains KEGG annotations
}
\description{
KEGG is a commercialized database. So EnrichGT can't pre-cache them locally. You can use this function to fetch KEGG database pathways and modules.
}
