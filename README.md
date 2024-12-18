# EnrichGT 0.5

**EnrichGT \<-** Fast, light weight enrichment analysis + explainable re-clustered HTML tables

![](https://zhimingye.github.io/EnrichGT/enrichGTTable.jpg)

``` mermaid

graph LR
    subgraph Enrichment Analysis
        A[egt_enrichment_analysis]
        B[egt_gsea_analysis]
    end

    subgraph Pathway Databases
        D[database_* funcs]
    end

    subgraph Visualize results
        P1[egt_plot_results]
        P2[egt_plot_umap]
    end

    subgraph egt_recluster_analysis
        K1[Pretty table]
        CC[cluster modules]
        MG[gene modules]
    end

    subgraph Pathway Act. and TF infer 
        
        I[egt_infer]
    end

    D --> A
    D --> B

    A --> C[Enriched Result]
    B --> C

    C --> CC
    C --> MG

    C --> P1

    CC --> K1
    MG --> K1

    CC --> P1
    CC --> P2

    MG --> I

```

Description and help pages: <https://zhimingye.github.io/EnrichGT/>

# Install EnrichGT

``` r
install.packages("pak")
pak::pkg_install("ZhimingYe/EnrichGT")
```

or

``` r
install.packages("devtools")
library(devtools)
install_github("ZhimingYe/EnrichGT")
```
