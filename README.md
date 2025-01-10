# EnrichGT 0.8

**EnrichGT \<-** Fast, light weight enrichment analysis + insightful re-clustering results make all results explainable + Pretty HTML tables, Just in **ONE** package, designed for researchers in wet-labs. 

Please see the package website for more info: <https://zhimingye.github.io/EnrichGT/>

> [!IMPORTANT]
> The primary goal of EnrichGT is to provide researchers in wet labs, who have been busy all day, with quick and insightful biological interpretations from dry lab data to support their experiments. Therefore, all computational methods employed are relatively straightforward and pragmatic. For example, C++ based ORA enrichment function, GSEA only using fgsea output without more analysis, use only vocabulary frequency matrix for re-enrichment instead of term's similarity, and more... Its purpose is to offer a “quick overview”. After this initial overview, you can use more widely recognized tools to generate statistically rigorous results (not necessary though). However, don’t misunderstand—this does not mean that the statistical processes in EnrichGT are incorrect. I’ve made every effort to ensure their accuracy; they are just less refined. Additionally, the tool has been tested in most typical scenarios, but extreme cases cannot be entirely ruled out.

![](https://zhimingye.github.io/EnrichGT/enrichGTTable.jpg)


### Install

``` r
install.packages("pak")
pak::pkg_install("ZhimingYe/EnrichGT")
```

### WorkFlows

``` mermaid

graph LR
    
    M[genes]
    N[genes with weights]
    
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
    
    M --> A
    N --> B
    
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

