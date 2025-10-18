

# EnrichGT 


**EnrichGT \<-** Fast, light weight enrichment analysis + insightful re-clustering results make all results explainable + Pretty HTML tables, Just in **ONE** package, designed for researchers in wet-labs. Supported databases including GO, KEGG, Reactome, MsigDB + AI based (LLM) result annotations and more ... 


You can use ask AI for help in [![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/ZhimingYe/EnrichGT) . 

Please see the package website for more info: <https://zhimingye.github.io/EnrichGT/>. 


- Efficient functions for rapid enrichment analysis

- Re-clustering of enriched results provides clear and actionable insights

- User-friendly HTML output that is easy to read and interpret

- AI based (LLM) result annotations

- Do a series of things just in ONE package



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
        P2[egt_plot_gsea]
    end

    subgraph egt_recluster_analysis
        K1[Pretty table]
        K2[AI Annotation]
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

    CC --> K2
    MG --> K2

    CC --> P1
    C --> P2

    MG --> I


```
