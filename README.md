

# EnrichGT <a href="https://zhimingye.github.io/EnrichGT/"><img src="https://zhimingye.github.io/EnrichGT/EnrichGTlogo.png" align="right" height="138" alt="EnrichGT website" /></a>

**EnrichGT \<-** Fast, light weight enrichment analysis + insightful re-clustering results make all results explainable + Pretty HTML tables, Just in **ONE** package, designed for researchers in wet-labs. Supported databases including GO, KEGG, Reactome, MsigDB + AI based (LLM) result annotations and more ... 

![](https://zhimingye.github.io/EnrichGT/Introduction.jpeg)

Please see the package website for more info: <https://zhimingye.github.io/EnrichGT/>

- Efficient C++-based functions for rapid enrichment analysis

- Simple input format, empowering non-pro users

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