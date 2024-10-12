# EnrichGT

![](man/figures/EnrichGT_2.png)

### Installation

```{r}
library(devtools)
install_github("ZhimingYe/EnrichGT")
```

Description and help pages: <https://zhimingye.github.io/EnrichGT/>

### OverView

The main purpose of developing this package is to provide a **lightweight and practical solution** to the problems mentioned above. Specifically, this package aims to:

1.  **Cluster enrichment results** based on hit genes or core enrichment from GSEA using term frequency analysis (from the output of the powerful `clusterProfiler`). This provides a clearer view of biological relevance by focusing on the genes that matter most.

2.  **Light weight and ultra-fast**. See the [benchmark](https://zhimingye.github.io/EnrichGT/articles/c_benchmark.html) page for detail.

3.  **Focus on wet-labs**. Using `GOSemSim` to measure the semantic similarity between pathways is, theoretically, the best way to tackle redundancy. However, in practical cases—especially in experimental bioinformatics validation—researchers are more focused on the genes behind these pathways rather than the pathways themselves.

4.  **Visualize results** using the **Posit PBC's `gt` package**, which offers two key advantages:

    -   **Improved readability** compared to traditional tables, making it easier to interpret results.
    -   **HTML-based output** that allows for internal searching and filtering, making exploration more user-friendly.

5.  **Seamless Integration with Publishing Tools** : Another significant advantage of using the `gt` package is its seamless integration with modern publishing systems like **rmarkdown** and **quarto**. This ensures that your results can be easily shared and published in a professional format, while also allowing for interactive exploration during the research phase.

6.  **Infering Pathway and TF Activity**: derived from clustering-based meta-gene modules, infer regulatory transcription factors and pathway activity.

### Acknowledgement

`EnrichGT` is based on the enrichment result from the powerful `clusterProfiler`. When using this, you should cite `clusterProfiler` and `ReactomePA`:

1.  T Wu<sup>\#</sup>, E Hu<sup>\#</sup>, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo<sup>\*</sup>, **G Yu**<sup>\*</sup>. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. ***The Innovation***. 2021, 2(3):100141. doi: [10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141)

2.  **G Yu**, QY He^\*^. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. ***Molecular BioSystems***. 2016, 12(2):477-479. doi: [10.1039/C5MB00663E](https://doi.org/10.1039/C5MB00663E)
