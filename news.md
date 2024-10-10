## EnrichGT 0.2.7

### Package major update

-   Add `CompareGT()` function (beta)

## EnrichGT 0.2.5

### Package major update

-   Bug fix

-   `EnrichGT()` supports a list containing multiple source enrichment results (e.g., GO and KEGG enrichment for a same DEG list), and fusing them. User can gain [Metascape like multi-DB enrichment](https://metascape.org/) result by using the fusing function (See Articles).

## EnrichGT 0.2.0

### Package major update

-   Bug fix in generating results.

-   Add `EnrichGT_obj` S4 object. It not only returns `gt` HTML tables, users can extract enriched tables, meta-gene modules, pathway modules in it.

### Others

-   Split `main.R` to multiple files, for better managing of codes.

-   Add package vignettes.

## EnrichGT 0.1.7

### Package major update

-   add `doGO()`, `doKEGG()`, `doRA()` wrapper for ORA enrichment.

-   add `BuildMultigroupDEGlist()` and `Ranked.GS()` for simplifying procedure of multi-group enrichment and `GSEA`.

## EnrichGT 0.1.2

### Package major update

-   Support GSEA results.

-   Give up `ggplot2` figures in `gt` object for speed up. Though the result is much ugly than before :)

-   Intelligent sensing the cluster numbers.

-   Provide self-check of column names, if `setReadable()` before, etc.

### Others

-   Auto check package dependencies when start-up.

## EnrichGT 0.1.0

-   Initial codes for generating `gt` objects for ORA objects.

-   Initial codes for enrichment, use `text2vec` package to parse core enriched genes.
