library(Rcpp)

sourceCpp("../src/ora.cpp")

background_genes <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF", "GeneG", "GeneH", "GeneI", "GeneJ")
input_genes <- c("GeneA", "GeneC", "GeneF", "GeneH")
gene_sets <- list(
  c("GeneA", "GeneB", "GeneC"),
  c("GeneD", "GeneE", "GeneF", "GeneG"),
  c("GeneH", "GeneI", "GeneJ"),
  c("GeneZ", "GeneX", "GeneI")
)

result <- ora_analysis(gene_sets, input_genes, background_genes)
print(result)
