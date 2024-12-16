#include <Rcpp.h>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <unordered_set>

using namespace Rcpp;

// Calculation of logarithmic factorial
double log_factorial(int n) {
  static std::unordered_map<int, double> cache; // Cache results
  if (n <= 1) return 0.0;
  if (cache.find(n) != cache.end()) return cache[n];

  double result = log_factorial(n - 1) + std::log(n);
  cache[n] = result;
  return result;
}

// Logarithmic combination calculation
double log_combination(int n, int k) {
  if (k > n || k < 0) return -std::numeric_limits<double>::infinity();
  return log_factorial(n) - log_factorial(k) - log_factorial(n - k);
}

// Hypergeometric distribution P(X >= k)
double hypergeometric_cdf_upper(int k, int N, int K, int n) {
  double log_p = -std::numeric_limits<double>::infinity();
  for (int i = k; i <= std::min(K, n); ++i) {
    double log_term = log_combination(K, i) + log_combination(N - K, n - i) - log_combination(N, n);
    log_p = std::log(std::exp(log_p) + std::exp(log_term));
  }
  return std::exp(log_p);
}

// [[Rcpp::export]]
DataFrame ora_analysis(
    std::vector<std::vector<std::string>> gene_sets, // List of gene sets
    std::vector<std::string> input_genes,           // Input gene list
    std::vector<std::string> background_genes       // Background gene set
) {
  int N = background_genes.size();  // Number of background genes
  int n = input_genes.size();       // Number of input genes

  // Build a hash table for the background gene set
  std::unordered_set<std::string> background_set(background_genes.begin(), background_genes.end());

  // Build a hash table for the input genes (for faster lookup)
  std::unordered_set<std::string> input_set(input_genes.begin(), input_genes.end());

  std::vector<int> overlaps(gene_sets.size());
  std::vector<double> p_values(gene_sets.size());

  // Iterate over each gene set
  for (size_t i = 0; i < gene_sets.size(); ++i) {
    const auto& gene_set = gene_sets[i];

    // Build a hash table for the current gene set
    std::unordered_set<std::string> gene_set_set(gene_set.begin(), gene_set.end());

    // Calculate the overlap between the input genes and the gene set
    int K = 0;
    int k = 0;
    for (const auto& gene : gene_set) {
      if (background_set.find(gene) != background_set.end()) {
        ++K; // Number of background genes in the gene set
        if (input_set.find(gene) != input_set.end()) {
          ++k; // Number of overlapping genes
        }
      }
    }

    // Skip calculation if K = 0
    if (K == 0) {
      overlaps[i] = 0;
      p_values[i] = 1.0; // Not significant
      continue;
    }

    // Calculate the hypergeometric distribution P-value
    overlaps[i] = k;
    p_values[i] = hypergeometric_cdf_upper(k, N, K, n);
  }

  // Return results as a DataFrame
  return DataFrame::create(
    Named("GeneSetIndex") = seq(1, gene_sets.size()),
    Named("Overlap") = overlaps,
    Named("PValue") = p_values
  );
}
