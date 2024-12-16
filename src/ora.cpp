#include <Rcpp.h>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <string>

using namespace Rcpp;

// 计算组合数 C(n, k)
long long combination(int n, int k) {
  if (k > n) return 0;
  if (k == 0 || k == n) return 1;
  long long result = 1;
  for (int i = 1; i <= k; ++i) {
    result *= (n - i + 1);
    result /= i;
  }
  return result;
}

// 超几何分布的 P(X >= k)
double hypergeometric_cdf_upper(int k, int N, int K, int n) {
  double p = 0.0;
  for (int i = k; i <= std::min(K, n); ++i) {
    p += (double)(combination(K, i) * combination(N - K, n - i)) / combination(N, n);
  }
  return p;
}

// [[Rcpp::export]]
DataFrame ora_analysis(
    std::vector<std::vector<std::string>> gene_sets, // 基因集列表
    std::vector<std::string> input_genes,           // 输入基因列表
    std::vector<std::string> background_genes       // 背景基因集
) {
  int N = background_genes.size();  // 总基因数
  int n = input_genes.size();       // 输入基因数

  std::vector<int> overlaps;
  std::vector<double> p_values;

  // 遍历每个基因集
  for (size_t i = 0; i < gene_sets.size(); ++i) {
    const auto& gene_set = gene_sets[i];
    int K = gene_set.size();  // 当前基因集中基因数

    // 计算输入基因列表中属于该基因集的基因数
    int k = 0;
    std::unordered_set<std::string> gene_set_set(gene_set.begin(), gene_set.end());
    for (const auto& gene : input_genes) {
      if (gene_set_set.find(gene) != gene_set_set.end()) {
        ++k;
      }
    }

    // 计算 P 值
    double p_value = hypergeometric_cdf_upper(k, N, K, n);

    overlaps.push_back(k);
    p_values.push_back(p_value);
  }

  // 返回结果为一个 DataFrame
  return DataFrame::create(
    Named("GeneSetIndex") = seq(1, gene_sets.size()),
    Named("Overlap") = overlaps,
    Named("PValue") = p_values
  );
}
