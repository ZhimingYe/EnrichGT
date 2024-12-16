#include <Rcpp.h>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <unordered_set>

using namespace Rcpp;

// 对数阶乘的计算
double log_factorial(int n) {
  static std::unordered_map<int, double> cache; // 缓存结果
  if (n <= 1) return 0.0;
  if (cache.find(n) != cache.end()) return cache[n];

  double result = log_factorial(n - 1) + std::log(n);
  cache[n] = result;
  return result;
}

// 对数组合数
double log_combination(int n, int k) {
  if (k > n || k < 0) return -std::numeric_limits<double>::infinity();
  return log_factorial(n) - log_factorial(k) - log_factorial(n - k);
}

// 超几何分布的 P(X >= k)
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
    std::vector<std::vector<std::string>> gene_sets, // 基因集列表
    std::vector<std::string> input_genes,           // 输入基因列表
    std::vector<std::string> background_genes       // 背景基因集
) {
  int N = background_genes.size();  // 背景基因数
  int n = input_genes.size();       // 输入基因数

  // 构建背景基因集的哈希表
  std::unordered_set<std::string> background_set(background_genes.begin(), background_genes.end());

  // 构建输入基因的哈希表（加速查找）
  std::unordered_set<std::string> input_set(input_genes.begin(), input_genes.end());

  std::vector<int> overlaps(gene_sets.size());
  std::vector<double> p_values(gene_sets.size());

  // 遍历每个基因集
  for (size_t i = 0; i < gene_sets.size(); ++i) {
    const auto& gene_set = gene_sets[i];

    // 构建当前基因集的哈希表
    std::unordered_set<std::string> gene_set_set(gene_set.begin(), gene_set.end());

    // 计算输入基因与基因集的重叠数
    int K = 0;
    int k = 0;
    for (const auto& gene : gene_set) {
      if (background_set.find(gene) != background_set.end()) {
        ++K; // 基因集的背景基因数量
        if (input_set.find(gene) != input_set.end()) {
          ++k; // 重叠基因数量
        }
      }
    }

    // 如果 K = 0，跳过计算
    if (K == 0) {
      overlaps[i] = 0;
      p_values[i] = 1.0; // 不显著
      continue;
    }

    // 计算超几何分布的 P 值
    overlaps[i] = k;
    p_values[i] = hypergeometric_cdf_upper(k, N, K, n);
  }

  // 返回结果为一个 DataFrame
  return DataFrame::create(
    Named("GeneSetIndex") = seq(1, gene_sets.size()),
    Named("Overlap") = overlaps,
    Named("PValue") = p_values
  );
}
