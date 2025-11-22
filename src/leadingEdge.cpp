#include <Rcpp.h>
using namespace Rcpp;

// stats Ranked Num vec
// hitIndices: hitted 1 based
// exponent: gseaParam
//
// [[Rcpp::export]]
List gseaScores2(NumericVector stats,
                         IntegerVector hitIndices,
                         double exponent = 1.0) {

  int N = stats.size();
  LogicalVector hits(N);
  int Nh = 0;

  for (int i = 0; i < hitIndices.size(); ++i) {
    int idx = hitIndices[i] - 1;
    if (idx >= 0 && idx < N && !hits[idx]) {
      hits[idx] = true;
      ++Nh;
    }
  }
  NumericVector runningES(N);

  if (Nh == 0) {
    IntegerVector position(N);
    IntegerVector x = seq_len(N);
    DataFrame df = DataFrame::create(
      _["x"]           = x,
      _["runningScore"] = runningES,
      _["position"]     = position
    );
    return List::create(
      _["ES"]        = 0.0,
      _["runningES"] = df
    );
  }

  double NR = 0.0;
  for (int i = 0; i < N; ++i) {
    if (hits[i]) {
      NR += std::pow(std::fabs(stats[i]), exponent);
    }
  }

  double missIncrement = 1.0 / static_cast<double>(N - Nh);
  double cumHit = 0.0, cumMiss = 0.0;
  double maxES = -1e30, minES =  1e30;
  for (int i = 0; i < N; ++i) {
    if (hits[i]) {
      cumHit += std::pow(std::fabs(stats[i]), exponent) / NR;
    } else {
      cumMiss += missIncrement;
    }
    runningES[i] = cumHit - cumMiss;
    if (runningES[i] > maxES) maxES = runningES[i];
    if (runningES[i] < minES) minES = runningES[i];
  }

  double ES = (std::fabs(maxES) > std::fabs(minES)) ? maxES : minES;
  IntegerVector position(N);
  for (int i = 0; i < N; ++i) {
    position[i] = hits[i] ? 1 : 0;
  }
  IntegerVector x = seq_len(N);
  DataFrame df = DataFrame::create(
    _["x"]            = x,
    _["runningScore"] = runningES,
    _["position"]     = position
  );
  return List::create(
    _["ES"]        = ES,
    _["runningES"] = df
  );
}
