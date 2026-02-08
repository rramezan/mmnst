#include <Rcpp.h>
using namespace Rcpp;

// Reverse comparison - this will result in the vector being sorted in reverse order
// i - 1/j - 1 are being used since C++ uses 0-indexing
bool comp_element (NumericVector x, int i, int j) { return (x[i - 1] >= x[j - 1]); }

// [[Rcpp::export]]
IntegerVector top_K_indices (NumericVector x, int K) {
  // Base conditions - return integer(0)
  if ((x.length() == 0) || (K == 0)) {
    return IntegerVector(0);
  }

  // Vector of the indices of x
  IntegerVector idx(seq(1, x.length()));


  auto lambda = [&x](int i, int j) {
      return comp_element(x, i, j);
  };

  // Sorts the first half of idx, so that the first K elements are the indices of
  //    the top K elements in x
  std::nth_element(idx.begin(), idx.begin() + K - 1, idx.end(), lambda);
  
  // Output only the first K elements
  return idx[Range(0, K - 1)];
}
