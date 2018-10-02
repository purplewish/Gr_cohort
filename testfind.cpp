#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::uvec find1(arma::vec indexy, int unid)
{

  uvec indexi = find(indexy == unid);
  return(indexi);
}