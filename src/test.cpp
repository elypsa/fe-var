#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List sim_var(const arma::mat& X,
                   const arma::mat& bet_bs,
                   const arma::ivec& rr,
                   const arma::mat& resb,
                   const arma::ivec& fe_t,
                   int Ny,
                   int lags) {

  int n = X.n_rows;
  int p = X.n_cols;

  arma::mat Xb(n, Ny * lags, arma::fill::value(NA_REAL)); // no const
  arma::mat Yb(n, Ny, arma::fill::value(NA_REAL));

  int per = 0;
  for (int idx = 0; idx < fe_t.n_elem; ++idx) {
    int ind = fe_t(idx);
    Xb.row(per) = X.row(per).cols(0, Ny * lags - 1);

    for (int j = per; j < per + ind; ++j) {
      arma::rowvec yb_row = Xb.row(j) * bet_bs + resb.row(j);
      Yb.row(j) = yb_row;

      if (j < per + ind - 1) {
        // Xb[j + 1, ] <- c(Yb[j, ], Xb[j, 1:(Ny * (lags - 1))])
        arma::rowvec new_xb(Ny * lags);
        new_xb.cols(0, Ny - 1) = Yb.row(j);
        new_xb.cols(Ny, Ny * lags - 1) = Xb.row(j).cols(0, Ny * (lags - 1) - 1);
        Xb.row(j + 1) = new_xb;
      }
    }
    per += ind;
  }

  return Rcpp::List::create(
    Rcpp::Named("Xb") = Xb,
    Rcpp::Named("Yb") = Yb
  );
}

// for (ind in fe_t) {
//   Xb[per, ] <- X[per, 1:(Ny * lags)]
//   for (j in per:(per + ind - 1)) {
//     Yb[j, ] <- Xb[j, ] %*% bet_bs + resb[j, ]
//     if (j < (per + ind - 1)) {
//       Xb[j + 1, ] <- c(Yb[j, ], Xb[j, 1:(Ny * (lags - 1))])
//     }
//   }
//   per <- per + ind
// }
