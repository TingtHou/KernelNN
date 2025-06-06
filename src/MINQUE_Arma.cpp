#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List MINQUE0(Rcpp::List KList, arma::vec y) {
  int VCs = KList.size();
  arma::mat C_mat(VCs, VCs, arma::fill::zeros);
  arma::vec RightY(VCs, arma::fill::zeros);

  // Compute C_mat and RightY efficiently
  for (int i = 0; i < VCs; i++) {
    arma::mat Ki = Rcpp::as<arma::mat>(KList[i]);
    RightY(i) = as_scalar(y.t() * Ki * y);

    for (int j = i; j < VCs; j++) {
      arma::mat Kj = Rcpp::as<arma::mat>(KList[j]);
      double val = arma::accu(Ki % Kj);  // Element-wise multiplication and sum
      C_mat(i, j) = val;
      C_mat(j, i) = val;
    }
  }

  // Solve C_mat * est_vcs = RightY
  arma::vec est_vcs = arma::solve(C_mat, RightY);

  return Rcpp::List::create(Rcpp::Named("vcs") = est_vcs);
}
// [[Rcpp::export]]
Rcpp::List MINQUE(Rcpp::List KList, arma::vec y, arma::vec prior) {
  int VCs = KList.size();
  // Convert R list of matrices into a std::vector<arma::mat>
  std::vector<arma::mat> K;
  K.reserve(VCs);
  for (int i = 0; i < VCs; ++i) {
    arma::mat Ki = Rcpp::as<arma::mat>(KList[i]);
    K.push_back(Ki);
  }

  int n = y.n_rows;
  // Build V = sum_{i=1}^VCs prior[i] * K[i]
  arma::mat V = arma::zeros<arma::mat>(n, n);
  for (int i = 0; i < VCs; ++i) {
    V += prior(i) * K[i];
  }

  // Invert V (with check)
  arma::mat Q;
  bool okV = arma::inv(Q, V);
  if (!okV || !Q.is_finite()) {
    return Rcpp::List::create(Rcpp::Named("err") = -1);
  }

  arma::vec Qy = Q * y;
  // Precompute K[i] %*% Q for each i
  std::vector<arma::mat> KQ(VCs);
  for (int i = 0; i < VCs; ++i) {
    KQ[i] = K[i] * Q;
  }

  // Build C.mat and RightY
  arma::mat Cmat = arma::zeros<arma::mat>(VCs, VCs);
  arma::vec RightY = arma::zeros<arma::vec>(VCs);

  for (int i = 0; i < VCs; ++i) {
    for (int j = 0; j < VCs; ++j) {
      // element-wise multiplication then sum
      Cmat(i, j) = accu(KQ[i].t() % KQ[j]);
    }
    // right-hand side: t(Qy) %*% K[i] %*% Qy
    RightY(i) = as_scalar(Qy.t() * K[i] * Qy);
  }

  // Invert C.mat
  arma::mat Cinv;
  bool okC = arma::inv(Cinv, Cmat);
  if (!okC || !Cinv.is_finite()) {
    return Rcpp::List::create(Rcpp::Named("err") = -1);
  }

  arma::vec est_vcs = Cinv * RightY;
  return Rcpp::List::create(Rcpp::Named("vcs") = est_vcs);
}
// [[Rcpp::export]]
Rcpp::List IMINQUE(Rcpp::List KList,
                 arma::vec y,
                 arma::vec prior,
                 int epoch = 100,
                 double threshold = 1e-4,
                 bool echo=false) {
  int VCs = KList.size();
  if (prior.is_empty()) {
    // initialize prior to a vector of 1s if not provided
    prior = arma::ones<arma::vec>(VCs);
  }
  arma::vec prior_est = prior;
  int it = 1;

  while (it < epoch) {
    // Ensure non-negative
    // prior_est = arma::abs(prior_est);

    // Call MINQUE_cpp
    Rcpp::List est = MINQUE(KList, y, prior_est);
    if (est.containsElementNamed("err")) {
      return Rcpp::List::create(Rcpp::Named("err") = -1);
    }

    arma::vec est_vcs = Rcpp::as<arma::vec>(est["vcs"]);
    double diff = norm(est_vcs - prior_est, 2);
    //
    // // Print iteration info
    if(echo)
    {
      Rcpp::Rcout << "IT: " << it << ", est: ";
      for (int i = 0; i < VCs; ++i) {
        Rcpp::Rcout << est_vcs(i);
         if (i < VCs - 1) Rcpp::Rcout << ",";
      }
      Rcpp::Rcout << ", diff: " << diff << std::endl;
    }

    if (diff < threshold) {
      break;
    } else {
      prior_est = est_vcs;
      ++it;
    }
  }
  return Rcpp::List::create(Rcpp::Named("vcs") = prior_est,Rcpp::Named("it") = it);
}
