#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List MNQTest0_Overall(
  Rcpp::List KList,
  arma::vec vcs,
//  Rcpp::IntegerVector index_interest,
  arma::vec wgt) {

  int nVCs = vcs.size();
  int Nn = Rcpp::as<arma::mat>(KList[0]).n_rows;

  // Compute C_mat
  arma::mat C_mat(nVCs, nVCs, arma::fill::zeros);
  for (int i = 0; i < nVCs; ++i) {
    arma::mat Ki = Rcpp::as<arma::mat>(KList[i]);
    for (int j = i; j < nVCs; ++j) {
      arma::mat Kj = Rcpp::as<arma::mat>(KList[j]);
      double val = arma::accu(Ki % Kj);
      C_mat(i, j) = val;
      C_mat(j, i) = val;
    }
  }

  // Invert C_mat
  arma::mat C_mat_inv;
  bool success = arma::inv(C_mat_inv, C_mat);
  if (!success) {
    Rcpp::Rcout << "C_mat is not invertible:\n" << C_mat;
    return Rcpp::List::create(Rcpp::Named("err") = -1);
  }

  // Compute QList
  std::vector<arma::mat> QList(nVCs, arma::zeros<arma::mat>(Nn, Nn));
  for (int i = 0; i < nVCs; ++i) {
    for (int j = 0; j < nVCs; ++j) {
      QList[i] += C_mat_inv(i, j) * Rcpp::as<arma::mat>(KList[j]);
    }
  }

  // Compute ratio and sigma_total
  double ratio = 0.0;
  arma::mat sigma_total = arma::zeros<arma::mat>(Nn, Nn);
  for (int i = 0; i < nVCs-1; ++i) {
    ratio += vcs(i) * wgt(i);
    sigma_total += QList[i] * wgt(i);
  }
//  for (int idx : index_interest) {
//    idx=idx-1;
    // Rcpp::Rcout<<idx<<"\n";
//    ratio += vcs(idx) * wgt(idx);
//    sigma_total += QList[idx] * wgt(idx);
//  }

  ratio /= (wgt(nVCs - 1) * vcs(nVCs - 1));
  sigma_total -= ratio * QList[nVCs - 1] * wgt[nVCs - 1];

  // Eigenvalues
  arma::vec eigen_test = arma::eig_sym(sigma_total);

  // Call R's davies function
  Rcpp::Function davies("davies", Rcpp::Environment::namespace_env("CompQuadForm"));
  Rcpp::List davies_result = davies(0, Rcpp::Named("lambda", eigen_test));
  double pvalue = davies_result["Qq"];

  return Rcpp::List::create(
    Rcpp::Named("z") = ratio,
    Rcpp::Named("pvalue") = pvalue
  );
}

// [[Rcpp::export]]
Rcpp::List MNQTest0_Component(Rcpp::List KList,
                                              arma::vec vcs,
                                              arma::vec vcs_h0,
                                              Rcpp::IntegerVector index_interest,
                                              arma::vec wgt) {
  int nVCs = KList.size();
  int Nn = Rcpp::as<arma::mat>(KList[0]).n_rows;
  arma::mat Cmat(nVCs, nVCs, arma::fill::zeros);

  // Fill C matrix
  for (int i = 0; i < nVCs; ++i) {
    arma::mat Ki =  Rcpp::as<arma::mat>(KList[i]);
    for (int j = i; j < nVCs; ++j) {
      arma::mat Kj =  Rcpp::as<arma::mat>(KList[j]);
      double val = accu(Ki % Kj);
      Cmat(i, j) = val;
      Cmat(j, i) = val;
    }
  }

  arma::mat Cinv;
  bool success = inv(Cinv, Cmat);
  if (!success || Cinv.has_nan() || any(diagvec(Cinv) < 0)) {
    return Rcpp::List::create(Rcpp::Named("err") = -1);
  }

  std::vector<arma::mat> QList(nVCs, arma::zeros<arma::mat>(Nn, Nn));
  for (int i = 0; i < nVCs; ++i) {
    for (int j = 0; j < nVCs; ++j) {
      QList[i] += Cinv(i, j) * Rcpp::as<arma::mat>(KList[j]);
    }
  }

  // Compute sigma_total from multiple indices
  double ratio = 0.0;
  arma::mat sigma_total = arma::zeros<arma::mat>(Nn, Nn);
  std::set<int> test_index_set;

  for (int k = 0; k < index_interest.size(); ++k) {
    int idx = index_interest[k] - 1; // R to C++ index
    test_index_set.insert(idx);
    sigma_total += wgt[idx] * QList[idx];
    ratio += wgt[idx] * vcs[idx];
  }


  arma::mat V_sigma_hat_h0 = arma::zeros<arma::mat>(Nn, Nn);
  int h0_idx=0;
  for (int i = 0; i < nVCs; ++i) {
    if (test_index_set.count(i) == 0) {
      V_sigma_hat_h0 += vcs_h0[h0_idx++] *  Rcpp::as<arma::mat>(KList[i]);
    }
  }

  arma::mat M=sigma_total*V_sigma_hat_h0;
  // Eigenvalues of M
  arma::vec eigen_test=arma::real(arma::eig_gen(M));
  // Compute p-value using Davies method
  Rcpp::Function davies("davies", Rcpp::Environment::namespace_env("CompQuadForm"));
  Rcpp::List davies_result = davies(ratio, Rcpp::Named("lambda", eigen_test));
  double pvalue = davies_result["Qq"];

  return Rcpp::List::create(
    Rcpp::Named("pvalue") = pvalue
  );
}


// [[Rcpp::export]]
Rcpp::List IMNQTest_Normal(Rcpp::List KList,
                           arma::vec vcs,
                           Rcpp::IntegerVector index_interest)
{
  int nVCs = vcs.size();
  int nN = Rcpp::as<arma::mat>(KList[0]).n_rows;

  // Calculate V_sigma = sum_i KList[[i]] * est_vcs[i]
  arma::mat V_sigma = arma::zeros(nN, nN);
  for (int i = 0; i < nVCs; ++i) {
    V_sigma += Rcpp::as<arma::mat>(KList[i]) * vcs[i];
  }

  arma::mat V_sigma_inv;
  bool success = arma::inv(V_sigma_inv, V_sigma);
  if (!success || V_sigma_inv.has_nan()) {
    return Rcpp::List::create(Rcpp::Named("err") = -1);
  }

  // Compute SigmaVi
  std::vector<arma::mat> SigmaVi(nVCs);
  for (int i = 0; i < nVCs; ++i) {
    SigmaVi[i] = V_sigma_inv * Rcpp::as<arma::mat>(KList[i]);
  }

  // Compute C.mat
  arma::mat C_mat(nVCs, nVCs, arma::fill::zeros);
  for (int i = 0; i < nVCs; ++i) {
    for (int j = 0; j < nVCs; ++j) {
      C_mat(i, j) = arma::accu(SigmaVi[i].t() % SigmaVi[j]) * 0.5;
    }
  }

  arma::mat C_mat_inv;
  success = arma::inv(C_mat_inv, C_mat);
  if (!success || C_mat_inv.has_nan()) {
    return Rcpp::List::create(Rcpp::Named("err") = -1);
  }

  // Standardize est_vcs[test_index]
  int test_len = index_interest.size();
  arma::uvec test_uvec(test_len);
  for (int k = 0; k < test_len; ++k) test_uvec[k] = index_interest[k] - 1; // Convert R 1-based to C++ 0-based

  arma::mat subC = C_mat_inv.submat(test_uvec, test_uvec);
  arma::vec vcs_sub = vcs.elem(test_uvec);

  // Matrix square root
  arma::mat sqrt_subC;
  bool sqrt_ok = arma::sqrtmat_sympd(sqrt_subC, subC);
  if (!sqrt_ok) return Rcpp::List::create(Rcpp::Named("err") = -1);

  arma::vec vcs_standard = vcs; // Copy original
  vcs_standard.elem(test_uvec) = arma::solve(sqrt_subC, vcs_sub);

  // Check for NaN or Inf
  if (vcs_standard.has_nonfinite()) {
    return Rcpp::List::create(Rcpp::Named("err") = -1);
  }

  // Chi-squared and zstats
  double chi = arma::dot(vcs_standard.elem(test_uvec), vcs_standard.elem(test_uvec));
  Rcpp::NumericVector zstat(nVCs), pvalue(nVCs);
  for (int i = 0; i < nVCs; ++i) {
    double wald = vcs[i] / std::sqrt(C_mat_inv(i, i));
    zstat[i] = wald;
    pvalue[i] = R::pnorm(wald, 0, 1, false, false);
  }

  double sum_standard = arma::accu(vcs_standard.elem(test_uvec));
  double chiP = (sum_standard > 0) ?
  R::pchisq(chi, test_len, false, false) / 2.0 : 1.0;

  return Rcpp::List::create(
    Rcpp::Named("z") = zstat,
    Rcpp::Named("components") = pvalue,
    Rcpp::Named("overall") = chiP,
    Rcpp::Named("df") = test_len
  );
}
