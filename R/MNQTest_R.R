#' Overall MINQUE0 Test (R implementation)
#'
#' R implementation of the MINQUE0-based overall variance component test.
#' It constructs a covariance kernel combination, inverts the MINQUE
#' coefficient matrix, and uses quadratic form theory with eigenvalues
#' to obtain a p-value via \code{CompQuadForm::davies()}.
#'
#' @param KList A list of \eqn{n \times n} kernel (covariance) matrices.
#'   All matrices must have the same dimension.
#' @param vcs A numeric vector of variance component estimates under the
#'   full model (one per kernel in \code{KList}).
#' @param index_interest Integer or numeric vector of indices (1-based)
#'   indicating which variance components are included in the overall test.
#' @param wgt A numeric vector of weights with length equal to
#'   \code{length(KList)}. Defaults to \code{rep(1, length(KList))}.
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item \code{z} The test statistic (ratio of weighted components relative
#'     to the last variance component).
#'   \item \code{pvalue} P-value from the Davies method applied to the
#'     eigenvalues of the adjusted covariance matrix.
#'   \item \code{err} If present with value \code{-1}, indicates that a matrix
#'     inversion failed or produced invalid values.
#' }
#'
#' @details
#' This function is the pure R counterpart to \code{MNQTest0_Overall_arm()},
#' and is used internally by \code{MNQTest0_ChiR()} for the overall test.
#'
#' @examples
#' \dontrun{
#' n <- 30
#' set.seed(1)
#' K1 <- diag(n)
#' K2 <- tcrossprod(matrix(rnorm(n * 3), n, 3))
#' vcs <- c(1, 0.5, 0.2)
#' MINQUE0_Overall_R(KList = list(K1, K2, diag(n)),
#'                   vcs = vcs,
#'                   index_interest = 1:2)
#' }
#'
#' @keywords internal
MINQUE0_Overall_R<-function(KList, vcs,index_interest,wgt=rep(1,length(KList)))
{
  # index_interest<-c(2,3)
  VCs<-length(vcs)
  Nn<-nrow(KList[[1]])
  C.mat<-matrix(NA,nrow=VCs,ncol=VCs)
  for (i in 1:VCs)
  {
    for (j in i:VCs)
    {
      C.mat[j,i]<-C.mat[i,j]<-sum(t( KList[[i]])* KList[[j]])
    }
  }
  C.mat<-tryCatch(solve(C.mat),
                  error = function(e){
                    print(e)
                    return(as.matrix(-999))
                  }
  )
  if (C.mat[1,1]==-999 || any(is.na(C.mat)))
  {
    return(list(err=-1))
  }
  if (any(diag(C.mat)<0)) {
    return(list(err=-1))
  }

  # sigmaH0   <- Reduce(f = "+", x = mapply(FUN = "*", KList,C.mat[VCs,], SIMPLIFY = FALSE)) # to calculate Psi=Phi(Sigma)
  QList<-list()
  for (i in 1:VCs)
  {
    sigmaH_i   <- Reduce(f = "+", x = mapply(FUN = "*", KList,C.mat[i,], SIMPLIFY = FALSE)) # to calculate Psi=Phi(Sigma)
    QList[[i]]<-sigmaH_i
  }

  ###################################################
  ####overall ###############################
  ######
  # wgt<-1/sqrt(diag(C.mat))
  ratio<-sum((vcs*wgt)[index_interest])
  sigma_total<-Reduce(f = "+", x = mapply(FUN = "*", QList[index_interest],wgt[index_interest], SIMPLIFY = FALSE))
  # for (i in 2:(VCs-1))
  # {
  #   ratio<-ratio+est_vcs[i]*wgt[j]
  #   sigma_total<-sigma_total+QList[[j]]*wgt[j]
  #   j<-j+1
  # }
  ratio<-ratio/(wgt[VCs]*vcs[VCs])
  sigma_total<-sigma_total-ratio*QList[[VCs]]*wgt[VCs]
  # eigen_test<-Get_Lambda(sigma_total)

  eigen_test<-eigen(sigma_total,symmetric=T,only.values = TRUE)$values
  # delta_nonCentered<-t(mu) %*% sigma_total %*% mu
  # delta_nonCentered<-ifelse(abs(delta_nonCentered)<1e-6,0,delta_nonCentered)
  # delta_nonCentered<-rep(delta_nonCentered,length(eigen_test))
  pvalue<-davies(q=0,lambda = eigen_test)$Qq
  return(list(z=ratio,pvalue=pvalue))
}



# MINQUE1_Overall_R<-function(KList, est_vcs, prior, index_interest=c(2,3),wgt=rep(1,length(KList)),covar=matrix(1,ncol =1, nrow = nrow(KList[[1]])))
# {
#   # index_interest<-c(2,3)
#
#   VCs<-length(est_vcs)
#   C.mat<-matrix(NA,nrow=VCs,ncol=VCs)
#   V.sigma   <- Reduce(f = "+", x = mapply(FUN = "*", KList,prior, SIMPLIFY = FALSE)) # to calculate Psi=Phi(Sigma)
#   V.sigma.inv<-solve(V.sigma)
#   Q<-V.sigma.inv-V.sigma.inv%*%covar%*%solve(t(covar)%*%V.sigma.inv%*%covar)%*%t(covar)%*%V.sigma.inv
#   # e<-eigen(KList[[1]])
#   # truncated<-max(e$values)
#   # a_point<--est_vcs[2]/truncated
#   # V.sigma<-solve(V.sigma)
#   SigmaVi<-list()
#   for (i in 1:VCs) {
#     SigmaVi[[i]]<-KList[[i]]%*%Q
#   }
#   for (i in 1:VCs)
#   {
#     for (j in 1:VCs)
#     {
#       C.mat[i,j]<-sum(t( SigmaVi[[i]])* SigmaVi[[j]])
#     }
#   }
#   C.mat<-tryCatch(solve(C.mat, tol = 1e-16),
#                   error = function(e){
#                     print(e)
#                     return(as.matrix(-999))
#                   }
#   )
#   if (C.mat[1,1]==-999 || any(is.na(C.mat)))
#   {
#     return(list(err=-1))
#   }
#   if (any(diag(C.mat)<0)) {
#     return(list(err=-1))
#   }
#
#   QList<-list()
#   for (i in 1:VCs)
#   {
#     sigmaH_i   <- Reduce(f = "+", x = mapply(FUN = "*", KList,C.mat[i,], SIMPLIFY = FALSE)) # to calculate Psi=Phi(Sigma)
#     QList[[i]]<-sigmaH_i
#   }
#
#   ###################################################
#   ####overall ###############################
#   ######
#
#   ratio<-sum((est_vcs*wgt)[index_interest])
#   sigma_total<-Reduce(f = "+", x = mapply(FUN = "*", QList[index_interest],wgt[index_interest], SIMPLIFY = FALSE))
#   ratio<-ratio/(wgt[VCs]*est_vcs[VCs])
#   sigma_total<-sigma_total-ratio*QList[[VCs]]*wgt[VCs]
#   # eigen_test<-Get_Lambda(sigma_total)
#   eigen_test<-eigen(sigma_total,symmetric=T,only.values = TRUE)$values
#   pvalue<-imhof(q=0,lambda = eigen_test)$Qq
#   if (est_vcs[VCs]<0) {
#     pvalue<-1-pvalue
#   }
#   return(list(z=ratio,pvalue=pvalue))
# }

#' MINQUE0 Component Test (R implementation)
#'
#' R implementation of the MINQUE0-based variance component test for one
#' or more components. It compares a full-model variance component vector
#' \code{vcs} to a null-model vector \code{vcs_h0} for the complement of
#' the tested indices.
#'
#' @param KList A list of \eqn{n \times n} kernel matrices.
#' @param vcs A numeric vector of full-model variance component estimates.
#' @param vcs_h0 A numeric vector of variance component estimates under
#'   the null model, corresponding to the kernels in \code{KList[-index_interest]}.
#' @param index_interest Integer vector (1-based) indicating which components
#'   are being tested jointly.
#' @param wgt A numeric vector of weights of length \code{length(KList)},
#'   default \code{rep(1, length(KList))}.
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item \code{pvalue} P-value computed from the eigenvalues of the product
#'     of the test covariance matrix and the null covariance matrix using
#'     \code{CompQuadForm::davies()}.
#'   \item \code{err} If present with value \code{-1}, indicates a failure
#'     in matrix inversion or invalid entries.
#' }
#'
#' @details
#' This function is used internally by \code{MNQTest0_ChiR()} to compute
#' component-wise or group-wise variance component tests, given the null
#' variance components \code{vcs_h0}.
#'
#' @examples
#' \dontrun{
#' n <- 30
#' set.seed(1)
#' K1 <- diag(n)
#' K2 <- tcrossprod(matrix(rnorm(n * 3), n, 3))
#' KList <- list(K1, K2)
#' vcs    <- c(1, 0.5)
#' vcs_h0 <- c(1)   # null model drops component 2
#' MINQUE0_Component_R(KList, vcs, vcs_h0,
#'                     index_interest = 2,
#'                     wgt = c(1, 1))
#' }
#'
#' @keywords internal
MINQUE0_Component_R<-function(KList, vcs,vcs_h0,index_interest,wgt=rep(1,length(KList)))
{
  # index_interest<-c(2,3)
  VCs<-length(KList)
  C.mat<-matrix(NA,nrow=VCs,ncol=VCs)
  # V.sigma   <- Reduce(f = "+", x = mapply(FUN = "*", KList,est_vcs, SIMPLIFY = FALSE))

  for (i in 1:VCs)
  {
    for (j in i:VCs)
    {
      C.mat[j,i]<-C.mat[i,j]<-sum(t( KList[[i]])* KList[[j]])
    }
  }
  C.mat<-tryCatch(solve(C.mat),
                  error = function(e){
                    print(e)
                    return(as.matrix(-999))
                  }
  )
  if (C.mat[1,1]==-999 || any(is.na(C.mat)))
  {
    return(list(err=-1))
  }
  if (any(diag(C.mat)<0)) {
    return(list(err=-1))
  }
  QList<-list()
  for (i in 1:VCs)
  {
    sigmaH_i   <- Reduce(f = "+", x = mapply(FUN = "*", KList,C.mat[i,], SIMPLIFY = FALSE)) # to calculate Psi=Phi(Sigma)
    QList[[i]]<-sigmaH_i
  }
  ratio<-sum((vcs*wgt)[index_interest])
  sigma_total<-Reduce(f = "+", x = mapply(FUN = "*", QList[index_interest],wgt[index_interest], SIMPLIFY = FALSE))
  V.sigma_hat_h0   <- Reduce(f = "+", x = mapply(FUN = "*", KList[-index_interest],vcs_h0, SIMPLIFY = FALSE))
  eigen_test<-Re(eigen(sigma_total%*%V.sigma_hat_h0, only.values = TRUE)$values)
  pvalue<-davies(q=ratio,lambda = eigen_test)$Qq



  return(list(pvalue=pvalue))
}


#' Iterative MINQUE-Based Normal Approximation Test (R implementation)
#'
#' R implementation of normal-approximation variance component tests based on
#' the iterative MINQUE covariance structure. It constructs Wald-type
#' statistics for each component and a joint chi-square test for a subset
#' of components.
#'
#' @param KList A list of \eqn{n \times n} kernel matrices.
#' @param vcs A numeric vector of variance component estimates.
#' @param test_index Integer vector of indices (1-based) specifying which
#'   components are included in the joint chi-square test (default \code{c(2, 3)}).
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item \code{z} Numeric vector of Wald z-statistics for each component.
#'   \item \code{components} Numeric vector of one-sided p-values
#'     (upper-tail) for each component.
#'   \item \code{chi} The chi-square statistic for the components specified
#'     in \code{test_index}.
#'   \item \code{overall} P-value for the joint chi-square test, with a
#'     50:50 mixture adjustment applied when the sum of standardized
#'     components is positive.
#'   \item \code{df} Degrees of freedom of the joint chi-square test
#'     (length of \code{test_index}).
#'   \item \code{err} If present with value \code{-1}, indicates failure
#'     in matrix inversion or non-numeric standardized estimates.
#' }
#'
#' @details
#' This is the R analogue of \code{IMNQTest_Normal_arm()}, and is used to
#' compute approximate normal and chi-square tests for variance components
#' based on the MINQUE covariance structure. The function internally builds
#' a MINQUE covariance matrix, inverts it, constructs the \code{C.mat}
#' information matrix, and then standardizes selected components.
#'
#' @examples
#' \dontrun{
#' n <- 40
#' set.seed(1)
#' K1 <- diag(n)
#' K2 <- tcrossprod(matrix(rnorm(n * 3), n, 3))
#' KList <- list(K1, K2)
#' vcs <- c(1, 0.5)
#' IMNQTest_Normal_R(KList, vcs, test_index = 1:2)
#' }
#'
#' @keywords internal
IMNQTest_Normal_R<-function(KList, vcs,test_index=c(2,3))
{
  VCs<-length(vcs)
  vcs.standard<-vcs
  nN<-nrow(KList[[1]])
  C.mat<-matrix(NA,nrow=VCs,ncol=VCs)
  V.sigma   <- Reduce(f = "+", x = mapply(FUN = "*", KList, vcs, SIMPLIFY = FALSE)) # to calculate Psi=Phi(Sigma)
  V.sigma<-tryCatch(solve(V.sigma,tol=1e-16),
                    error = function(e){
                      print("cannot inversed")
                      return(as.matrix(-999))
                    }
  )
  if (V.sigma[1,1]==-999 || any(is.na(V.sigma)))
  {
    return(list(err=-1))
  }
  SigmaVi<-list()
  for (i in 1:VCs) {
    SigmaVi[[i]]<-V.sigma%*%KList[[i]]
  }
  for (i in 1:VCs)
  {
    for (j in 1:VCs)
    {
      C.mat[i,j]<-sum(t( SigmaVi[[i]])* SigmaVi[[j]])*0.5
    }
  }

  C.mat<-tryCatch(solve(C.mat,tol=1e-16),
                  error = function(e){
                    print("cannot inversed")
                    return(as.matrix(-999))
                  }
  )
  if (C.mat[1,1]==-999 || any(is.na(C.mat)))
  {
    return(list(err=-1))
  }

  if (length(test_index)==1) {
    vcs.standard[test_index]<-vcs[test_index]/sqrt(C.mat[test_index,test_index])
  }else
  {
    vcs.standard[test_index]<-solve(sqrtm(C.mat[test_index,test_index]),tol=1e-16)%*%vcs[test_index]
  }

  if (!all(is.numeric(vcs.standard))) {
    return(list(err=-1))
  }
  chi<-t(vcs.standard[test_index])%*%vcs.standard[test_index]
  pvalue<-c()
  zstat<-c()
  for (id in 1:VCs)
  {
    wald<-vcs[id]/sqrt(C.mat[id,id])

    pvalue_i<-pnorm(q=wald, mean = 0, sd = 1, lower.tail = F)
    pvalue<-c(pvalue,pvalue_i)
    zstat<-c(zstat,wald)
  }
  # chiP<-pchisq(chi,df = 2,lower.tail = F)
  chiP<-ifelse(sum(vcs.standard[test_index])>0,pchisq(chi,length(test_index),lower.tail = F)/2,1)
  return(list(z=zstat,components=pvalue,chi=chi,overall=chiP,df= length(test_index)))
}
