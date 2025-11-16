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



IMNQTest_Iterate_wald<-function(KList, est_vcs)
{
  VCs<-length(est_vcs)
  nN<-nrow(KList[[1]])
  C.mat<-matrix(NA,nrow=VCs,ncol=VCs)
  V.sigma   <- Reduce(f = "+", x = mapply(FUN = "*", KList, est_vcs, SIMPLIFY = FALSE)) # to calculate Psi=Phi(Sigma)
  V.sigma<-solve(V.sigma)
  V.sigma_sqrt<-sqrtm(V.sigma)
  SigmaVi<-list()
  for (i in 1:VCs) {
    SigmaVi[[i]]<-V.sigma%*%KList[[i]]
  }
  for (i in 1:VCs)
  {
    for (j in i:VCs)
    {
      C.mat[j,i]<- C.mat[i,j]<-sum(t( SigmaVi[[i]])* SigmaVi[[j]])
    }
  }
  # if(abs(det(C.mat))<1e-10)
  # {
  #   return(list(err=-1))
  # }

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

  eigen_test<-matrix(nrow = nN,ncol =VCs-2 )
  for (i in 2:seq_len(VCs-1))
  {
    sigmaH   <- Reduce(f = "+", x = mapply(FUN = "*", KList,C.mat[i,], SIMPLIFY = FALSE)) # to calculate Psi=Phi(Sigma)

    KL<-V.sigma_sqrt%*%sigmaH%*%V.sigma_sqrt
    eigen_test[,i]<-eigen(KL,symmetric=T)$values
  }


  pvalue<-c()
  zstat<-c()
  for (id in 2:seq_len(VCs-1))
  {

    pvalue_i<-davies(q=est_vcs[id],lambda = eigen_test[,id])$Qq
    pvalue<-c(pvalue,pvalue_i)
    zstat<-c(zstat,est_vcs[id])
  }
  return(list(z=zstat,pvalue=c(NA,pvalue)))

}
