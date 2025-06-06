
polynomial<-function(Kernels,order=2,prior=c(1,1,1))
{

  nN<-nrow(Kernels[[1]])
  J<-matrix(1,nrow = nN,ncol = nN)
  I <- diag(rep(1, nN));
  K_id<-1
  Kernel_O<-list()
  KernelList<-list()
  Kernel_O[[K_id]]<-J
  K_id<-K_id+1

  for (i in 1:length(Kernels)) {
    Kernel_O[[K_id]]<-Kernels[[i]]
    K_id<-K_id+1
  }

  ### generate different combinations of kernels and coefs
  expand.k<-multinomial(length(Kernel_O),order)
  mul_ceof<-expand.k$coef
  exp_coef<-expand.k$exp

  ### generate new kernels, centered except J matrix
  mkernel<-c()
  K_id<-1
  for (k in 1:nrow(exp_coef))
  {
    mut<-matrix(1,nrow = nN,ncol = nN)
    for (vc in seq_len(length(Kernel_O)))
    {
      mut<-mut*(Kernel_O[[vc]]^(exp_coef[k,vc]))
    }

    KernelList[[K_id]]<-mut
    K_id<-K_id+1
  }


  KernelList[[K_id]]<-I

  return(list(KList=KernelList))
}

Linear<-function(Kernels)
{
  nN<-nrow(Kernels[[1]])
  I <- diag(rep(1, nN));
  Kernels[[length(Kernels)+1]]<-I
  return(list(KList=Kernels))
}
