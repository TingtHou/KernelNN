multinomial<-function(ele, deg)
{
  exp_coef<-list()
  for (i in 1:ele) {
    exp_coef[[i]]<-c(0:deg)
  }
  all<-expand.grid(exp_coef)
  all$sum<-rowSums(all)
  exp_c<-all[which(all$sum==deg),-ncol(all)]
  # exp_c<-exp_c%>%arrange(desc(vars(colnames(exp_c))))
  coef_c<-apply(exp_c, 1, multichoose)
  list(exp=exp_c,coef=coef_c)
}

KernelGenerator <- function(DataMatrix, kernelName) {
  grm<-list()
  ### load multiple kernels from grm files

  KernelList<-list()
  mkernel<-c()
  if (!is.list(DataMatrix)) {
    DataMatrix<-list(DataMatrix)
  }
  nVCs<-length(DataMatrix)
  for (vc in seq_len(nVCs)) {
    G<-DataMatrix[[vc]]
    G<-as.matrix(G)
    grm.kernel<-findKernel(kernelName, geno = G);
    KernelList[[vc]]<-grm.kernel
  }

  return(KernelList)
}
