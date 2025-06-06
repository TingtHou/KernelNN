
KernelNN<-function(y,Data_Matrix,activeFun="linear", MINQUE="MINQUE0", Kernel="Product", Testing=F, OverallOnly=T,TestingID=NULL,iteration=100, threshold=1e-5)
{
  KernelList<-KernelGenerator(DataMatrix = Data_Matrix,kernelName = Kernel)
  HiddenKernel<-switch(activeFun,
                       linear = {Linear(KernelList)},
                       polynomial = {polynomial(KernelList);},
                       stop("Invalied Kernel Name!")
  )
  if(MINQUE=="MINQUE0")
  {
    vcs.result<-MINQUE0(KList = HiddenKernel, y = y)
    VCs <-vcs.result$vcs
    if (!all(is.numeric(VCs))) {
      return(rep(NA, length(KList)))
    }
    if (Testing)
    {
      if (OverallOnly) {
        TestingID<-NULL

      }else
      {
        if (is.null(TestingID)) {
          TestingID=seq_len(length(HiddenKernel))
        }
      }
      result<-MNQTest0_Chi(y,KList = HiddenKernel,vcs = VCs,TestingID)
    }
  }else if(MINQUE=="IMINQUE")
  {
    vcs.result<-IMINQUE(KList = HiddenKernel, y = y)
    VCs <-vcs.result$vcs
    if (!all(is.numeric(VCs))) {
      return(rep(NA, length(KList)))
    }
    if (Testing)
    {
      if (OverallOnly) {
        TestingID<-NULL

      }else
      {
        if (is.null(TestingID)) {
          TestingID=seq_len(length(HiddenKernel))
        }
      }
      result<-IMNQTest_Normal(y,KList = HiddenKernel,vcs = VCs,TestingID)
    }
  }

}
