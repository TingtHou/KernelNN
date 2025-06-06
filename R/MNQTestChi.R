Rcpp::sourceCpp("src/MINQUE_Arma.cpp")
Rcpp::sourceCpp("src/MNQTestChi_Arma.cpp")
SingleVarianceComp_fullRank<-function(phe,KernelList,TestingIndex=2)
{

  It<-1
  vcs.result<-MINQUE0(KList = KernelList[-TestingIndex],y = phe)
  VCs<-vcs.result$vcs
  return(list(vcs=VCs,Iterate=It))
}

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
     vcs.result<-MINQUE0_arma(KList = HiddenKernel, y = y)
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
       MNQTest0_Chi(y,KList = HiddenKernel,vcs = VCs,TestingID)
     }
   }

}


MNQTest0_Chi <- function(y,KList, vcs ,TestingID) {
  tryCatch({

    # Overall test
    Overall <- MNQTest0_Overall(KList = KList, vcs  = vcs,
                                             index_interest = c(1:(length(KList) - 1)),
                                             wgt = rep(1,length(KList)))
    if ("err" %in% names(Overall)) {
      return(rep(NA, length(KList)))
    }

    # Store p-values
    components <- matrix(nrow=length(TestingID),ncol=2)
  #  row_result[1] <- Overall$pvalue

    # Additional tests
    for (kid in seq_along(TestingID)) {
      FullRankTesting<-SingleVarianceComp_fullRank(phe = y, KernelList = KList,TestingIndex = TestingID[[kid]])
      MNQ_Single_plugin_par <- MNQTest0_Component(KList = KList, vcs  = vcs,
                                                                 vcs_h0 = FullRankTesting$vcs,
                                                                  index_interest  =  TestingID[[kid]],
                                                                  wgt = rep(1,length(KList)))
      components[kid,] <-c(TestingID[kid] ,
                            if ("err" %in% names(MNQ_Single_plugin_par)) NA else MNQ_Single_plugin_par$pvalue)
    }
    result<-list(overall=Overall$pvalue,components=components)
    return(result)

  }, error = function(e) {
    cat(sprintf(conditionMessage(e)))
    return(rep(NA, length(KList)))  # Return NA for that row
  })
}
