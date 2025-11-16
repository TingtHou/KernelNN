SingleVarianceComp_fullRank<-function(phe,KernelList,TestingIndex=2)
{

  It<-1
  vcs.result<-MINQUE0(KList = KernelList[-TestingIndex],y = phe)
  VCs<-vcs.result$vcs
  return(list(vcs=VCs,Iterate=It))
}


MNQTest0_Chi <- function(KList, vcs ,ComponentID=NULL, y=NULL,OverallID=seq_len(length(KList)-1)) {
  tryCatch({
    if (!is.null(ComponentID) && is.null(y)) {
      stop("Error: When testing a specific variance component, please ensure that both ComponentID and the response vector y are provided.\n")
    }
    # Overall test
    Overall <- MNQTest0_Overall_arm(KList = KList, vcs  = vcs, wgt = rep(1,length(KList)),index_interest =OverallID )
    if ("err" %in% names(Overall)) {
      return(rep(NA, length(KList)))
    }

    # Store p-values
    components <-rep(NA,length(ComponentID)) #matrix(nrow=length(TestingID),ncol=2)
  #  row_result[1] <- Overall$pvalue

    # Additional tests
    for (kid in seq_along(ComponentID)) {
      FullRankTesting<-SingleVarianceComp_fullRank(phe = y, KernelList = KList,TestingIndex = ComponentID[[kid]])
      MNQ_Single_plugin_par <- MNQTest0_Component_arm(KList = KList, vcs  = vcs,
                                                                 vcs_h0 = FullRankTesting$vcs,
                                                                  index_interest  =  ComponentID[[kid]],
                                                                  wgt = rep(1,length(KList)))
      components[kid] <- if ("err" %in% names(MNQ_Single_plugin_par)) NA else MNQ_Single_plugin_par$pvalue
    }
    result<-list(overall=Overall$pvalue,components=components)
    return(result)

  }, error = function(e) {
    cat(sprintf(conditionMessage(e)))
    result<-list(overall=NA,components=rep(NA,length(ComponentID)) )
    return(rep(NA, length(KList)))  # Return NA for that row
  })
}


IMNQTest_Normal <- function(KList, vcs, ComponentID=NULL,OverallID=seq_len(length(KList)-1))
{
  tryCatch({
  ztest<-IMNQTest_Normal_arm(KList,vcs,OverallID)
  Overall<-ztest$overall
  components <-ztest$components[ComponentID]
  result<-list(overall=Overall,components=components)
  return(result)
  }, error = function(e) {
    cat(sprintf(conditionMessage(e)))
    result<-list(overall=NA,components=rep(NA,length(ComponentID)) )
  })

}
