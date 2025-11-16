SingleVarianceComp_fullRank<-function(phe,KernelList,TestingIndex=2)
{

  It<-1
  vcs.result<-MINQUE0(KList = KernelList[-TestingIndex],y = phe)
  VCs<-vcs.result$vcs
  return(list(vcs=VCs,Iterate=It))
}

#' MINQUE0-Based Variance Component Test (C++ Implementation)
#'
#' Performs an overall chi-square type test for multiple variance components
#' and optional single-component tests using MINQUE0-based quadratic forms
#' (C++ implementation).
#'
#' @param KList A list of kernel (covariance) matrices corresponding to
#'   variance components in the model.
#' @param vcs A numeric vector of variance component estimates under the
#'   full model, typically obtained from \code{MINQUE0()} or \code{IMINQUE()}.
#' @param ComponentID Optional integer vector giving the indices of
#'   components to test individually. If \code{NULL}, only the overall test
#'   is performed.
#' @param y Optional numeric response vector. Required if
#'   \code{ComponentID} is not \code{NULL}, because single-component tests
#'   construct a null model by re-estimating variance components with the
#'   target component removed.
#' @param OverallID Integer vector of indices specifying which components
#'   are included in the overall test (default: all components except the
#'   last, \code{seq_len(length(KList) - 1)}).
#'
#' @return
#' If successful, a list with components:
#' \itemize{
#'   \item \code{overall} Overall test p-value.
#'   \item \code{components} A numeric vector of p-values for each component
#'     in \code{ComponentID}, in the same order.
#' }
#'
#' If the underlying C++ routine fails (returns \code{"err"}) or an error is
#' thrown, a vector of \code{NA} of length \code{length(KList)} is returned.
#'
#' @details
#' The overall test is computed via \code{MNQTest0_Overall_arm()}, and each
#' single-component test uses \code{MNQTest0_Component_arm()} with a null
#' model estimated by \code{SingleVarianceComp_fullRank()}. All weights
#' are currently set to 1.
#'
#' @examples
#' \dontrun{
#' n <- 50
#' set.seed(1)
#' y  <- rnorm(n)
#' K1 <- tcrossprod(matrix(rnorm(n * 5), n, 5))
#' K2 <- tcrossprod(matrix(rnorm(n * 5), n, 5))
#' vcs <- MINQUE0(list(K1, K2, diag(n)), y)$vcs
#' MNQTest0_Chi(KList = list(K1, K2, diag(n)),
#'              vcs = vcs,
#'              ComponentID = c(1, 2),
#'              y = y)
#' }
#'
#' @export
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

#' MINQUE0-Based Variance Component Test (R Implementation)
#'
#' R implementation of the MINQUE0-based overall and component-wise
#' variance component tests, analogous to \code{MNQTest0_Chi()} but using
#' R versions of the underlying MINQUE0 routines.
#'
#' @inheritParams MNQTest0_Chi
#'
#' @return
#' If successful, a list with components:
#' \itemize{
#'   \item \code{overall} Overall test p-value.
#'   \item \code{components} A numeric vector of p-values for each component
#'     in \code{ComponentID}.
#' }
#'
#' If an error occurs or the underlying R routines return \code{"err"}, a
#' vector of \code{NA} of length \code{length(KList)} is returned.
#'
#' @details
#' This function uses \code{MINQUE0_Overall_R()} and \code{MINQUE0_Component_R()}
#' instead of the C++-based \code{*_arm} functions, but is otherwise
#' parallel in structure to \code{MNQTest0_Chi()}.
#'
#' @examples
#' \dontrun{
#' n <- 50
#' set.seed(1)
#' y  <- rnorm(n)
#' K1 <- tcrossprod(matrix(rnorm(n * 5), n, 5))
#' K2 <- tcrossprod(matrix(rnorm(n * 5), n, 5))
#' vcs <- MINQUE0(list(K1, K2, diag(n)), y)$vcs
#' MNQTest0_ChiR(KList = list(K1, K2, diag(n)),
#'               vcs = vcs,
#'               ComponentID = c(1, 2),
#'               y = y)
#' }
#'
#' @export
MNQTest0_ChiR <- function(KList, vcs ,ComponentID=NULL, y=NULL,OverallID=seq_len(length(KList)-1)) {
  tryCatch({
    if (!is.null(ComponentID) && is.null(y)) {
      stop("Error: When testing a specific variance component, please ensure that both ComponentID and the response vector y are provided.\n")
    }
    # Overall test
    Overall <- MINQUE0_Overall_R(KList = KList, vcs  = vcs, wgt = rep(1,length(KList)),index_interest =OverallID )
    if ("err" %in% names(Overall)) {
      return(rep(NA, length(KList)))
    }

    # Store p-values
    components <-rep(NA,length(ComponentID)) #matrix(nrow=length(TestingID),ncol=2)
    #  row_result[1] <- Overall$pvalue

    # Additional tests
    for (kid in seq_along(ComponentID)) {
      FullRankTesting<-SingleVarianceComp_fullRank(phe = y, KernelList = KList,TestingIndex = ComponentID[[kid]])
      MNQ_Single_plugin_par <- MINQUE0_Component_R(KList = KList, vcs  = vcs,
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

#' Iterative MINQUE-Based Normal Approximation Test
#'
#' Performs variance component tests based on an iterative MINQUE (IMINQUE)
#' framework with normal approximation, using C++ backend
#' \code{IMNQTest_Normal_arm()}.
#'
#' @param KList A list of kernel (covariance) matrices.
#' @param vcs A numeric vector of variance component estimates under the
#'   full model.
#' @param ComponentID Optional integer vector of component indices for which
#'   normal-approximation test statistics/p-values are extracted. If
#'   \code{NULL}, only the overall result is meaningful.
#' @param OverallID Integer vector of indices specifying components to be
#'   included in the overall test (default: all except the last).
#'
#' @return
#' If successful, a list with components:
#' \itemize{
#'   \item \code{overall} The overall test statistic or p-value (depending on
#'     the output of \code{IMNQTest_Normal_arm()}).
#'   \item \code{components} A subset of component-wise results indexed by
#'     \code{ComponentID}.
#' }
#'
#' If an error occurs, a list with \code{overall = NA} and
#' \code{components = NA} (for each requested component) is returned.
#'
#' @details
#' This is a wrapper around the C++ function \code{IMNQTest_Normal_arm()},
#' which returns a list with elements \code{overall} and \code{components}.
#' The \code{components} element is then subset by \code{ComponentID}.
#'
#' @examples
#' \dontrun{
#' n <- 50
#' set.seed(1)
#' K1 <- tcrossprod(matrix(rnorm(n * 5), n, 5))
#' K2 <- tcrossprod(matrix(rnorm(n * 5), n, 5))
#' vcs <- rep(1, 3)  # placeholder, typically from IMINQUE
#' IMNQTest_Normal(KList = list(K1, K2, diag(n)),
#'                 vcs = vcs,
#'                 ComponentID = c(1, 2))
#' }
#'
#' @export
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
