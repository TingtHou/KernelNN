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

#' Generate Kernel Matrices from an Input Data Matrix or List
#'
#' Constructs a list of kernel (covariance) matrices from an input data
#' matrix or a list of matrices, using a specified kernel type.
#'
#' @param DataMatrix A numeric matrix or a \code{list} of numeric matrices.
#'   Each matrix is typically an \eqn{n \times p} feature/genotype matrix,
#'   where rows correspond to samples.
#' @param kernelName A character string specifying the kernel type to use.
#'   Passed to \code{findKernel()}. Supported options include (but are not
#'   limited to):
#'   \itemize{
#'     \item \code{"CAR"} – conditional autoregressive kernel
#'     \item \code{"identity"} – identity kernel
#'     \item \code{"product"} – product kernel
#'     \item \code{"polynomial"} – polynomial kernel
#'     \item \code{"Gaussian"} – Gaussian (RBF) kernel
#'   }
#'
#' @return
#' A list of kernel matrices, one for each input matrix in \code{DataMatrix},
#' where each element is the result of \code{findKernel(kernelName, geno = G)}.
#'
#' @details
#' If \code{DataMatrix} is not already a list, it is wrapped into a single-element
#' list internally. Each element is coerced to a numeric matrix via \code{as.matrix()}
#' before being passed to \code{findKernel()}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' G  <- matrix(rnorm(50 * 10), 50, 10)
#' Ks <- KernelGenerator(DataMatrix = G, kernelName = "Gaussian")
#' length(Ks)
#' dim(Ks[[1]])
#' }
#'
#' @export
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
