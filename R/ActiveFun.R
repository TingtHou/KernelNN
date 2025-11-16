#' Polynomial Expansion of Kernel Matrices
#'
#' Constructs a list of kernel matrices based on polynomial combinations
#' of the input kernels up to a specified order.
#'
#' The function first augments the input kernel list with:
#' \itemize{
#'   \item a constant kernel \eqn{J}, an \eqn{n \times n} matrix of ones, and
#'   \item later an identity matrix \eqn{I_n} appended at the end.
#' }
#' Then it generates all multinomial combinations of these base kernels of
#' a given polynomial order, and for each combination, forms a kernel by
#' taking the element-wise product of the base kernels raised to the
#' corresponding exponents.
#'
#' @param Kernels A list of \eqn{n \times n} kernel (covariance) matrices.
#'   All matrices must have the same dimension.
#' @param order Integer; the polynomial order used to construct combinations
#'   of the base kernels (default \code{2}).
#' @param prior A numeric vector of prior variance component values
#'   (currently not used inside the function, but kept for compatibility
#'   with other MINQUE-related interfaces).
#'
#' @return
#' A list with one element:
#' \itemize{
#'   \item \code{KList} — a list of kernel matrices including all polynomial
#'   combinations and the identity matrix at the end.
#' }
#'
#' @details
#' This function relies on an internal helper \code{multinomial(m, order)}
#' that returns all exponent patterns for multinomial terms of the specified
#' order, along with their coefficients. It assumes that \code{multinomial()}
#' returns a list with components \code{$coef} and \code{$exp}, where
#' \code{exp} is a matrix whose rows are exponent vectors.
#'
#' @examples
#' \dontrun{
#' n <- 5
#' K1 <- diag(n)
#' K2 <- matrix(0.5, n, n); diag(K2) <- 1
#' poly_kernels <- Polynomial(list(K1, K2), order = 2)
#' length(poly_kernels$KList)
#' }
#'
#' @export
Polynomial<-function(Kernels,order=2,prior=c(1,1,1))
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

#' Linear Kernel Expansion with Identity
#'
#' Takes a list of kernel matrices and appends the identity matrix as an
#' additional kernel component.
#'
#' This is a simple special case of kernel construction where the model
#' includes all original kernels plus an independent residual term
#' represented by the identity matrix.
#'
#' @param Kernels A list of \eqn{n \times n} kernel matrices with matching
#'   dimensions.
#'
#' @return
#' A list named \code{KList} containing the original kernels with an
#' identity matrix appended as the last element.
#'
#' @examples
#' n <- 5
#' K1 <- diag(n)
#' K2 <- matrix(0.5, n, n); diag(K2) <- 1
#' lin_kernels <- Linear(list(K1, K2))
#' length(lin_kernels$KList)  # original 2 + 1 identity
#'
#' @export
Linear<-function(Kernels)
{
  nN<-nrow(Kernels[[1]])
  I <- diag(rep(1, nN));
  Kernels[[length(Kernels)+1]]<-I
  return(KList=Kernels)
}
