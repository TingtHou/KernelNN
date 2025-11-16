#' MINQUE0 Variance Component Estimation (C++)
#'
#' Computes MINQUE(0) estimates of variance components using a list of kernel
#' matrices and a response vector \code{y}.
#'
#' @param KList A list of \eqn{n \times n} symmetric kernel (covariance) matrices.
#' @param y A numeric vector of length \eqn{n}.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{vcs} A numeric vector of estimated variance components.
#'   \item \code{it} An integer, currently always 1 (for consistency with IMINQUE).
#' }
#'
#' @examples
#' \dontrun{
#' n <- 30
#' set.seed(1)
#' y <- rnorm(n)
#' K1 <- diag(n)
#' K2 <- tcrossprod(matrix(rnorm(n*5), n, 5))
#' MINQUE0(list(K1, K2), y)
#' }
#'
#' @export
#' @name MINQUE0
NULL


#' MINQUE Variance Component Estimation with Priors (C++)
#'
#' Computes MINQUE estimates of variance components for a linear mixed model,
#' given a list of kernel matrices, a response vector \code{y}, and a vector of
#' prior variance component values. The prior vector defines the working
#' covariance matrix \eqn{V = \sum_i \text{prior}_i K_i}.
#'
#' @param KList A list of \eqn{n \times n} symmetric kernel matrices.
#' @param y A numeric vector of length \eqn{n}.
#' @param prior A numeric vector of prior variance component values
#'   (same length as \code{KList}).
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{vcs} Estimated variance components.
#'   \item \code{it} An integer, currently always 1.
#'   \item \code{err} If present (value \code{-1}), indicates failure due to
#'     non-invertible matrices.
#' }
#'
#' @examples
#' \dontrun{
#' n <- 30
#' set.seed(1)
#' y <- rnorm(n)
#' K1 <- diag(n)
#' K2 <- tcrossprod(matrix(rnorm(n*5), n, 5))
#' MINQUE(list(K1, K2), y, prior = c(1, 0.5))
#' }
#'
#' @export
#' @name MINQUE
NULL


#' Iterative MINQUE (IMINQUE) Variance Component Estimation (C++)
#'
#' Performs iterative MINQUE estimation. At each iteration, the variance
#' component estimates from \code{MINQUE()} are used as the new priors until
#' convergence or until reaching \code{epoch} iterations.
#'
#' @param KList A list of \eqn{n \times n} symmetric kernel matrices.
#' @param y A numeric response vector.
#' @param prior Initial prior values (numeric vector). If empty, it defaults
#'   to a vector of ones of length \code{length(KList)}.
#' @param epoch Maximum number of iterations (default: \code{100}).
#' @param threshold Convergence tolerance for the Euclidean norm of successive
#'   estimates (default: \code{1e-4}).
#' @param echo Logical; if \code{TRUE}, prints iteration details.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{vcs} Final variance component estimates.
#'   \item \code{it} Number of iterations performed.
#'   \item \code{err} If present (value \code{-1}), indicates failure in the
#'     underlying \code{MINQUE()} call.
#' }
#'
#' @examples
#' \dontrun{
#' n <- 40
#' set.seed(1)
#' y <- rnorm(n)
#' K1 <- diag(n)
#' K2 <- tcrossprod(matrix(rnorm(n*5), n, 5))
#' IMINQUE(list(K1, K2), y, prior = c(1,1), echo = FALSE)
#' }
#'
#' @export
#' @name IMINQUE
NULL

#' Overall MINQUE0 Test (arm-level C++ backend)
#'
#' Low-level C++ backend for the MINQUE0-based overall variance component test.
#' This function constructs quadratic forms from the kernel matrices, inverts
#' the covariance structure, and computes a test statistic and p-value using
#' the Davies method implemented in \pkg{CompQuadForm}.
#'
#' @param KList A list of \eqn{n \times n} kernel (covariance) matrices.
#' @param vcs A numeric vector of variance component estimates under the
#'   full model (one entry per kernel in \code{KList}).
#' @param index_interest A numeric or integer vector of indices (1-based) of
#'   components to be included in the overall test.
#' @param wgt A numeric vector of weights of the same length as \code{vcs}.
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item \code{z} Numeric value of the test statistic (ratio).
#'   \item \code{pvalue} P-value computed via \code{CompQuadForm::davies()}.
#'   \item \code{err} If present with value \code{-1}, indicates that the
#'     internal matrix inversion failed.
#' }
#'
#' @details
#' This is an internal backend called by higher-level R functions such as
#' \code{MNQTest0_Chi()}, and is not intended to be used directly by most users.
#'
#' @examples
#' \dontrun{
#' n <- 20
#' set.seed(1)
#' K1 <- diag(n)
#' K2 <- tcrossprod(matrix(rnorm(n * 3), n, 3))
#' vcs <- c(1, 0.5)
#' MNQTest0_Overall_arm(KList = list(K1, K2),
#'                      vcs = vcs,
#'                      index_interest = c(1, 2),
#'                      wgt = rep(1, 2))
#' }
#'
#' @keywords internal
#' @name MNQTest0_Overall_arm
NULL


#' Single/Multiple Component MINQUE0 Test (arm-level C++ backend)
#'
#' Low-level C++ backend for single or multiple variance component tests
#' under MINQUE0. It compares a full model (with variance components
#' \code{vcs}) to a null model with variance components \code{vcs_h0},
#' using kernel matrices in \code{KList}.
#'
#' @param KList A list of \eqn{n \times n} kernel matrices.
#' @param vcs A numeric vector of full-model variance component estimates.
#' @param vcs_h0 A numeric vector of null-model variance component estimates
#'   (typically produced by \code{SingleVarianceComp_fullRank()}).
#' @param index_interest An integer vector (1-based) indicating which
#'   components are being jointly tested.
#' @param wgt A numeric weight vector (same length as \code{vcs}).
#'
#' @return
#' A list with component:
#' \itemize{
#'   \item \code{pvalue} P-value for the component test.
#'   \item \code{err} If present with value \code{-1}, indicates a failure
#'     in the internal matrix inversions or numerical issues.
#' }
#'
#' @details
#' This function is used internally by \code{MNQTest0_Chi()} for each
#' tested component or group of components.
#'
#' @examples
#' \dontrun{
#' n <- 20
#' set.seed(1)
#' K1 <- diag(n)
#' K2 <- tcrossprod(matrix(rnorm(n * 3), n, 3))
#' KList <- list(K1, K2)
#' vcs    <- c(1, 0.5)
#' vcs_h0 <- c(1)  # e.g., null model excludes the second component
#' MNQTest0_Component_arm(KList, vcs, vcs_h0,
#'                        index_interest = 2L,
#'                        wgt = rep(1, 2))
#' }
#'
#' @keywords internal
#' @name MNQTest0_Component_arm
NULL


#' Iterative MINQUE-Based Normal Approximation Test (arm-level C++ backend)
#'
#' Low-level C++ routine implementing normal-approximation tests for variance
#' components based on an iterative MINQUE covariance structure. It computes
#' Wald-type z-statistics for each component as well as a chi-square type
#' overall test for a subset of components.
#'
#' @param KList A list of \eqn{n \times n} kernel matrices.
#' @param vcs A numeric vector of variance component estimates.
#' @param index_interest Integer vector (1-based) indicating which components
#'   are included in the joint overall test.
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item \code{z} A numeric vector of Wald z-statistics, one per component.
#'   \item \code{components} A numeric vector of (one-sided) p-values for each
#'     component, obtained from the standard normal distribution.
#'   \item \code{overall} P-value for the joint chi-square test of the
#'     components specified by \code{index_interest}.
#'   \item \code{df} Degrees of freedom for the overall chi-square test
#'     (length of \code{index_interest}).
#'   \item \code{err} If present with value \code{-1}, indicates failure
#'     in matrix inversion or decomposition.
#' }
#'
#' @details
#' This is the backend for \code{IMNQTest_Normal()}, which provides a more
#' user-friendly R interface. Most users should call \code{IMNQTest_Normal()}
#' instead of this arm-level function.
#'
#' @examples
#' \dontrun{
#' n <- 20
#' set.seed(1)
#' K1 <- diag(n)
#' K2 <- tcrossprod(matrix(rnorm(n * 3), n, 3))
#' KList <- list(K1, K2)
#' vcs <- c(1, 0.5)
#' IMNQTest_Normal_arm(KList, vcs, index_interest = 1:2)
#' }
#'
#' @keywords internal
#' @name IMNQTest_Normal_arm
NULL
