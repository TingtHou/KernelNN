#' MINQUE0 Variance Component Estimator
#'
#' Computes MINQUE(0) estimates of variance components for a linear mixed model
#' given a list of kernel (covariance) matrices and a response vector.
#'
#' The method solves a quadratic estimating equation based on the trace of
#' products of kernel matrices and quadratic forms in \code{y}.
#'
#' @param KList A list of length \eqn{q} of \eqn{n \times n} symmetric
#'   (semi-)positive definite kernel matrices, one for each variance component.
#' @param y A numeric vector of length \eqn{n} containing the response.
#'
#' @returns
#' A list with component:
#' \item{vcs}{A numeric vector of length \code{length(KList)} containing
#'   the estimated variance components (in the same order as \code{KList}).}
#'
#' @examples
#' n <- 50
#' set.seed(1)
#' y <- rnorm(n)
#' K1 <- diag(n)
#' K2 <- tcrossprod(matrix(rnorm(n * 5), n, 5))
#' fit0 <- MINQUE0_R(list(K1, K2), y)
#' fit0$vcs
#'
#' @export
MINQUE0_R <- function(KList, y)
{
  VCs <- length(KList)
  C.mat <- matrix(NA, nrow = VCs, ncol = VCs)
  RightY <- c()
  for (i in 1:VCs)
  {
    for (j in i:VCs)
    {
      C.mat[j, i] <- C.mat[i, j] <- sum(t(KList[[i]]) * KList[[j]])
    }
    right <- t(y) %*% KList[[i]] %*% y
    RightY <- c(RightY, right)
  }
  C.mat <- solve(C.mat)
  est.vcs <- C.mat %*% RightY
  return(list(vcs = est.vcs,it=1))
}

#' MINQUE Variance Component Estimator with Prior Weights
#'
#' Computes MINQUE estimates of variance components for a linear mixed model
#' given a list of kernel matrices, a response vector, and prior variance
#' component values. The prior values are used to construct the working
#' covariance matrix \eqn{V = \sum_k \text{prior}_k K_k}.
#'
#' @param KList A list of length \eqn{q} of \eqn{n \times n} symmetric
#'   kernel (covariance) matrices.
#' @param y A numeric vector of length \eqn{n} containing the response.
#' @param prior A numeric vector of length \code{length(KList)} giving prior
#'   values for the variance components (defaults to all ones).
#'
#' @returns
#' A list with component:
#' \item{vcs}{A numeric vector of length \code{length(KList)} containing
#'   the estimated variance components.}
#'
#' If the method fails due to a non-invertible matrix, the function returns
#' \code{list(err = -1)}.
#'
#' @examples
#' n <- 50
#' set.seed(1)
#' y <- rnorm(n)
#' K1 <- diag(n)
#' K2 <- tcrossprod(matrix(rnorm(n * 5), n, 5))
#' prior <- c(1, 0.5)
#' fit <- MINQUE_R(list(K1, K2), y, prior = prior)
#' if (!"err" %in% names(fit)) fit$vcs
#'
#' @export
MINQUE_R <- function(KList, y, prior = rep(1, length(KList)))
{
  VCs <- length(KList)
  C.mat <- matrix(NA, nrow = VCs, ncol = VCs)
  V <- Reduce(f = "+", x = mapply(FUN = "*", KList, prior, SIMPLIFY = FALSE))

  Q <- tryCatch(
    solve(V, tol = 1e-6),
    error = function(e) {
      print("cannot inversed")
      return(as.matrix(-999))
    }
  )
  if (Q[1, 1] == -999 || any(is.na(Q)))
  {
    return(list(err = -1))
  }
  KQ <- list()
  Qy <- Q %*% y
  RightY <- c()
  for (i in 1:VCs)
  {
    KQ[[i]] <- KList[[i]] %*% Q
  }
  for (i in 1:VCs)
  {
    for (j in 1:VCs)
    {
      C.mat[i, j] <- sum(as.double(t(KQ[[i]]) * KQ[[j]]))
    }
    right <- t(Qy) %*% KList[[i]] %*% Qy
    RightY <- c(RightY, right)
  }

  C.mat.inv <- tryCatch(
    solve(C.mat, tol = 1e-6),
    error = function(e) {
      print("cannot inversed")
      return(as.matrix(-999))
    }
  )
  if (C.mat.inv[1, 1] == -999 || any(is.na(C.mat.inv)))
  {
    return(list(err = -1))
  }
  est.vcs <- C.mat.inv %*% RightY
  return(list(vcs = est.vcs))
}

#' Iterative MINQUE (IMINQUE) Estimator
#'
#' Computes iterative MINQUE (IMINQUE) estimates of variance components by
#' repeatedly applying \code{MINQUE_R()} until convergence or a maximum
#' number of iterations is reached.
#'
#' At each iteration, the current variance component estimates are used as
#' the new prior vector in \code{MINQUE_R()}, and the algorithm stops when
#' the Euclidean norm between successive estimates is below \code{threshold}
#' or when \code{epoch} iterations have been reached.
#'
#' @param KList A list of \eqn{n \times n} symmetric kernel matrices.
#' @param y A numeric response vector of length \eqn{n}.
#' @param prior A numeric vector of initial prior values for the variance
#'   components. Defaults to a vector of ones of length \code{length(KList)}.
#' @param epoch Integer; maximum number of iterations (default \code{100}).
#' @param threshold Numeric; convergence threshold for the Euclidean norm
#'   between successive estimates (default \code{1e-4}).
#' @param echo Logical; if \code{TRUE}, prints iteration details including
#'   the current estimates and the difference (default \code{FALSE}).
#'
#' @returns
#' A list with components:
#' \item{vcs}{A numeric vector of the final variance component estimates.}
#' \item{it}{The number of iterations performed.}
#'
#' If the underlying call to \code{MINQUE_R()} fails, the function returns
#' \code{list(err = -1)}.
#'
#' @examples
#' n <- 50
#' set.seed(1)
#' y <- rnorm(n)
#' K1 <- diag(n)
#' K2 <- tcrossprod(matrix(rnorm(n * 5), n, 5))
#' fit_im <- IMINQUE_R(list(K1, K2), y, echo = FALSE)
#' if (!"err" %in% names(fit_im)) {
#'   fit_im$vcs
#'   fit_im$it
#' }
#'
#' @export
IMINQUE_R <- function(KList, y, prior = rep(1, length(KList)),
                      epoch = 100, threshold = 1e-4, echo = FALSE)
{
  prior_est <- prior
  it <- 1
  while (it < epoch)
  {
    prior_est <- abs(prior_est)

    est <- MINQUE(KList, y, prior = prior_est)
    if ("err" %in% names(est))
    {
      return(list(err = -1))
    }
    est.vcs <- est$vcs
    diff <- norm((est.vcs - prior_est), type = "2")
    if (echo) {
      cat("IT: ", it, ", est: ", paste(est.vcs, collapse = ","), ", diff: ", diff, "\n")
    }

    if (diff < threshold)
    {
      break
    } else {
      prior_est <- est.vcs
      it <- it + 1
    }
  }
  return(list(vcs = prior_est, it = it))
}
