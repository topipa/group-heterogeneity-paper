#' Kernel matrices and their derivatives

#' @rdname kernels
#' @export
#' @description A function for computing the squared exponential kernel matrix
#' and its first and second derivatives with respect to the input data.
#' @param X First matrix of observation inputs (N1 x D).
#' @param Z Second matrix of observation inputs. By default it is NULL
#' and treated as equal to X (N2 x D).
#' @param alpha Marginal standard deviation hyperparameter. The kernel
#' is multiplied by the square of this.
#' @param lengthscale Lengthscale hyperparameter. The quadratic terms of
#' the kernel are divided by the square of this.
#' @return Returns a list with the squared exponential kernel matrix and
#' its first and second derivatives with respect to X (N1 x N2).
squared_exponential_cov <- function(X, Z=NULL, alpha=1, lengthscale=1) {
  if (is.null(Z))
    Z <- X
  kk <- new(sqr_exponential_kernel, c(alpha, lengthscale), X, Z)
  dims <- c(nrow(Z), nrow(X), ncol(X), ncol(X))
  dK <- array(kk$dK_dX, dim=dims[-4])
  dK <- aperm(dK,perm = c(2,1,3))
  d2K <- array(kk$d2K_dXdX, dim=dims)
  d2K <- aperm(d2K,perm = c(2,1,3,4))
  return (list(K=kk$Kernel, dK=dK, d2K=d2K))
}

#' @rdname kernels
#' @export
#' @description A function for computing the ARD squared exponential kernel matrix
#' and its first and second derivatives with respect to the input data.
#' @param X First matrix of observation inputs (N1 x D).
#' @param Z Second matrix of observation inputs. By default it is NULL
#' and treated as equal to X (N2 x D).
#' @param alpha Marginal standard deviation hyperparameter. The kernel
#' is multiplied by the square of this.
#' @param lengthscale A vector (D x 1) of lengthscale hyperparameters.
#' The quadratic terms of the kernel are divided by the squares of these.
#' If a single number is given, all lengthscales are treated equal
#' (consider using squared_exponential_cov()).
#' @return Returns a list with the squared exponential kernel matrix and
#' its first and second derivatives with respect to X (N1 x N2).
squared_exponential_ard_cov <- function(X, Z=NULL, alpha=1, lengthscale=1) {
  if (is.null(Z))
    Z <- X
  d <- ncol(X)
  
  if (length(lengthscale) == 1)
    lengthscale <- rep(lengthscale,d)
  if (length(lengthscale) != d)
    stop('Lengthscale must be a single number or a vector with length equal to columns of X')
  
  kk <- new(sqr_exponential_ard_kernel, c(alpha),c(lengthscale), X, Z)
  dims <- c(nrow(Z), nrow(X), ncol(X), ncol(X))
  dK <- array(kk$dK_dX, dim=dims[-4])
  dK <- aperm(dK,perm = c(2,1,3))
  d2K <- array(kk$d2K_dXdX, dim=dims)
  d2K <- aperm(d2K,perm = c(2,1,3,4))
  return (list(K=kk$Kernel, dK=dK, d2K=d2K))
}


#' @rdname kernels
#' @export
#' @description A function for computing the linear kernel matrix
#' and its first and second derivatives with respect to the input data.
#' @param X First matrix of observation inputs (N1 x D).
#' @param Z Second matrix of observation inputs. By default it is NULL
#' and treated as equal to X (N2 x D).
#' @param bias bias hyperparameter.
#' @param deviation A vector (D x 1) of hyperparameters describing the deviation of each
#' input variable, sort of equivalent to regression coefficient.
#' @return Returns a list with the linear kernel matrix and
#' its first and second derivatives with respect to X (N1 x N2).
linear_cov <- function(X, Z=NULL, bias=1, deviation=1) {
  if (is.null(Z))
    Z <- X
  d <- ncol(X)
  
  if (length(deviation) == 1)
    deviation <- rep(deviation,d)
  if (length(deviation) != d)
    stop('Deviation must be a single number or a vector with length equal to columns of X')
  
  kk <- new(linear_kernel, c(bias),c(deviation), X, Z)

  dims <- c(nrow(Z), nrow(X), ncol(X), ncol(X))
  dK <- array(kk$dK_dX, dim=dims[-4])
  dK <- aperm(dK,perm = c(2,1,3))
  d2K <- array(kk$d2K_dXdX, dim=dims)
  d2K <- aperm(d2K,perm = c(2,1,3,4))
  return (list(K=kk$Kernel, dK=dK, d2K=d2K))
}
