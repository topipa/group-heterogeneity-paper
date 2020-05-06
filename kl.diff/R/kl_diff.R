#' KL-diff and KL-diff2

#' @rdname KLdiff
#' @export
#' @description A function for computing the KL-diff values for all observations
#' and input variables for a GP with Gaussian observation model.
#' @param y Vector (N1) of observation outputs.
#' @param X Matrix of observation inputs (N1 x D).
#' @param alpha Marginal standard deviation hyperparameter. The kernel
#' is multiplied by the square of this.
#' @param lengthscale A vector (D x 1) of lengthscale hyperparameters.
#' The quadratic terms of the kernel are divided by the squares of these.
#' If a single number is given, all lengthscales are treated equal.
#' @param sigma Noise variance hyperparameter.
#' (NOTE: variance, not standard deviation)
#' @param pointwise Logical indicating whether to return KL-diff values
#' for each observation, or the mean over all observations.
#' @return Returns a matrix (N1 x D) or vector (D) of KL-diff values.
KL_diff_gaussian <- function(y, X, alpha = 1, lengthscale = 1,
                             sigma = 1, pointwise = FALSE) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  k_X_X <- squared_exponential_ard_cov(X,X,alpha,lengthscale)
  V_f <- k_X_X$K - t(k_X_X$K) %*% solve(k_X_X$K + diag(n)*sigma, k_X_X$K)
  V_f = diag(V_f)
  V_y <- V_f + sigma
  
  WB_vec <- solve(k_X_X$K + diag(n)*sigma,y)
  WB_mat <- solve(k_X_X$K + diag(n)*sigma)
  dE_f <- sapply(1:d,dE,k_X_X$dK,WB_vec)
  dV_f <- sapply(1:d,dV,k_X_X$K,k_X_X$dK,WB_mat)
  KLdiff <- sapply(1:d,KL_diff_gaussian_i,dE_f,dV_f,V_y)
  
  
  if (pointwise)
    KLdiff
  else
    colMeans(KLdiff)
}

#' @rdname KLdiff
#' @export
#' @description A function for computing the KL-diff2 values for all observations
#' and input variables for a GP with Gaussian observation model.
#' @param y Vector (N1) of observation outputs.
#' @param X Matrix of observation inputs (N1 x D).
#' @param subset_variables Optional vector that represents a subset
#' of the D variables for which the KL-diff2 values are computed
#' @param alpha Marginal standard deviation hyperparameter. The kernel
#' is multiplied by the square of this.
#' @param lengthscale A vector (D x 1) of lengthscale hyperparameters.
#' The quadratic terms of the kernel are divided by the squares of these.
#' If a single number is given, all lengthscales are treated equal.
#' @param sigma Noise variance hyperparameter.
#' (NOTE: variance, not standard deviation)
#' @param pointwise Logical indicating whether to return KL-diff values
#' for each observation, or the mean over all observations.
#' @return Returns an array (N1 x D x D) or matrix (D x D) of KL-diff2 values.
KL_diff2_gaussian <- function(y, X, subset_variables = NULL, alpha = 1, lengthscale = 1,
                              sigma = 1, pointwise = FALSE) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  if (is.null(subset_variables)) {
    subset_variables <- seq(d)
  }
  
  k_X_X <- squared_exponential_ard_cov(X,X,alpha,lengthscale)
  V_f <- k_X_X$K - t(k_X_X$K) %*% solve(k_X_X$K + diag(n)*sigma, k_X_X$K)
  V_f = diag(V_f)
  V_y <- V_f + sigma
  
  WB_vec <- solve(k_X_X$K + diag(n)*sigma,y)
  WB_mat <- solve(k_X_X$K + diag(n)*sigma)
  
  d2E_f <- array(0, c(n,d,d))
  d2V_f <- array(0, c(n,d,d))
  KLdiff2 <- array(0, c(n,d,d))
  for (j in intersect(subset_variables, (2:d))) {
    for (i in 1:(j-1)) {
      d2E_f[,i,j] <- d2E(i,j,k_X_X$d2K,WB_vec)
      d2V_f[,i,j] <- d2V(i,j,k_X_X$K,k_X_X$dK,k_X_X$d2K,WB_mat)
      KLdiff2[,i,j] <- KL_diff2_gaussian_ij(i,j,d2E_f,d2V_f,V_y)
    }
  }
  
  if (pointwise)
    KLdiff2
  else
    means.along(KLdiff2,1)
}

#' @rdname KLdiff
#' @export
#' @description A function for computing the KL-diff2 values for all observations
#' and input variables for a GP with Gaussian observation model when the
#' model of interest is a linear model. This version uses a linear kernel
#' in the GP in addition to the squared exponential kernel.
#' @param y Vector (N1) of observation outputs.
#' @param X Matrix of observation inputs (N1 x D). Includes the dummy variables.
#' @param Xlin Matrix of observation inputs (N1 x D_numeric) excluding
#' the dummy variables. The linear kernel operates only on these.
#' @param subset_variables Optional vector that represents a subset
#' of the D variables for which the KL-diff2 values are computed
#' @param alpha Marginal standard deviation hyperparameter. The kernel
#' is multiplied by the square of this.
#' @param lengthscale A vector (D x 1) of lengthscale hyperparameters.
#' The quadratic terms of the kernel are divided by the squares of these.
#' If a single number is given, all lengthscales are treated equal.
#' @param sigma Noise variance hyperparameter.
#' (NOTE: variance, not standard deviation)
#' @param bias bias hyperparameter of the linear kernel.
#' @param deviation A vector (D x 1) of hyperparameters describing the deviation
#' of the linear kernel of each
#' input variable, sort of equivalent to regression coefficient.
#' @param pointwise Logical indicating whether to return KL-diff values
#' for each observation, or the mean over all observations.
#' @return Returns an array (N1 x D x D) or matrix (D x D) of KL-diff2 values.
KL_diff2_gaussian_linear <- function(y, X, Xlin, subset_variables = NULL, alpha = 1, lengthscale = 1,
                                     sigma = 1, bias = 1, deviation = 1, pointwise = FALSE) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  if (is.null(subset_variables)) {
    subset_variables <- seq(d)
  }
  dlin <- dim(Xlin)[2]
  
  k_X_X <- squared_exponential_ard_cov(X,X,alpha,lengthscale)
  k_X_X2 <- linear_cov(Xlin,Xlin,bias, deviation)
  
  
  k_X_X$K <- k_X_X$K + k_X_X2$K
  k_X_X$dK[,,1:dlin] <- k_X_X$dK[,,1:dlin] + k_X_X2$dK
  k_X_X$d2K[,,1:dlin,1:dlin] <- k_X_X$d2K[,,1:dlin,1:dlin] + k_X_X2$d2K
  
  
  V_f <- k_X_X$K - t(k_X_X$K) %*% solve(k_X_X$K + diag(n)*sigma, k_X_X$K)
  V_f = diag(V_f)
  V_y <- V_f + sigma
  
  WB_vec <- solve(k_X_X$K + diag(n)*sigma,y)
  WB_mat <- solve(k_X_X$K + diag(n)*sigma)
  
  d2E_f <- array(0, c(n,d,d))
  d2V_f <- array(0, c(n,d,d))
  KLdiff2 <- array(0, c(n,d,d))
  for (j in intersect(subset_variables, (2:d))) {
    for (i in 1:(j-1)) {
      d2E_f[,i,j] <- d2E(i,j,k_X_X$d2K,WB_vec)
      d2V_f[,i,j] <- d2V(i,j,k_X_X$K,k_X_X$dK,k_X_X$d2K,WB_mat)
      KLdiff2[,i,j] <- KL_diff2_gaussian_ij(i,j,d2E_f,d2V_f,V_y)
    }
  }
  
  if (pointwise)
    KLdiff2
  else
    means.along(KLdiff2,1)
}

KL_diff2_bernoulli_laplace_probit <- function(y, X, f, subset_variables = NULL, alpha = 1, lengthscale = 1,
                              pointwise = FALSE) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  if (is.null(subset_variables)) {
    subset_variables <- seq(d)
  }

  jitter <- diag(rep(1e-5,n))
  k_X_X <- squared_exponential_ard_cov(X,X,alpha,lengthscale)
  k_X_X$K <- k_X_X$K + jitter
  W <- (dnorm(f)/pnorm(y*f))^2 + y * f * dnorm(f)/pnorm(y*f)
  Winv <- diag(1/W)

  WB_vec <- solve(k_X_X$K,f)
  WB_mat <- solve(k_X_X$K + Winv)

  E_f <- t(k_X_X$K) %*% WB_vec
  V_f <- k_X_X$K - t(k_X_X$K) %*% solve(k_X_X$K + Winv, k_X_X$K)
  V_f = diag(V_f)

  dE_f <- array(0, c(n,d))
  dV_f <- array(0, c(n,d))
  for (i in 1:d) {
    dE_f[,i] <- dE(i,k_X_X$dK,WB_vec)
    dV_f[,i] <- dV(i,k_X_X$K,k_X_X$dK,WB_mat)
  }

  d2E_f <- array(0, c(n,d,d))
  d2V_f <- array(0, c(n,d,d))
  KLdiff2 <- array(0, c(n,d,d))
  for (j in intersect(subset_variables, (2:d))) {
    for (i in 1:(j - 1)) {
      d2E_f[,i,j] <- d2E(i,j,k_X_X$d2K,WB_vec)
      d2V_f[,i,j] <- d2V(i,j,k_X_X$K,k_X_X$dK,k_X_X$d2K,WB_mat)
      KLdiff2[,i,j] <- KL_diff2_bernoulli_ij(i,j,E_f,V_f,dE_f,dV_f,d2E_f,d2V_f)
    }
  }

  if (pointwise)
    KLdiff2
  else
    means.along(KLdiff2,1)
}

KL_diff2_bernoulli_laplace_probit_pseudo <- function(y, X, subset_variables = NULL, alpha = 1, lengthscale = 1,
                              pointwise = FALSE) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  if (is.null(subset_variables)) {
    subset_variables <- seq(d)
  }

  sigma <- 1e-3

  k_X_X <- squared_exponential_ard_cov(X,X,alpha,lengthscale)
  V_f <- k_X_X$K - t(k_X_X$K) %*% solve(k_X_X$K + diag(n)*sigma, k_X_X$K)
  V_f = diag(V_f)
  V_y <- V_f

  WB_vec <- solve(k_X_X$K + diag(n)*sigma,y)
  WB_mat <- solve(k_X_X$K + diag(n)*sigma)

  d2E_f <- array(0, c(n,d,d))
  d2V_f <- array(0, c(n,d,d))
  KLdiff2 <- array(0, c(n,d,d))
  for (j in intersect(subset_variables, (2:d))) {
    for (i in 1:(j-1)) {
      d2E_f[,i,j] <- d2E(i,j,k_X_X$d2K,WB_vec)
      d2V_f[,i,j] <- d2V(i,j,k_X_X$K,k_X_X$dK,k_X_X$d2K,WB_mat)
      KLdiff2[,i,j] <- KL_diff2_gaussian_ij(i,j,d2E_f,d2V_f,V_y)
    }
  }

  if (pointwise)
    KLdiff2
  else
    means.along(KLdiff2,1)
}

#' @description A function for computing the derivative of the mean prediction
#' with respect to a single input variable for a GP with Gaussian observation
#' model.
#' @param d which input variable.
#' @param dK Array (N2 x N1 x D) with the first derivatives of a Gaussian process
#' kernel matrix with respect to all D inputs.
#' @param WB_vec (K_X_X + sigma^2 I)^(-1) * y
#' @return Returns a vector with the derivatives of the mean prediction.
dE <- function(d,dK,WB_vec) {
  t(dK[,,d]) %*% WB_vec
}

#' @description A function for computing the derivative of the predictive variance
#' with respect to a single input variable for a GP with Gaussian observation
#' model.
#' @param d which input variable.
#' @param K Kernel matrix (N2 x N1).
#' @param dK Array (N2 x N1 x D) with the first derivatives of a Gaussian process
#' kernel matrix with respect to all D inputs.
#' @param WB_mat (K_X_X + sigma^2 I)^(-1)
#' @return Returns a vector with the derivatives of the predictive variance.
dV <- function(d,K,dK,WB_mat) {
  diag(dK[,,d] - t(dK[,,d]) %*% WB_mat %*% K + K %*% WB_mat %*% dK[,,d])
}

#' @description A function for computing the pointwise KL-diff values for a
#' single input variable for a GP with Gaussian observation model.
#' @param d which input variable is evaluated.
#' @param dE A vector with the pointwise derivatives of the mean prediction.
#' @param dV A vector with the pointwise derivatives of the predictive variance.
#' @param V_y A vector with the pointwise predictive variances.
#' @return Returns a vector with the pointwise KL-diff values for input variable d.
KL_diff_gaussian_i <- function(d,dE,dV,V_y) {
  sqrt(dE[,d] * dE[,d] / V_y + 0.5 * dV[,d] * dV[,d] / (V_y * V_y))
}

#' @description A function for computing the second derivative of the mean
#' prediction with respect to a two input variables for a GP with Gaussian
#' observation model.
#' @param d First input variable.
#' @param e Second input variable.
#' @param d2K Array (N2 x N1 x D x D) with the second derivatives of a
#' Gaussian process kernel matrix with respect to all D inputs.
#' @param WB_vec (K_X_X + sigma^2 I)^(-1) * y
#' @return Returns a vector with the derivatives of the mean prediction.
d2E <- function(d,e,d2K,WB_vec) {
  t(d2K[,,d,e]) %*% WB_vec
}

#' @description A function for computing the second derivative of the predictive
#' variance with respect to two input variables for a GP with Gaussian observation
#' model.
#' @param d First input variable.
#' @param e Second input variable.
#' @param K Kernel matrix (N2 x N1).
#' @param dK Array (N2 x N1 x D) with the first derivatives of a Gaussian process
#' kernel matrix with respect to all D inputs.
#' @param d2K Array (N2 x N1 x D x D) with the second derivatives of a
#' Gaussian process kernel matrix with respect to all D inputs.
#' @param WB_mat (K_X_X + sigma^2 I)^(-1)
#' @return Returns a vector with the derivatives of the predictive variance.
d2V <- function(d,e,K,dK,d2K,WB_mat) {
  diag(d2K[,,d,e] - t(d2K[,,d,e]) %*% WB_mat %*% K
       + K %*% WB_mat %*% d2K[,,d,e]
       - t(dK[,,d]) %*% WB_mat %*% dK[,,e]
       - t(dK[,,e]) %*% WB_mat %*% dK[,,d])
}

#' @description A function for computing the pointwise KL-diff2 values for two
#' input variables for a GP with Gaussian observation model.
#' @param d First input variable.
#' @param e Second input variable.
#' @param d2E A vector with the pointwise second derivatives of the mean prediction.
#' @param d2V A vector with the pointwise second derivatives of the predictive variance.
#' @param V_y A vector with the pointwise predictive variances.
#' @return Returns a vector with the pointwise KL-diff2 values for
#' input variables d and e.
KL_diff2_gaussian_ij <- function(d,e,d2E,d2V,V_y) {
  sqrt(d2E[,d,e] * d2E[,d,e] / V_y
       + 0.5 * d2V[,d,e] * d2V[,d,e] / (V_y * V_y))
}


KL_diff2_bernoulli_ij <- function(d,e,E_f,V_f,dE_f,dV_f,d2E_f,d2V_f) {


  EV <- E_f/sqrt(1 + V_f)
  V1 <- sqrt(1 + V_f)
  pis <- pnorm(EV)
  N <- dnorm(EV)

  term1 <- dE_f[,d]/V1 - 0.5 * EV/(V1 * V1) * dV_f[,d]
  term2 <- dE_f[,e]/V1 - 0.5 * EV/(V1 * V1) * dV_f[,e]

  abs(  -EV * N * term1 * term2 + N * (d2E_f[,d,e]/V1 - dE_f[,d] * dV_f[,e] * 0.5/(V1 * V1 * V1) - d2V_f[,d,e] * E_f * 0.5/(V1*V1*V1) - dV_f[,d]*(dE_f[,e]*0.5/(V1*V1*V1) - 0.75*E_f*dV_f[,e]/(V1*V1*V1*V1*V1) ) ) )*sqrt(pis*(1.0 - pis))


}

# ------------------------------
# internal helper functions

# computes the mean of a multidimensional array a
# with respect to dimension i
means.along <- function(a, i) {
  n <- length(dim(a))
  b <- aperm(a, c(seq_len(n)[-i], i))
  rowMeans(b, dims = n - 1)
}
