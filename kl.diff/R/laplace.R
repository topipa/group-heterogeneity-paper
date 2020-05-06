library(optimx)
library(mvtnorm)
#library(invgamma)

#' @export
#' @description A function for computing the squared exponential kernel during
#' laplace approximation.
#' @param X Matrix of observation inputs (N1 x D).
#' @param Z Matrix of observation inputs (N2 x D).
#' @param pars hyperparameters of the GP in logarithms. 1st parameter is the bias
#' term, 2nd parameter is the magnitude term, and the rest are length scales.
#' All are in non-squared form.
se_kernel <- function(X, Z, pars) {
  D <- ncol(X)
  N <- nrow(X)
  M <- nrow(Z)
  if (length(pars) != (D + 2)) {
    warning('wrong number of hyperparameters')
    break
  }

  pars <- exp(pars)
  b <- pars[1]
  alpha <- pars[2]
  ls <- pars[-c(1, 2)]

  X_ <- t(t(X) / ls)
  Z_ <- t(t(Z) / ls)

  X2 <- matrix(apply(X_ * X_, 1, sum), N, M)
  Z2 <- t(matrix(apply(Z_ * Z_, 1, sum), M, N))

  XZ <- X_ %*% t(Z_)

  inner <- (X2 - 2 * XZ + Z2)
  return(b^2 + alpha^2 * exp(-1/2 * inner))
}

#' @export
#' @description Derivative of squared exponential kernel with respect to
#' bias hyperparameter
#' @param X Matrix of observation inputs (N1 x D).
#' @param Z Matrix of observation inputs (N2 x D).
#' @param pars hyperparameters of the GP in logarithms. 1st parameter is the bias
#' term, 2nd parameter is the magnitude term, and the rest are length scales.
#' All are in non-squared form.
se_kernel_d_bias <- function(X, Z, pars) {
  N <- nrow(X)
  M <- nrow(Z)
  b <- exp(pars[1])
  return(matrix(2 * b,N,M))
}

#' @export
#' @description Derivative of squared exponential kernel with respect to
#' magnitude hyperparameter
#' @param X Matrix of observation inputs (N1 x D).
#' @param Z Matrix of observation inputs (N2 x D).
#' @param pars hyperparameters of the GP in logarithms. 1st parameter is the bias
#' term, 2nd parameter is the magnitude term, and the rest are length scales.
#' All are in non-squared form.
se_kernel_d_alpha <- function(X, Z, pars) {
  N <- nrow(X)
  M <- nrow(Z)
  pars <- exp(pars)
  alpha <- pars[2]
  ls <- pars[-c(1, 2)]

  X_ <- t(t(X) / ls)
  Z_ <- t(t(Z) / ls)

  X2 <- matrix(apply(X_ * X_, 1, sum), N, M)
  Z2 <- t(matrix(apply(Z_ * Z_, 1, sum), M, N))

  XZ <- X_ %*% t(Z_)

  inner <- (X2 - 2 * XZ + Z2)
  return(2 * alpha * exp(-1/2 * inner))
}

#' @export
#' @description Derivative of squared exponential kernel with respect to
#' length scale hyperparameters
#' @param X Matrix of observation inputs (N1 x D).
#' @param Z Matrix of observation inputs (N2 x D).
#' @param pars hyperparameters of the GP in logarithms. 1st parameter is the bias
#' term, 2nd parameter is the magnitude term, and the rest are length scales.
#' All are in non-squared form.
se_kernel_d_ls <- function(X, Z, pars, dim) {
  N <- nrow(X)
  M <- nrow(Z)
  pars <- exp(pars)
  alpha <- pars[2]
  ls <- pars[-c(1, 2)]

  X_ <- t(t(X) / ls)
  Z_ <- t(t(Z) / ls)

  X2 <- matrix(apply(X_ * X_, 1, sum), N, M)
  Z2 <- t(matrix(apply(Z_ * Z_, 1, sum), M, N))

  XZ <- X_ %*% t(Z_)

  inner <- (X2 - 2 * XZ + Z2)
  return(alpha^2 * exp(-1/2 * inner) *
         ((matrix(rep(X[,dim],M),N,M)
           - t(matrix(rep(Z[,dim],N),M,N)))^2 / ls[dim]^3))
}

log_lik <- function(y, f, fam) {
  return(log(pnorm(y * f)))
}

log_prior <- function(f, K) {
  return(dmvnorm(f, sigma=K, log=TRUE))
}

log_post <- function(y, X, f, pars, fam, kernel) {
  return(sum(log_lik(y, f, fam)) + log_prior(f, kernel(X, X, pars)))
}

log_marginal <- function(y, K, f, W, fam, a) {
  N <- ncol(W)
  L <- chol(diag(N) + sqrt(W) %*% K %*% sqrt(W))
  fKf <- as.numeric(-0.5 * t(a) %*% f)
  ll <- sum(log_lik(y, f, fam))
  B <- sum(log(diag(L)))
  return(fKf + ll - B)
}

# find_mode <- function(y, X, f, pars, fam, kernel) {
#   log_post_optim <- function(f) {
#     return(log_post(y, X, f, pars, fam, kernel))
#   }
#   sol <- optim(f, log_post_optim, method="L-BFGS-B",
#                hessian=TRUE, control=list(fnscale=-1, maxit=10))
#   return(list(f=sol$par, W=sol$hessian))
# }

#' @export
#' @description Function for optimizing latent values f in Laplace approximation.
#' @param y vector of outputs (-1 or 1)
#' @param X Matrix of observation inputs (N1 x D).
#' @param pars hyperparameters of the GP in logarithms. 1st parameter is the bias
#' term, 2nd parameter is the magnitude term, and the rest are length scales.
#' All are in non-squared form.
#' @param kernel kernel function, currently only se_kernel is supported
#' @param maxit maximum number of iterations
#' @param tol stopping tolerance
find_mode_manual <- function(y, X, f, pars, fam, kernel, maxit = 50, tol = 1e-8) {
  N <- nrow(X)
  K <- kernel(X, X, pars)
  obj <- -10000000
  for (i in seq(maxit)) {
    W <- (dnorm(f)/pnorm(y*f))^2 + y * f * dnorm(f)/pnorm(y*f)
    W <- diag(W)
    L <- chol(diag(N) + sqrt(W) %*% K %*% sqrt(W))
    L <- t(L)
    grad <- y * dnorm(f)/pnorm(y*f)
    b <- W %*% f + grad
    a <- sqrt(W) %*% K %*% b
    a <- solve(L,a)
    a <- solve((t(L)),a)
    a <- b - sqrt(W) %*% a
    f <- c(K %*% a)

    objective <- -0.5 * t(a) %*% f + sum(log(pnorm(y * f)))
    if (abs(obj - objective) < tol) {
      obj <- objective
      break
    }
    obj <- objective
  }
  list(f = f, W = W, a = a)
}

# optim_hyper <- function(y, X, f, pars, fam, kernel, maxit=100) {
#   D <- ncol(X)
#   for (i in 1:maxit) {
#     #print(i)
#     #print(f)
#     mode <- find_mode_manual(y, X, f, pars, fam, kernel)
#     f <- mode$f
#     W <- mode$W
#     logZ <- function(pars) {
#       K <- kernel(X, X, pars)
#       ## log_reg_pars <- dinvgamma(pars[-c(1, 2)], shape=1, scale=1, log=TRUE)
#       log_ml <- log_marginal(y, K, f, W, fam, a)
#       return(log_ml)
#     }
#     opt <- optim(pars, logZ, method="L-BFGS-B",
#                  control=list(fnscale=-1, maxit=1))
#     if (mean((pars - opt$par)^2) <= 1e-6) {
#       pars <- opt$par
#       break
#     }
#     pars <- opt$par
#     #print(pars)
#   }
#   mode <- find_mode_manual(y, X, f, pars, fam, kernel)
#   f <- mode$f
#   return(list(f=f, pars=pars))
# }

dlogZ_par <- function(y, X, f, K, s2, a, R, kernel_d_par, pars) {
  C <- kernel_d_par(X,X,pars)
  s1 <- 0.5 * t(a) %*% C %*% a - 0.5 * sum(diag(R %*% C))
  b <- C %*% (y * dnorm(f)/pnorm(y*f))
  s3 <- b - K %*% R %*% b
  return(s1 + t(s2) %*% s3)
}


#' @export
#' @description Function for optimizing GP hyperparameters in Laplace approximation.
#' @param y vector of outputs (-1 or 1)
#' @param X Matrix of observation inputs (N1 x D).
#' @param pars hyperparameters of the GP in NON-LOG form. 1st parameter is the bias
#' term, 2nd parameter is the magnitude term, and the rest are length scales.
#' All are in non-squared form.
#' @param kernel kernel function, currently only se_kernel is supported
#' @param maxit maximum number of iterations
#' @param stoptol stopping tolerance
#' @param optimizer which optimizer is used. homemade seems to work best
#' @param lr learning rate. Only used for 'homemade' optimizer
#' @return optimized latent values f and hyperparameters in non-log form
optim_hyper_grad <- function(y, X, f, pars, fam, kernel, maxit=200, stoptol = 1e-9, optimizer = 'homemade', lr=5e-2) {
  pars <- log(pars)
  D <- ncol(X)
  N <-  nrow(X)
  lml.last <- -999
  for (i in 1:maxit) {
    mode <- find_mode_manual(y, X, f, pars, fam, kernel)
    f <- mode$f
    W <- mode$W
    a <- mode$a

    logZ <- function(pars) {
      K <- kernel(X, X, pars)
      log_ml <- log_marginal(y, K, f, W, fam, a)
      return(log_ml)
    }

    dlogZ <- function(pars) {
      K <- kernel(X, X, pars)
      L <- chol(diag(N) + sqrt(W) %*% K %*% sqrt(W))
      L <- t(L)

      R <- sqrt(W)
      R <- solve(L,R)
      R <- solve(t(L),R)
      R <- sqrt(W) %*% R

      C <- solve(L,(sqrt(W) %*% K))

      nn <- dnorm(f)
      tt <- pnorm(f*y)
      deriv3 <- 2 * f * nn * nn / (tt * tt) +
        2 * y * (nn / tt)^3 - y * nn / tt + y * f * f * nn / tt + f * (nn / tt)^2

      s2 <- -0.5 * diag(diag(K) - diag(t(C) %*% C)) %*% deriv3

      npars <- length(pars)
      ret <- rep(0,npars)

      ret[1] <- dlogZ_par(y, X, f, K, s2, a, R, se_kernel_d_bias, pars)
      ret[2] <- dlogZ_par(y, X, f, K, s2, a, R, se_kernel_d_alpha, pars)

      for (npar in 3:npars)
        ret[npar] <- dlogZ_par(y, X, f, K, s2, a, R, function(X1, X2, pars)
          se_kernel_d_ls(X1, X2, pars, npar - 2),
          pars)

      return(ret)
    }

    if (optimizer == 'optim') {
      opt <- optim(pars, logZ, method = "L-BFGS-B",
                   gr = dlogZ,
                   control = list(fnscale = -1, maxit = 1))

      pars <- opt$par
      print(pars)

      lml <- logZ(pars)
      print(lml)
      if (abs((lml - lml.last)^2 / lml.last^2) <= stoptol) {
        break
      }
      lml.last <- lml
    }
    else if (optimizer == 'optimx') {
      opt <- optimx(pars, logZ,  gr = dlogZ,
                    method = "L-BFGS-B",
                    itnmax = 1,
                    control = list(maximize=TRUE, starttests = FALSE, kkt = FALSE, maxit = 1))

      best <- which(opt$value == max(opt$value))
      pars <- unlist(lapply(1:length(pars), function(i) opt[[paste0("p", i)]][best]))
      print(pars)
      lml <- logZ(pars)
      print(lml)

      if (abs((lml - lml.last)^2 / lml.last^2) <= stoptol) {
        break
      }
      lml.last <- lml
    }
    else if (optimizer == 'homemade') {
      grad <- dlogZ(pars)
      pars <- pars + lr * grad
      # print(pars)
      lml <- logZ(pars)

      # print(lml)

      if (abs((lml - lml.last)^2 / lml.last^2) <= stoptol) {
        break
      }
      lml.last <- lml
    }
    else {
      warning('unknown optimizer')
      break
    }

  }

  mode <- find_mode_manual(y, X, f, pars, fam, kernel)
  f <- mode$f
  return(list(f=f, pars=exp(pars)))
}

