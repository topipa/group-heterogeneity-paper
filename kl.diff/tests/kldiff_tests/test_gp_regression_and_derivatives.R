library(Rcpp)
library(kl.diff)




# generate toy data
set.seed(124)
n <- 100
n = as.integer(n)
k = as.integer(1)
x = mvtnorm::rmvnorm(n = n,mean = rep(0.0,k),sigma = diag(rep(1.0,k)))
w = c(1)
noise = rnorm(n = n, mean = 0.0,sd = 1)
y = rowSums(sweep(x,MARGIN=2,w,`*`)) + noise

x = normalize_matrix(x)
y = normalize_vector(y)


#lm = lm(y~ x)


sigma <- 1
sigma <- sigma * sigma

alpha <- 1
ls <- 1

xs <- matrix(seq(-4,4,length.out = 250),250,1)
k_x_x <- squared_exponential_cov(x,x,alpha,ls)
k_x_xs <- squared_exponential_cov(x,xs,alpha,ls)
# DO NOT COMPUTE THIS
# USE t(k_x_xs)
#k_xs_x <- squared_exponential_cov(xs,x,alpha,ls)
k_xs_xs <- squared_exponential_cov(xs,xs,alpha,ls)

E_f <- t(k_x_xs$K) %*% solve(k_x_x$K + diag(n)*sigma,y)
V_f <- k_xs_xs$K - t(k_x_xs$K) %*% solve(k_x_x$K + diag(n)*sigma, k_x_xs$K)
V_f = diag(V_f)
V_y <- V_f + sigma


plot(x,y, xlim = c(-3.7,3.7), ylim = c(-3,3))
points(xs,E_f,col = 'red')
points(xs,E_f + 2* sqrt(V_f),col = 'blue')
points(xs,E_f - 2* sqrt(V_f),col = 'blue')
points(xs,E_f + 2* sqrt(V_y),col = 'blue')
points(xs,E_f - 2* sqrt(V_y),col = 'blue')

plot(x,y)
ddE_f <- (E_f[2:250,1] - E_f[1:249,1])/(8/250)
xxs <- matrix(seq(-4,4,length.out = 249),249,1)
points(xxs,ddE_f,col = 'green')

dE_f <- t(k_x_xs$dK[,,1]) %*% solve(k_x_x$K + diag(n)*sigma,y)
points(xs,dE_f,col = 'green')


dd2E_f <- (dE_f[2:249,1] - dE_f[1:248,1])/(8/249)
xxxs <- matrix(seq(-4,4,length.out = 248),248,1)
points(xxxs,dd2E_f,col = 'orange')


d2E_f <- t(k_x_xs$d2K[,,1,1]) %*% solve(k_x_x$K + diag(n)*sigma,y)
points(xs,d2E_f,col = 'red')

