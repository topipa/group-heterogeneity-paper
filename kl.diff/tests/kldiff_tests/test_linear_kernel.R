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

bias <- 10
dev <- 1

xs <- matrix(seq(-6,6,length.out = 250),250,1)
k_x_x <- linear_cov(x,x,bias,dev)
k_x_xs <- linear_cov(x,xs,bias,dev)
# DO NOT COMPUTE THIS
# USE t(k_x_xs)
#k_xs_x <- squared_exponential_cov(xs,x,alpha,ls)
k_xs_xs <- linear_cov(xs,xs,bias,dev)

E_f <- t(k_x_xs$K) %*% solve(k_x_x$K + diag(n)*sigma,y)
V_f <- k_xs_xs$K - t(k_x_xs$K) %*% solve(k_x_x$K + diag(n)*sigma, k_x_xs$K)
V_f = diag(V_f)
V_y <- V_f + sigma


plot(x,y, xlim = c(-6,6), ylim = c(-5,5))
points(xs,E_f,col = 'red')
points(xs,E_f + 2* sqrt(V_f),col = 'blue')
points(xs,E_f - 2* sqrt(V_f),col = 'blue')
points(xs,E_f + 2* sqrt(V_y),col = 'blue')
points(xs,E_f - 2* sqrt(V_y),col = 'blue')

plot(x,y)
ddE_f <- (E_f[2:250,1] - E_f[1:249,1])/(12/250)
xxs <- matrix(seq(-4,4,length.out = 249),249,1)
points(xxs,ddE_f,col = 'green')

dE_f <- t(k_x_xs$dK[,,1]) %*% solve(k_x_x$K + diag(n)*sigma,y)
points(xs,dE_f,col = 'red')


dd2E_f <- (dE_f[2:249,1] - dE_f[1:248,1])/(12/249)
xxxs <- matrix(seq(-4,4,length.out = 248),248,1)
points(xxxs,dd2E_f,col = 'orange')


d2E_f <- t(k_x_xs$d2K[,,1,1]) %*% solve(k_x_x$K + diag(n)*sigma,y)
points(xs,d2E_f,col = 'red')






# test sum of lin and se kernels



# generate toy data
set.seed(124)
n <- 100
n = as.integer(n)
k = as.integer(1)
x = mvtnorm::rmvnorm(n = n,mean = rep(0.0,k),sigma = diag(rep(1.0,k)))
w = c(1)
noise = rnorm(n = n, mean = 0.0,sd = 0.3)
y = rowSums(sweep(x,MARGIN=2,w,`*`)) + sin(5*x) + noise


x = normalize_matrix(x)
y = normalize_vector(y)


# sigma <- 1
# sigma <- sigma * sigma
# bias <- 0
# dev <- 1.5
# alpha <- 2
# ls <- 0.4


# optimized hyperparams
sigma <- 0.25
sigma <- sigma * sigma
bias <- 0
dev <- 0.75
alpha <- 0.6
ls <- 0.346




xs <- matrix(seq(-6,6,length.out = 250),250,1)
k_x_x <- linear_cov(x,x,bias,dev) 
k_x_x2 <- squared_exponential_cov(x,x,alpha,ls)
k_x_x$K <-k_x_x$K + k_x_x2$K
k_x_x$dK <-k_x_x$dK+ k_x_x2$dK
k_x_x$d2K <-k_x_x$d2K + k_x_x2$d2K

k_x_xs <- linear_cov(x,xs,bias,dev)
k_x_xs2 <- squared_exponential_cov(x,xs,alpha,ls)
k_x_xs$K <-k_x_xs$K + k_x_xs2$K
k_x_xs$dK <-k_x_xs$dK+ k_x_xs2$dK
k_x_xs$d2K <-k_x_xs$d2K + k_x_xs2$d2K
# DO NOT COMPUTE THIS
# USE t(k_x_xs)
#k_xs_x <- squared_exponential_cov(xs,x,alpha,ls)
k_xs_xs <- linear_cov(xs,xs,bias,dev)
k_xs_xs2 <- squared_exponential_cov(xs,xs,alpha,ls)
k_xs_xs$K <-k_xs_xs$K + k_xs_xs2$K
k_xs_xs$dK <-k_xs_xs$dK+ k_xs_xs2$dK
k_xs_xs$d2K <-k_xs_xs$d2K + k_xs_xs2$d2K

E_f <- t(k_x_xs$K) %*% solve(k_x_x$K + diag(n)*sigma,y)
V_f <- k_xs_xs$K - t(k_x_xs$K) %*% solve(k_x_x$K + diag(n)*sigma, k_x_xs$K)
V_f = diag(V_f)
V_y <- V_f + sigma


plot(x,y, xlim = c(-6,6), ylim = c(-15,15))
points(xs,E_f,col = 'red')
points(xs,E_f + 2* sqrt(V_f),col = 'blue')
points(xs,E_f - 2* sqrt(V_f),col = 'blue')
points(xs,E_f + 2* sqrt(V_y),col = 'blue')
points(xs,E_f - 2* sqrt(V_y),col = 'blue')






