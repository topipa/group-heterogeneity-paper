

# generate toy data
set.seed(124)
n <- 500
n = as.integer(n)
k = as.integer(3)
x = mvtnorm::rmvnorm(n = n,mean = rep(0.0,k),sigma = diag(rep(1.0,k)))
w = c(0,1,0)
noise = rnorm(n = n, mean = 0.0,sd = 1)
y = rowSums(sweep(x,MARGIN=2,w,`*`))
y <- y + x[,1]*x[,3] + noise

normalize_matrix <- function(a) {
  b = sweep(a,MARGIN=2,(apply(a,MARGIN = 2,FUN=mean)),`-`)
  return(sweep(b,MARGIN=2,(apply(b,MARGIN = 2,FUN=sd)),`/`))
}

normalize_vector <- function(a) {
  b = a - mean(a)
  return(b/sd(b))
}

x = normalize_matrix(x)
y = normalize_vector(y)



sigma <- 1
sigma <- sigma * sigma

alpha <- 1.0
ls <- 1.5


KL_diff_gaussian(y,x,alpha,ls,sigma,pointwise = FALSE)

KL_diff2_gaussian(y,x,alpha,ls,sigma,pointwise = FALSE)





KLdiff_pw <- KL_diff_gaussian(y,x,alpha,ls,sigma,pointwise = TRUE)
plot(x[,1],y)
points(x[,1],KLdiff_pw[,1],col = 'red')


plot(x[,2],y)
points(x[,2],KLdiff_pw[,2],col = 'red')

plot(x[,3],y)
points(x[,3],KLdiff_pw[,3],col = 'red')



