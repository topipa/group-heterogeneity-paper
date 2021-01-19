
library(rstan)
library(kl.diff)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
source("../helper_simulations.R")


bikedata = read.csv('data/day.csv', header = TRUE)

y <- bikedata$cnt
x <- data.frame(time = bikedata$instant,
                temp = bikedata$temp, wind = bikedata$windspeed, hum = bikedata$hum,
                weather = bikedata$weathersit, holiday = bikedata$holiday,
                weekday = bikedata$weekday, month = bikedata$mnth,
                season = bikedata$season)
D <- 4
Ddummy <- 5
x[,1:D] <- normalize_matrix(x[,1:D])
y <- normalize_vector(y)
N <- length(y)

gp_model <- stan_model("../stanmodels/gp_ard_constrained.stan")


SEED <- 124

gp_data <- list(x = as.matrix(x), y = y, N = N, D = D, Ddummy = Ddummy)
opt_fit <- optimizing(gp_model, data = gp_data, seed = SEED,
                       hessian = FALSE)

sigma <- opt_fit$par[D + Ddummy + 2]
sigma2 <- sigma * sigma
alpha <- opt_fit$par[D + Ddummy + 1]
ls <- sqrt(opt_fit$par[1:(D + Ddummy)])

kld2 <- KL_diff2_gaussian(y,as.matrix(x),NULL,alpha,ls,sigma2,pointwise = FALSE)



grouping_importance <- colSums(kld2[2:4,5:9])
names(grouping_importance) <- colnames(x)[5:9]
# most important grouping is month

# which predictors vary most by month
month_predictor_importance <- kld2[2:4,8]
names(month_predictor_importance) <- colnames(x)[2:4]
# 1. temperature
# 2. humidity
# 3. wind speed





library(rstanarm)
df <- x
df$y <- y

# form a base model
model1 <- stan_glmer(y ~ wind + temp + hum + time + (1 | month) + (1 | weekday) + (1 | weather) + (1 | holiday) + (1 | season), data = df, chains = 4, seed = SEED, cores = 4)
loo1 <- loo(model1, save_psis = TRUE)
pred1 <- posterior_predict(model1)
plot(y, colMeans(pred1))
points(y,y,col = 'red')


# then add single varying slopes to month indicator
model51 <- stan_glmer(y ~ wind + temp + hum + time + (1 + temp | month) + (1 | weekday) + (1 | weather) + (1 | holiday) + (1 | season), data = df, chains = 4, seed = SEED, cores = 4)
loo51 <- loo(model51, save_psis = TRUE)
model52 <- stan_glmer(y ~ wind + temp + hum + time + (1 + hum | month) + (1 | weekday) + (1 | weather) + (1 | holiday) + (1 | season), data = df, chains = 4, seed = SEED, cores = 4)
loo52 <- loo(model52, save_psis = TRUE)
model53 <- stan_glmer(y ~ wind + temp + hum + time + (1 + wind | month) + (1 | weekday) + (1 | weather) + (1 | holiday) + (1 | season), data = df, chains = 4, seed = SEED, cores = 4)
loo53 <- loo(model53, save_psis = TRUE)


# elpd differences
loo_compare(loo51, loo1)
loo_compare(loo52, loo1)
loo_compare(loo53, loo1)

# slope variability
sd(model51$coefficients[c(6,8,10,12,14,16,18,20,22,24,26,28)+1])
sd(model52$coefficients[c(6,8,10,12,14,16,18,20,22,24,26,28)+1])
sd(model53$coefficients[c(6,8,10,12,14,16,18,20,22,24,26,28)+1])


# MSE
pred1 <- (posterior_predict(model1) - t(matrix(y,N,4000)))^2
pred51 <- (posterior_predict(model51) - t(matrix(y,N,4000)))^2
pred52 <- (posterior_predict(model52) - t(matrix(y,N,4000)))^2
pred53 <- (posterior_predict(model53) - t(matrix(y,N,4000)))^2
mse1 <- loo::E_loo(pred1, loo1$psis_object)
mse51 <- loo::E_loo(pred51, loo51$psis_object)
mse52 <- loo::E_loo(pred52, loo52$psis_object)
mse53 <- loo::E_loo(pred53, loo53$psis_object)

sum(mse51$value - mse1$value)
sum(mse52$value - mse1$value)
sum(mse53$value - mse1$value)

sd(mse51$value - mse1$value) * sqrt(N)
sd(mse52$value - mse1$value) * sqrt(N)
sd(mse53$value - mse1$value) * sqrt(N)


model_final <- stan_glmer(y ~ wind + temp + hum + time + (1 + temp + hum | month) + (1 | weekday) + (1 | weather) + (1 | holiday) + (1 | season), data = df, chains = 4, seed = SEED, cores = 4)
loo_final <- loo(model_final, save_psis = TRUE)
