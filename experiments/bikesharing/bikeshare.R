library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
library(kl.diff)
source("../helper_simulations.R")


bikedata = read.csv("data/day.csv", header = TRUE)

SEED <- 124
set.seed(SEED)

gp_model2 <- stan_model("../gaussian_simulated/gaussian_gp_ard_constrained.stan")
gp_model4 <- stan_model("../gaussian_simulated/gaussian_gp_linear_and_ard_priors.stan")



y <- bikedata$cnt
x <- data.frame(temp = bikedata$temp, wind = bikedata$windspeed, hum = bikedata$hum,
                weather = bikedata$weathersit, holiday = bikedata$holiday,
                weekday = bikedata$weekday, month = bikedata$mnth,
                season = bikedata$season)

D <- 3
Ddummy <- 5
x[,1:D] <- normalize_matrix(x[,1:D])
y <- normalize_vector(y)
N <- length(y)

gp_data2 <- list(x = as.matrix(x), y = y, N = N, D = D, Ddummy = Ddummy)
opt_fit2 <- optimizing(gp_model2, data = gp_data2, seed = SEED,
                       hessian = FALSE)

sigma <- opt_fit2$par[D + Ddummy + 2]
sigma <- sigma * sigma
alpha <- opt_fit2$par[D + Ddummy + 1]
ls <- sqrt(opt_fit2$par[1:(D + Ddummy)])

rankints <- rank_interactions_gaussian(y,as.matrix(x),NULL,alpha,ls,sigma,pointwise = FALSE)
#rankintsb <- rank_interactions_gaussian_linear(y,as.matrix(x), Xlin = as.matrix(x[,1:3]),NULL,alpha,ls,sigma,pointwise = FALSE)
#kld1 <- rank_variables_gaussian(y,as.matrix(x),alpha,ls,sigma,pointwise = FALSE)


# [,1]      [,2]       [,3]      [,4]       [,5]      [,6]       [,7]       [,8]
# [1,]    0 0.1185256 0.28700710 0.2393158 0.08069783 0.2475762 0.32704364 0.17890637
# [2,]    0 0.0000000 0.09264108 0.0922663 0.02885629 0.1122320 0.12242904 0.04072951
# [3,]    0 0.0000000 0.00000000 0.2092002 0.06779481 0.2316296 0.26760411 0.12386179
# [4,]    0 0.0000000 0.00000000 0.0000000 0.06281771 0.2314011 0.18047486 0.10558779
# [5,]    0 0.0000000 0.00000000 0.0000000 0.00000000 0.1239234 0.08607992 0.01419978
# [6,]    0 0.0000000 0.00000000 0.0000000 0.00000000 0.0000000 0.22304086 0.09930635
# [7,]    0 0.0000000 0.00000000 0.0000000 0.00000000 0.0000000 0.00000000 0.19487921
# [8,]    0 0.0000000 0.00000000 0.0000000 0.00000000 0.0000000 0.00000000 0.00000000


rowSums(kld2[1:3,4:8])
colSums(kld2[1:3,4:8])
colnames(x)




df <- x
df$y <- y

library(rstanarm)


# form a base model

model4 <- stan_glmer(y ~ wind + temp + hum + (1 | month) + (1 | weekday) + (1 | weather) + (1 | holiday) + (1 | season), data = df)
loo4 <- loo(model4)

# then add single varying slopes

model51 <- stan_glmer(y ~ wind + temp + hum + (1 + temp | month) + (1 | weekday) + (1 | weather) + (1 | holiday) + (1 | season), data = df)
loo51 <- loo(model51)
model52 <- stan_glmer(y ~ wind + temp + hum + (1 + hum | month) + (1 | weekday) + (1 | weather) + (1 | holiday) + (1 | season), data = df)
loo52 <- loo(model52)
model53 <- stan_glmer(y ~ wind + temp + hum + (1 + wind | month) + (1 | weekday) + (1 | weather) + (1 | holiday) + (1 | season), data = df)
loo53 <- loo(model53)

loo4$estimates[1,1]
loo51$estimates[1,1]
loo52$estimates[1,1]
loo53$estimates[1,1]

loo51$estimates[1,1] - loo4$estimates[1,1]
loo52$estimates[1,1] - loo4$estimates[1,1]
loo53$estimates[1,1] - loo4$estimates[1,1]

sd(model51$coefficients[c(6,8,10,12,14,16,18,20,22,24,26,28)])
sd(model52$coefficients[c(6,8,10,12,14,16,18,20,22,24,26,28)])
sd(model53$coefficients[c(6,8,10,12,14,16,18,20,22,24,26,28)])

model54 <- stan_glmer(y ~ wind + temp + hum + (1 | month) + (1 + temp | weekday) + (1 | weather) + (1 | holiday) + (1 | season), data = df)
loo54 <- loo(model54)
model55 <- stan_glmer(y ~ wind + temp + hum + (1 | month) + (1 + hum | weekday) + (1 | weather) + (1 | holiday) + (1 | season), data = df)
loo55 <- loo(model55)
model56 <- stan_glmer(y ~ wind + temp + hum + (1 | month) + (1 + wind | weekday) + (1 | weather) + (1 | holiday) + (1 | season), data = df)
loo56 <- loo(model56)

loo54$estimates[1,1]
loo55$estimates[1,1]
loo56$estimates[1,1]

sd(model54$coefficients[c(18,20,22,24,26,28,30)])
sd(model55$coefficients[c(18,20,22,24,26,28,30)])
sd(model56$coefficients[c(18,20,22,24,26,28,30)])



model57 <- stan_glmer(y ~ wind + temp + hum + (1 | month) + (1 | weekday) + (1 + temp | weather) + (1 | holiday) + (1 | season), data = df)
loo57 <- loo(model57)
model58 <- stan_glmer(y ~ wind + temp + hum + (1 | month) + (1 | weekday) + (1 + hum | weather) + (1 | holiday) + (1 | season), data = df)
loo58 <- loo(model58)
model59 <- stan_glmer(y ~ wind + temp + hum + (1 | month) + (1 | weekday) + (1 + wind | weather) + (1 | holiday) + (1 | season), data = df)
loo59 <- loo(model59)

loo57$estimates[1,1]
loo58$estimates[1,1]
loo59$estimates[1,1]

sd(model57$coefficients[c(29,31,33)])
sd(model58$coefficients[c(29,31,33)])
sd(model59$coefficients[c(29,31,33)])
