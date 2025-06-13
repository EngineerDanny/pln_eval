library(glmnet)
set.seed(2020)
n <- 100
p <- 4
x <- matrix(runif(n * p, 5, 10), n)
y <- rpois(n, exp(rowMeans(x)))
# glm fit
glmfit <- glm(y ~ x - 1, family = poisson)
coef(glmfit)
