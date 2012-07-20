# load data and fit an LS regression model
data("mtcars")
fit <- lm(mpg ~ wt + cyl, data=mtcars)

# perform bootstrap prediction error estimation
bootPE(fit, R = 10, bootType = "0.632", seed = 1234)
bootPE(fit, R = 10, bootType = "out-of-bag", seed = 1234)
