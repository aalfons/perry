# load data and fit an LS regression model
data("mtcars")
fit <- lm(mpg ~ wt + cyl, data=mtcars)

# perform random splitting
repRS(fit, m = 6, R = 10, seed = 1234)
