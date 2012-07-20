# load data and fit an LS regression model
data("mtcars")
fit <- lm(mpg ~ wt + cyl, data=mtcars)

# perform cross-validation
repCV(fit, K = 5, R = 10, seed = 1234)  # K-fold CV
repCV(fit, K = nrow(mtcars))            # leave-one-out CV
