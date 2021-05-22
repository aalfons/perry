library("perryExamples")
data("coleman")
set.seed(1234)  # set seed for reproducibility

## set up folds for cross-validation
folds <- cvFolds(nrow(coleman), K = 5, R = 10)

## compare LS, MM and LTS regression

# perform cross-validation for an LS regression model
fitLm <- lm(Y ~ ., data = coleman)
cvLm <- perry(fitLm, splits = folds,
              cost = rtmspe, trim = 0.1)

# perform cross-validation for an MM regression model
fitLmrob <- lmrob(Y ~ ., data = coleman, maxit.scale = 500)
cvLmrob <- perry(fitLmrob, splits = folds,
                 cost = rtmspe, trim = 0.1)

# perform cross-validation for an LTS regression model
fitLts <- ltsReg(Y ~ ., data = coleman)
cvLts <- perry(fitLts, splits = folds,
               cost = rtmspe, trim = 0.1)

# combine results into one object
cv <- perrySelect(LS = cvLm, MM = cvLmrob, LTS = cvLts,
                  .selectBest = "min")
cv

# plot results for the MM regression model
plot(cvLmrob, which = "box")
plot(cvLmrob, which = "density")
plot(cvLmrob, which = "dot", seFactor = 2)

# plot combined results
plot(cv, which = "box")
plot(cv, which = "density")
plot(cv, which = "dot", seFactor = 2)
