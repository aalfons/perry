context("perryFit - repeated CV")


## load packages
library("perry", quietly=TRUE)

## set seed for reproducibility
set.seed(1234)

## generate data for tests
n <- 20
x <- rnorm(n)
y <- x + rnorm(n)
yy <- cbind(y1=y, y2=y+rnorm(n))
x <- as.matrix(x)
xy <- data.frame(x, y)

## fit univariate models via lmrob() and lts()
lmrobFit <- lmrob(y~x, data=xy)
ltsFit <- ltsReg(x, y, alpha=0.75)
ltsFit$call[[1]] <- as.name("ltsReg")

## fit multivariate models via lm()
lmFit <- lm(yy~x)

## set up cross-validation folds
K <- 5
R <- 2
folds <- cvFolds(n, K, R)


## run tests

test_that("univariate response yields correct \"perry\" object", {
        ## MM-regression
        lmrobCV <- perryFit(lmrobFit, data=xy, y=xy$y, 
            splits=folds, cost=rtmspe)
        
        expect_is(lmrobCV, "perry")
        lmrobPE <- lmrobCV$pe
        expect_is(lmrobPE, "numeric")
        expect_equal(length(lmrobPE), 1)
        lmrobReps <- lmrobCV$reps
        expect_is(lmrobReps, "matrix")
        expect_equal(dim(lmrobReps), c(R, 1))
        lmrobSE <- lmrobCV$se
        expect_is(lmrobSE, "numeric")
        expect_equal(length(lmrobSE), 1)
        expect_false(any(is.na(lmrobSE)))
        
        ## reweighted and raw LTS
        ltsCV <- perryFit(ltsFit, x=x, y=y, splits=folds, 
            predictArgs=list(fit="both"), cost=rtmspe)
        
        expect_is(ltsCV, "perry")
        ltsPE <- ltsCV$pe
        expect_is(ltsPE, "numeric")
        expect_equal(length(ltsPE), 2)
        expect_equal(ltsFit$reps, NULL)
        ltsReps <- ltsCV$reps
        expect_is(ltsReps, "matrix")
        expect_equal(dim(ltsReps), c(R, 2))
        ltsSE <- ltsCV$se
        expect_is(ltsSE, "numeric")
        expect_equal(length(ltsSE), 2)
        expect_false(any(is.na(ltsSE)))
    })

test_that("multivariate response yields correct \"perry\" object", {
        ## multivariate LS regression
        lmCV <- perryFit(lmFit, data=lmFit$model, y=yy, splits=folds)
        
        expect_is(lmCV, "perry")
        lmPE <- lmCV$pe
        expect_is(lmPE, "numeric")
        expect_equal(length(lmPE), 1)
        lmReps <- lmCV$reps
        expect_is(lmReps, "matrix")
        expect_equal(dim(lmReps), c(R, 1))
        lmSE <- lmCV$se
        expect_is(lmSE, "numeric")
        expect_equal(length(lmSE), 1)
        expect_false(any(is.na(lmSE)))
    })
