# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @export
perryTool <- function(call, data = NULL, x = NULL, y, cost = rmspe, splits, 
        names = NULL, predictArgs = list(), costArgs = list(),  
        envir = parent.frame()) {
    UseMethod("perryTool", splits)
}

#' @S3method perryTool cvFolds
perryTool.cvFolds <- function(call, data = NULL, x = NULL, y, cost = rmspe, 
        splits, names = NULL, predictArgs = list(), costArgs = list(), 
        envir = parent.frame()) {
    # define an expression that obtains predictions for each test data subset 
    # in one replication of cross-validation
    if(is.null(data)) {
        if(is.null(names)) names <- c("x", "y")
        cvExpr <- expression(
            lapply(subsets, cvXY, call=call, x=x, y=y, names=names, 
                predictArgs=predictArgs, envir=envir)
        )
    } else {
        if(is.null(names)) names <- "data"
        cvExpr <- expression(
            lapply(subsets, cvData, call=call, data=data, names=names, 
                predictArgs=predictArgs, envir=envir)
        )
    }
    # define a function the evaluates the expression to obtain the prediction 
    # and computes the prediction error for one replication of cross-validation
    cvFun <- function(r, keepSE) {
        subsets <- getIndices(splits, r)
        # instead of collecting the results from the folds in the original 
        # order of the observations, the response is re-ordered accordingly
        yHat <- eval(cvExpr)       # obtain predictions from the folds
        yHat <- combineData(yHat)  # combine predictions from the folds
        y <- dataSubset(y, unlist(subsets))  # re-order response
        # compute cost function for predicted values
        pe <- computeCost(cost, y, yHat, args=costArgs, keepSE=keepSE)
    }
    # perform (repeated) cross-validation
    R <- splits$R
    if(R == 1) {
        pe <- cvFun(R, keepSE=TRUE)
        if(!is.list(pe)) pe <- t(pe)
    } else {
        pe <- lapply(seq_len(R), cvFun, keepSE=FALSE)  # for each replication
        pe <- combineData(pe, drop=FALSE)              # combine results
    }
    # return prediction error (add default names if necessary)
    addNames(pe)
}

#' @S3method perryTool randomSplits
perryTool.randomSplits <- function(call, data = NULL, x = NULL, y, 
        cost = rmspe, splits, names = NULL, predictArgs = list(), 
        costArgs = list(), envir = parent.frame()) {
    # define an expression that obtains predictions for the test data in one 
    # replication of Monte Carlo cross-validation
    if(is.null(data)) {
        if(is.null(names)) names <- c("x", "y")
        cvExpr <- expression(
            cvXY(subset, call=call, x=x, y=y, names=names, 
                predictArgs=predictArgs, envir=envir)
        )
    } else {
        if(is.null(names)) names <- "data"
        cvExpr <- expression(
            cvData(subset, call=call, data=data, names=names, 
                predictArgs=predictArgs, envir=envir)
        )
    }
    # define a function the evaluates the expression to obtain the prediction 
    # and computes the prediction error for one replication of Monte Carlo 
    # cross-validation
    cvFun <- function(r, keepSE) {
        subset <- getIndices(splits, r)
        yHat <- eval(cvExpr)   # obtain predictions for test data
        y <- dataSubset(y, subset)  # response for test data
        # compute cost function for predicted values
        pe <- computeCost(cost, y, yHat, args=costArgs, keepSE=keepSE)
    }
    # perform Monte Carlo cross-validation
    R <- splits$R
    if(R == 1) {
        pe <- cvFun(R, keepSE=TRUE)
        if(!is.list(pe)) pe <- t(pe)
    } else {
        pe <- lapply(seq_len(R), cvFun, keepSE=FALSE)  # for each replication
        pe <- combineData(pe, drop=FALSE)              # combine results
    }
    # return prediction error (add default names if necessary)
    addNames(pe)
}

# fit the model for the training data and obtain predictions of the test data
# for functions that take the predictors and the response as separate arguments
cvXY <- function(i, call, x, y, names, predictArgs, envir) {
    # plug training data into function call
    call[[names[1]]] <- dataSubset(x, -i)
    call[[names[2]]] <- dataSubset(y, -i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predict, fit, dataSubset(x, i), args=predictArgs)
}

# fit the model for the training data and obtain predictions of the test data
# for functions that have one argument for all the data
cvData <- function(i, call, data, names, predictArgs, envir) {
    # plug training data into function call
    call[[names]] <- dataSubset(data, -i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predict, fit, dataSubset(data, i), args=predictArgs)
}

#' @S3method perryTool bootSamples
perryTool.bootSamples <- function(call, data = NULL, x = NULL, y, 
        cost = rmspe, splits, names = NULL, predictArgs = list(), 
        costArgs = list(), envir = parent.frame()) {
    ## compute out-of-bag prediction error
    # define an expression that obtains the out-of-bag predictions for one 
    # bootstrap sample
    if(is.null(data)) {
        if(is.null(names)) names <- c("x", "y")
        bootExpr <- expression(
            bootXY(s, call=call, x=x, y=y, names=names, 
                predictArgs=predictArgs, envir=envir)
        )
    } else {
        if(is.null(names)) names <- "data"
        bootExpr <- expression(
            bootData(s, call=call, data=data, names=names, 
                predictArgs=predictArgs, envir=envir)
        )
    }
    # define a function the evaluates the expression to obtain the out-of-bag 
    # predictions and computes the prediction error for one bootstrap sample 
    bootFun <- function(r, keepSE) {
        s <- getIndices(splits, r)
        yHat <- eval(bootExpr)    # obtain out-of-bag predictions
        y <- dataSubset(y, -s)    # out-of-bag response values
        # compute cost function for obtain out-of-bag predictions
        pe <- computeCost(cost, y, yHat, args=costArgs, keepSE=keepSE)
    }
    # perform bootstrap replications
    R <- splits$R
    if(R == 1 && splits$type != "0.632") {
        pe <- bootFun(R, keepSE=TRUE)
        if(!is.list(pe)) pe <- t(pe)
    } else {
        pe <- lapply(seq_len(R), bootFun, keepSE=FALSE)  # for each replication
        pe <- combineData(pe, drop=FALSE)                # combine results
    }
    ## if requested, compute the 0.632 estimator
    if(splits$type == "0.632") {
        # check if the apparent error has already been computed
        if(is.null(ae <- splits$ae)) {
            # fit the model from all observations
            if(is.null(data)) {
                call[[names[1]]] <- x
                call[[names[2]]] <- y
            } else {
                call[[names]] <- data
            }
            fit <- eval(call, envir)
            # predict the response for all observations
            if(is.null(data)) {
                yHat <- doCall(predict, fit, x, args=predictArgs)
            } else {
                yHat <- doCall(predict, fit, data, args=predictArgs)
            }
            # compute the apparent error with the supplied cost function
            ae <- computeCost(cost, y, yHat, args=costArgs, keepSE=FALSE)
        }
        # compute 0.632 estimator as combination of apparent error and 
        # out-of-bag prediction error
        pe <- sweep(0.632 * pe, 2, 0.368 * ae, "+", check.margin=FALSE)
    }
    # return prediction error (add default names if necessary)
    addNames(pe)
}

# fit the model for the bootstrap sample and obtain predictions for the 
# out-of-bag observations for functions that take the predictors and the 
# response as separate arguments
bootXY <- function(i, call, x, y, names, predictArgs, envir) {
    # plug training data into function call
    call[[names[1]]] <- dataSubset(x, i)
    call[[names[2]]] <- dataSubset(y, i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predict, fit, dataSubset(x, -i), args=predictArgs)
}

# fit the model for the bootstrap sample and obtain predictions for the 
# out-of-bag observations for functions that have one argument for all the data
bootData <- function(i, call, data, names, predictArgs, envir) {
    # plug training data into function call
    call[[names]] <- dataSubset(data, i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predict, fit, dataSubset(data, -i), args=predictArgs)
}
