# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @export
perryTuning <- function(object, ...) UseMethod("perryTuning")

#' @S3method perryTuning function
perryTuning.function <- function(object, formula, data = NULL, x = NULL, y, 
        tuning = list(), args = list(), cost = rmspe, splits = foldControl(), 
        names = NULL, predictArgs = list(), costArgs = list(), 
        selectBest = c("min", "hastie"), seFactor = 1, 
        envir = parent.frame(), seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perryTuning")
    call <- as.call(c(object, args))  # set up unevaluated function call
    haveFormula <- !missing(formula)
    if(haveFormula || !missing(data)) {
        if(is.null(names)) names <- c("formula", "data")
        if(haveFormula) call[[names[1]]] <- formula
        names <- names[-1]
        mf <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data"), names(mf), 0)
        mf <- mf[c(1, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1]] <- as.name("model.frame")
        data <- eval(mf, envir)
        if(is.empty.model(attr(data, "terms"))) stop("empty model")
        y <- model.response(data)  # extract response from model frame
    }
    ## call method for unevaluated function calls
    out <- perryTuning(call, data=data, x=x, y=y, tuning=tuning, cost=cost, 
        splits=splits, names=names, predictArgs=predictArgs, costArgs=costArgs, 
        selectBest=selectBest, seFactor=seFactor, envir=envir, seed=seed)
    out$call <- matchedCall
    out
}

#' @S3method perryTuning call
perryTuning.call <- function(object, data = NULL, x = NULL, y, 
        tuning = list(), cost = rmspe, splits = foldControl(), 
        names = NULL, predictArgs = list(), costArgs = list(), 
        selectBest = c("min", "hastie"), seFactor = 1, 
        envir = parent.frame(), seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perryTuning")
    n <- nobs(y)
    if(is.null(data)) {
        sx <- "x"
        nx <- nobs(x)
    } else {
        sx <- "data"
        nx <- nobs(data)
    }
    if(!isTRUE(n == nx)) stop(sprintf("'%s' must have %d observations", sx, nx))
    selectBest <- match.arg(selectBest)
    # create all combinations of tuning parameters
    tuning <- do.call("expand.grid", tuning)
    nTuning <- nrow(tuning)
    pTuning <- ncol(tuning)
    if(nTuning == 0 || pTuning == 0) {
        # use function perryFit() if no tuning parameters are supplied
        out <- perryFit(object, data, x, y, cost=cost, splits=splits, 
            names=names, predictArgs=predictArgs, costArgs=costArgs, 
            envir=envir, seed=seed)
        return(out)
    }
    # make sure that .Random.seed exists if no seed is supplied
    if(is.null(seed)) {
        if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) runif(1)
        seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
    } else set.seed(seed)
    ## compute data splits
    if(is.null(getS3method("perryTool", class(splits), optional=TRUE))) {
        splits <- perrySplits(n, control=splits)
    }
    ## estimate the prediction error for each combination of tuning parameters
    tuningNames <- names(tuning)
    pe <- lapply(seq_len(nTuning), 
        function(i) {
            # add tuning parameters to function call
            for(j in seq_len(pTuning)) {
                object[[tuningNames[j]]] <- tuning[i, j]
            }
            # estimate the prediction error
            perryTool(object, data, x, y, cost=cost, splits=splits, 
                names=names, predictArgs=predictArgs, costArgs=costArgs, 
                envir=envir)
        })
    ## prepare the prediction errors and standard errors
    if(is.list(pe[[1]])) {
        R <- 1  # only one replication
        se <- combineData(lapply(pe, "[[", 2), drop=FALSE)
        pe <- combineData(lapply(pe, "[[", 1), drop=FALSE)
    } else {
        R <- nrow(pe[[1]])
        pe <- combineData(pe, drop=FALSE)
        if(R <= 1) se <- matrix(NA, nrow(pe), ncol(pe), dimnames=dimnames(pe))
    }
    pe <- data.frame(Fit=rep(seq_len(nTuning), each=R), pe, row.names=NULL)
    if(R > 1) {
        # compute average prediction errors and standard errors
        reps <- pe
        pe <- aggregate(reps[, -1, drop=FALSE], reps[, 1, drop=FALSE], mean)
        se <- aggregate(reps[, -1, drop=FALSE], reps[, 1, drop=FALSE], sd)
    } else se <- data.frame(Fit=seq_len(nTuning), se, row.names=NULL)
    ## find optimal values of tuning parameters
    if(selectBest == "min") {
        seFactor <- NA
        best <- sapply(pe[, -1, drop=FALSE], selectMin)
    } else {
        seFactor <- rep(seFactor, length.out=1)
        best <- sapply(names(pe)[-1], 
            function(j) selectHastie(pe[, j], se[, j], seFactor=seFactor))
    }
    ## construct return object
    out <- list(splits=splits, tuning=tuning, best=best, pe=pe, se=se, 
        selectBest=selectBest, seFactor=seFactor)
    if(R > 1) out$reps <- reps
    out$seed <- seed
    out$call <- matchedCall
    class(out) <- c("perryTuning", "perrySelect")
    out
}
