# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @export
perrySelect <- function(..., .reshape = FALSE, 
    .selectBest = c("min", "hastie"), 
    .seFactor = 1) {
    ## initializations
    objects <- list(...)
    m <- length(objects)
    if(m == 0) stop("empty list of objects")
    # check class of objects
    isPerrySelect <- sapply(objects, inherits, "perrySelect")
    if(!all(sapply(objects, inherits, "perry") | isPerrySelect)) {
        stop("all objects must inherit from class \"perry\" or \"perrySelect\"")
    }
    .selectBest <- match.arg(.selectBest)
    # remove empty objects
    keep <- sapply(objects, function(x) npe(x) > 0 && !isTRUE(nfits(x) == 0))
    objects <- objects[keep]
    m <- length(objects)
    if(m == 0) stop("all objects are empty")
    isPerrySelect <- isPerrySelect[keep]
    # check if the same data splits have been used
    splits <- unique(lapply(objects, "[[", "splits"))
    if(length(splits) > 1) {
        stop("all object must be computed from the same data splits")
    } else splits <- splits[[1]]
    # check names for the supplied objects
    fits <- names(objects)
    if(is.null(fits)) {
        fits <- defaultFitNames(m)
    } else if(any(i <- fits == "")) fits[i] <- defaultFitNames(m)[i]
    names(objects) <- fits
    # check dimensions or reshape objects with more than one column
    d <- sapply(objects, npe)
    if(isTRUE(.reshape)) {
        .reshape <- which(d > 1)
        objects[.reshape] <- lapply(objects[.reshape], perryReshape)
        isPerrySelect[.reshape] <- TRUE
        d <- 1
    } else {
        d <- unique(d)
        if(length(d) > 1) stop("all objects must have the same dimension")
    }
    ## prepare objects of class "perrySelect"
    if(any(isPerrySelect)) {
        # prepare names
        fits <- as.list(fits)
        fits[isPerrySelect] <- mapply(function(f, x) paste(f, x$pe$Fit, sep="."), 
            fits[isPerrySelect], objects[isPerrySelect], SIMPLIFY=FALSE)
        fits <- unlist(fits)
        # prepare basic information
        objects[isPerrySelect] <- lapply(objects[isPerrySelect], 
            function(x) {
                # remove column specifying fit from results
                x$pe <- x$pe[, -1, drop=FALSE]
                x$se <- x$se[, -1, drop=FALSE]
                if(!is.null(x$reps)) x$reps <- x$reps[, -1, drop=FALSE]
                x
            })
    }
    ## combine prediction error results and standard errors
    pe <- lapply(objects, 
        function(x) {
            pe <- x$pe                                     # extract prediction errors
            if(is.null(dim(pe))) t(pe) else as.matrix(pe)  # return matrix
        })
    se <- lapply(objects, 
        function(x) {
            se <- x$se                                     # extract standard errors
            if(is.null(dim(se))) t(se) else as.matrix(se)  # return matrix
        })
    if(m > 1) {
        # check if names are the same for all objects
        peNames <- unique(lapply(pe, colnames))
        adjustNames <- length(peNames) > 1
        peNames <- if(adjustNames) defaultNames(d) else peNames[[1]]
    }
    pe <- combineData(pe, drop=FALSE)
    se <- combineData(se, drop=FALSE)
    if(m > 1 && adjustNames) colnames(pe) <- colnames(se) <- peNames
    pe <- data.frame(Fit=factor(fits, levels=fits), pe, row.names=NULL)
    se <- data.frame(Fit=factor(fits, levels=fits), se, row.names=NULL)
    ## combine repeated prediction error results
    # TODO: don't use information from 'splits'
    R <- splits$R
    haveReps <- R > 1
    if(haveReps) {
        reps <- lapply(objects, function(x) as.matrix(x$reps))
        reps <- do.call("rbind", reps)
        if(m > 1 && adjustNames) colnames(reps) <- peNames
        reps <- data.frame(Fit=factor(rep(fits, each=R), levels=fits), 
            reps, row.names=NULL)
    }
    ## find best model
    if(.selectBest == "min") {
        .seFactor <- NA
        best <- sapply(pe[, -1, drop=FALSE], selectMin)
    } else {
        .seFactor <- rep(.seFactor, length.out=1)
        best <- sapply(names(pe)[-1], 
            function(j) selectHastie(pe[, j], se[, j], seFactor=.seFactor))
    }
    ## construct return object
    out <- list(splits=splits, best=best, pe=pe, se=se, 
        selectBest=.selectBest, seFactor=.seFactor)
    if(haveReps) out$reps <- reps
    class(out) <- "perrySelect"
    out
}
