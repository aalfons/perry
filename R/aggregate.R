# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @S3method aggregate perry
aggregate.perry <- function(x, FUN = mean, select = NULL, ...) {
    if(is.null(PE <- x$reps)) PE <- t(x$pe)  # matrix is required
    if(!is.null(select)) PE <- PE[, select, drop=FALSE]
    apply(PE, 2, FUN=FUN, ...)
}


#' @S3method aggregate perrySelect
aggregate.perrySelect <- function(x, FUN = mean, select = NULL, ...) {
    by <- "Fit"
    if(is.null(PE <- x$reps)) PE <- x$pe
    peNames <- peNames(x)
    if(is.null(select)) {
        select <- peNames
    } else if(!is.character(select)) select <- peNames[select]
    aggregate(PE[, select, drop=FALSE], by=PE[, by, drop=FALSE], FUN=FUN, ...)
}


#' @S3method aggregate perryTuning
aggregate.perryTuning <- function(x, ...) {
    # call method for class "perrySelect"
    out <- aggregate.perrySelect(x, ...)
    # replace column specifying the fit by grid of tuning parameters
    cbind(x$tuning, out[, -1, drop=FALSE])
}
