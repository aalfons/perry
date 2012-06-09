# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @export
peNames <- function(x) UseMethod("peNames")

#' @S3method peNames perry
peNames.perry <- function(x) names(x$pe)

#' @S3method peNames perrySelect
peNames.perrySelect <- function(x) names(x$pe)[-1]


#' @export
"peNames<-" <- function(x, value) UseMethod("peNames<-")

#' @S3method peNames<- perry
"peNames<-.perry" <- function(x, value) {
    object <- x
    names(object$pe) <- names(object$se) <- value
    if(!is.null(x$reps)) colnames(object$reps) <- value
    eval.parent(substitute(x <- object))
}

#' @S3method peNames<- perrySelect
"peNames<-.perrySelect" <- function(x, value) {
    object <- x
    names(object$best) <- value
    value <- c("Fit", value)
    names(object$pe) <- names(object$se) <- value
    if(!is.null(x$reps)) names(object$reps) <- value
    eval.parent(substitute(x <- object))
}


#' @export
fits <- function(x) UseMethod("fits")

#' @S3method fits perry
fits.perry <- function(x) NULL

#' @S3method fits perrySelect
fits.perrySelect <- function(x) x$pe$Fit


#' @export
"fits<-" <- function(x, value) UseMethod("fits<-")

#' @S3method fits<- perry
"fits<-.perry" <- function(x, value) eval.parent(substitute(x))

#' @S3method fits<- perrySelect
"fits<-.perrySelect" <- function(x, value) {
    object <- x
    if(is.factor(value)) value <- factor(as.character(value), levels=value)
    object$pe$Fit <- object$se$Fit <- value
    if(!is.null(reps <- x$reps)) {
        indices <- match(reps$Fit, x$pe$Fit, nomatch=0)
        object$reps$Fit <- value[indices]
    }
    eval.parent(substitute(x <- object))
}


#' @export
npe <- function(x) UseMethod("npe")

#' @S3method npe perry
npe.perry <- function(x) length(x$pe)

#' @S3method npe perrySelect
npe.perrySelect <- function(x) ncol(x$pe) - 1


#' @export
nfits <- function(x) UseMethod("nfits")

#' @S3method nfits perry
#' @S3method nfits perrySelect
nfits.perry <- nfits.perrySelect <- function(x) nrow(x$pe)
