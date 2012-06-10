# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## workhorse function
plotPE <- function(x, 
        method = c("bwplot", "densityplot", "xyplot", "dotplot"), ...) {
    ## initializations
    if(x$splits$R == 1) {
        choices <- eval(formals(sys.function())[["method"]])
        if(identical(method, choices)) {
            method <- "xyplot"
        } else method <- match.arg(method, c("xyplot", "dotplot"))
    } else method <- match.arg(method)
    ## call plot function
    if(method == "bwplot") {
        bwplot(x, ...)
    } else if(method == "densityplot") {
        densityplot(x, ...)
    } else if(method == "xyplot") {
        xyplot(x, ...)
    } else dotplot(x, ...)
}


#' @S3method plot perry
plot.perry <- plotPE

#' @S3method plot perrySelect
plot.perrySelect <- plotPE
