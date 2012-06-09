# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @export
perryReshape <- function(x, ...) UseMethod("perryReshape")

#' @S3method perryReshape perry
perryReshape.perry <- function(x, 
    selectBest = c("min", "hastie"), seFactor = 1, ...) {
    # initializations
    if(npe(x) == 0 || isTRUE(nfits(x) == 0)) stop("empty object")
    peNames <- peNames(x)
    # create list of objects with one column
    peName <- defaultNames(1)
    objects <- lapply(peNames, 
        function(s) {
            xs <- subset(x, select=s)
            peNames(xs) <- peName
            xs
        })
    # substitute "PE" in default names by "Fit"
    if(identical(peNames, defaultNames(length(peNames)))) {
        fitName <- defaultFitNames(1)
        peNames <- gsub(peName, fitName, peNames, fixed=TRUE)
    }
    # call perrySelect() to combine the model fits
    names(objects) <- peNames
    objects$.selectBest <- selectBest
    objects$.seFactor <- seFactor
    do.call(perrySelect, objects)
}

#' @S3method perryReshape perrySelect
perryReshape.perrySelect <- perryReshape.perry
