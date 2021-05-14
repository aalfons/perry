# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' @export
setupPerryPlot <- function(object, ...) UseMethod("setupPerryPlot")

#' @export
setupPerryPlot.perry <- function(object, select = NULL, reps = NULL,
                                 seFactor = 1, ...) {
  # initializations
  if (is.null(reps)) reps <- object$splits$R > 1
  else reps <- isTRUE(reps)
  # extract subset of models
  object <- subset(object, select = select)
  if (reps) {
    PE <- object$reps
    if (is.null(PE)) stop("replications not available")
    else PE <- as.data.frame(PE)
  } else PE <- as.data.frame(t(object$pe))
  if (npe(object) == 0) stop("empty prediction error object")
  # stack selected results on top of each other
  fitName <- defaultFitNames(1)
  peName <- defaultNames(1)
  peNames <- peNames(object)
  n <- nrow(PE)
  Fit <- data.frame(Fit = rep.int(fitName, n))
  # no column for conditional plots if there is only one method with default
  # name
  if (isTRUE(peNames == peName)) {
    PE <- cbind(Fit, PE)
    facets <- NULL
  } else {
    PE <- lapply(peNames,
                 function(j) cbind(Fit, Name = rep.int(j, n), PE = PE[, j]))
    PE <- do.call(rbind, PE)
    names(PE) <- c("Fit", "Name", peName)
    facets <- ~ Name
  }
  # add data for error bars unless all replications are requested
  if (reps) includeSE <- NULL
  else {
    includeSE <- !is.null(seFactor) && !is.na(seFactor) &&
      !all(is.na(object$se))
    if (includeSE) {
      halflength <- seFactor * object$se
      PE$Lower <- PE[, peName] - halflength
      PE$Upper <- PE[, peName] + halflength
    }
  }
  # construct object to return
  out <- list(PE = PE, reps = reps)
  if (!is.null(includeSE)) out$includeSE <- includeSE
  if (!is.null(facets)) out$facets <- facets
  class(out) <- "setupPerryPlot"
  out
}

#' @export
setupPerryPlot.perrySelect <- function(object, subset = NULL, select = NULL,
                                       reps = NULL, seFactor = NULL, ...) {
  # initializations
  if (is.null(reps)) reps <- object$splits$R > 1
  else reps <- isTRUE(reps)
  if (is.null(seFactor)) seFactor <- object$seFactor
  # extract subset of models
  object <- subset(object, subset = subset, select = select)
  fits <- fits(object)
  if (reps) {
    PE <- object$reps
    if (is.null(PE)) stop("replications not available")
  } else PE <- object$pe
  if (nfits(object) == 0 || npe(object) == 0) {
    stop("empty prediction error object")
  }
  # ensure that models are shown in the correct order and drop unused levels
  # ensure that correct values are shown for a numeric tuning parameter
  if (!is.numeric(PE[, "Fit"])) {
    fits <- fits(object)
    PE$Fit <- factor(PE[, "Fit"], levels = fits)
  }
  # stack selected results on top of each other
  peName <- defaultNames(1)
  peNames <- peNames(object)
  n <- nrow(PE)
  # no column for conditional plots if there is only one column of results
  # with default name
  if (isTRUE(peNames == peName)) facets <- NULL
  else {
    Fit <- PE[, "Fit", drop = FALSE]
    PE <- lapply(peNames,
                 function(j) cbind(Fit, Name = rep.int(j, n), PE = PE[, j]))
    PE <- do.call(rbind, PE)
    names(PE) <- c("Fit", "Name", peName)
    facets <- ~ Name
  }
  # add data for error bars unless all replications are requested
  if (reps) includeSE <- NULL
  else {
    includeSE <- !is.null(seFactor) && !is.na(seFactor) &&
      !all(is.na(object$se[, peNames]))
    if (includeSE) {
      halflength <- seFactor * unlist(object$se[, peNames], use.names = FALSE)
      PE$Lower <- PE[, peName] - halflength
      PE$Upper <- PE[, peName] + halflength
    }
  }
  # construct object to return
  out <- list(PE = PE, reps = reps)
  if (!is.null(includeSE)) out$includeSE <- includeSE
  if (!is.null(facets)) out$facets <- facets
  class(out) <- "setupPerryPlot"
  out
}

#' @export
setupPerryPlot.perryTuning <- function(object, ...) {
  # adjust column specifying the model in case of only one tuning parameter
  if(ncol(object$tuning) == 1) fits(object) <- object$tuning[, 1]
  # call method for class "perrySelect"
  setupPerryPlot.perrySelect(object, ...)
}
