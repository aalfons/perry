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
    PE <- as.data.frame(object$reps)
    if (is.null(PE)) stop("replications not available")
  } else PE <- as.data.frame(t(object$pe))
  if (npe(object) == 0) stop("empty prediction error object")
  # stack selected results on top of each other
  fitName <- defaultFitNames(1)
  peName <- defaultNames(1)
  peNames <- peNames(object)
  n <- nrow(PE)
  Fit <- data.frame(Fit = rep.int(fitName, n))
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
    if (is.null(seFactor)) seFactor <- NA
    includeSE <- !is.na(seFactor) && !all(is.na(object$se))
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
