# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Reshape resampling-based prediction error results
#' 
#' Reshape resampling-based prediction error results into an object of class 
#' \code{"perrySelect"} with only one column of results.
#' 
#' @param x  an object inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contains prediction error results.
#' @param selectBest  a character string specifying a criterion for selecting 
#' the best model.  Possible values are \code{"min"} (the default) or 
#' \code{"hastie"}.  The former selects the model with the smallest prediction 
#' error.  The latter is useful for nested models or for models with a tuning 
#' parameter controlling the complexity of the model (e.g., penalized 
#' regression).  It selects the most parsimonious model whose prediction error 
#' is no larger than \code{seFactor} standard errors above the prediction error 
#' of the best overall model.  Note that the models are thereby assumed to be 
#' ordered from the most parsimonious one to the most complex one.  In 
#' particular a one-standard-error rule is frequently applied.
#' @param seFactor  a numeric value giving a multiplication factor of the 
#' standard error for the selection of the best model.  This is ignored if 
#' \code{selectBest} is \code{"min"}.
#' @param tuning  a list of tuning parameter values that correspond to the 
#' different prediction error results.  The names of the list components should 
#' thereby correspond to the argument names of the tuning parameters.  For each 
#' tuning parameter, a vector of values can be supplied.  A data frame 
#' containing all possible combinations of tuning parameter values is then 
#' added to the reshaped prediction error results.
#' @param \dots  additional arguments to be passed down.
#' 
#' @return  An object of class \code{"perrySelect"} (subclass 
#' \code{"perryTuning"} if a list of tuning parameters is supplied) with the 
#' following components:
#' @returnItem pe  a data frame containing the estimated prediction errors for 
#' the models.  In case of more than one resampling replication, those are 
#' average values over all replications.
#' @returnItem se  a data frame containing the estimated standard errors of the 
#' prediction loss for the models.
#' @returnItem reps  a data frame containing the estimated prediction errors 
#' for the models from all replications.  This is only returned in case of more 
#' than one resampling replication.
#' @returnItem splits  an object giving the data splits used to estimate the 
#' prediction error.
#' @returnItem y  the response.
#' @returnItem yHat  a list containing the predicted values for the 
#' models.  Each list component is again a list containing the corresponding 
#' predicted values from all replications.
#' @returnItem best  an integer giving the index of the model with the best 
#' prediction performance.
#' @returnItem selectBest  a character string specifying the criterion used for 
#' selecting the best model.
#' @returnItem seFactor  a numeric value giving the multiplication factor of 
#' the standard error used for the selection of the best model.
#' @returnItem tuning  a data frame containing the grid of tuning parameter 
#' values that correspond to the different prediction error results (only 
#' subclass \code{"perryTuning"}).
#' 
#' @author Andreas Alfons
#' 
#' @references 
#' Hastie, T., Tibshirani, R. and Friedman, J. (2009) \emph{The Elements of 
#' Statistical Learning: Data Mining, Inference, and Prediction}.  Springer, 
#' 2nd edition.
#' 
#' @seealso \code{\link{perryFit}}, \code{\link{perrySelect}}, 
#' \code{\link{perryTuning}}
#' 
#' @example inst/doc/examples/example-perryReshape.R
#' 
#' @keywords utilities
#' 
#' @export

perryReshape <- function(x, selectBest = c("min", "hastie"), seFactor = 1, 
                         tuning = list(), ...) {
  # initializations
  if(!inherits(x, c("perry", "perrySelect"))) 
    stop("object must inherit from class \"perry\" or \"perrySelect\"")
  if(npe(x) == 0 || isTRUE(nfits(x) == 0)) stop("empty object")
  peNames <- peNames(x)
  # create list of objects with one column
  peName <- defaultNames(1)
  objects <- lapply(peNames, function(s) {
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
  pe <- perrySelect(.list=objects, .selectBest=selectBest, .seFactor=seFactor)
  # add tuning parameters if requested
  if(length(tuning) > 0) {
    n <- nfits(pe)
    tuning <- do.call(expand.grid, tuning)
    if(nrow(tuning) != n) {
      stop(sprintf("'tuning' must correspond to %d models", n))
    }
    pe$tuning <- tuning
    class(pe) <- c("perryTuning", "perrySelect")
  }
  # return object
  pe
}
