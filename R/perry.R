# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Resampling-based prediction error for fitted models
#' 
#' Generic function to estimate the prediction error of a fitted model via 
#' (repeated) \eqn{K}-fold cross-validation, (repeated) random splitting (also 
#' known as random subsampling or Monte Carlo cross-validation), or the 
#' bootstrap.
#' 
#' The idea is that developers write easy-to-use methods for end users to 
#' leverage the prediction error estimation framework for their models.  A 
#' typical \code{perry} method consists of the following two parts: first the 
#' data are extracted from the model, then function \code{\link{perryFit}} is 
#' called to perform prediction error estimation.  The programming effort of 
#' implementing prediction error estimation for a certain model is thus greatly 
#' reduced.
#' 
#' Examples for methods are available in package perryExamples (see 
#' \code{\link[perryExamples]{perry-methods}}).
#' 
#' @param object  the fitted model for which to estimate the prediction error.
#' @param \dots  additional arguments to be passed down to methods.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perryFit}}, \code{\link[perryExamples]{perry-methods}}
#' 
#' @keywords utilities
#' 
#' @export

perry <- function(object, ...) UseMethod("perry")
