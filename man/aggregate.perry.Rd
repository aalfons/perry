\name{aggregate.perry}
\alias{aggregate.perry}
\alias{aggregate.perrySelect}
\alias{aggregate.perryTuning}
\title{Aggregate resampling-based prediction error results}
\usage{
  \method{aggregate}{perry} (x, FUN = mean, select = NULL,
    ...)

  \method{aggregate}{perrySelect} (x, FUN = mean,
    select = NULL, ...)

  \method{aggregate}{perryTuning} (x, ...)
}
\arguments{
  \item{x}{an object inheriting from class \code{"perry"}
  or \code{"perrySelect"} that contains prediction error
  results (note that the latter includes objects of class
  \code{"perryTuning"}).}

  \item{FUN}{a function to compute the summary statistics.}

  \item{select}{a character, integer or logical vector
  indicating the columns of prediction error results for
  which to compute the summary statistics.}

  \item{\dots}{for the \code{"perryTuning"} method,
  additional arguments to be passed to the
  \code{"perrySelect"} method.  Otherwise additional
  arguments to be passed to \code{FUN}.}
}
\value{
  The \code{"perry"} method returns a vector or matrix of
  aggregated prediction error results, depending on whether
  \code{FUN} returns a single value or a vector.

  For the other methods, a data frame containing the
  aggregated prediction error results for each model is
  returned.  In the case of the \code{"perryTuning"}
  method, the data frame contains the combinations of
  tuning parameters rather than a column describing the
  models.
}
\description{
  Compute summary statistics of resampling-based prediction
  error results.
}
\note{
  Duplicate indices in \code{subset} or \code{select} are
  removed such that all models and prediction error results
  are unique.
}
\examples{
data("coleman")
set.seed(1234)  # set seed for reproducibility

## set up folds for cross-validation
folds <- cvFolds(nrow(coleman), K = 5, R = 10)

## compare raw and reweighted LTS estimators for 
## 50\% and 75\% subsets

# 50\% subsets
fit50 <- ltsReg(Y ~ ., data = coleman, alpha = 0.5)
cv50 <- repCV(fit50, folds = folds, fit = "both", 
    cost = rtmspe, trim = 0.1)

# 75\% subsets
fit75 <- ltsReg(Y ~ ., data = coleman, alpha = 0.75)
cv75 <- repCV(fit75, folds = folds, fit = "both", 
    cost = rtmspe, trim = 0.1)

# combine results into one object
cv <- perrySelect("0.5" = cv50, "0.75" = cv75)
cv

# summary of the results with the 50\% subsets
aggregate(cv50, summary)
# summary of the combined results
aggregate(cv, summary)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}}, \code{\link{perrySelect}},
  \code{\link{perryTuning}}, \code{\link[stats]{aggregate}}
}
\keyword{utilities}

