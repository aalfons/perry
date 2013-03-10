\name{accessors}
\alias{accessors}
\alias{fits}
\alias{fits<-}
\alias{nfits}
\alias{npe}
\alias{peNames}
\alias{peNames<-}
\title{Access or set information on resampling-based prediction error results}
\usage{
  peNames(x)

  peNames(x) <- value

  fits(x)

  fits(x) <- value

  npe(x)

  nfits(x)
}
\arguments{
  \item{x}{an object inheriting from class \code{"perry"}
  or \code{"perrySelect"} that contains prediction error
  results.}

  \item{value}{a vector of replacement values.}
}
\value{
  \code{peNames} returns the names of the prediction error
  results.  The replacement function thereby returns them
  invisibly.

  \code{fits} returns the identifiers of the models for
  objects inheriting from class \code{"perrySelect"} and
  \code{NULL} for objects inheriting from class
  \code{"perry"}.  The replacement function thereby returns
  those values invisibly.

  \code{npe} returns the number of prediction error
  results.

  \code{nfits} returns the number of models included in
  objects inheriting from class \code{"perrySelect"} and
  \code{NULL} for objects inheriting from class
  \code{"perry"}.
}
\description{
  Retrieve or set the names of resampling-based prediction
  error results, retrieve or set the identifiers of the
  models, or retrieve the number of prediction error
  results or included models.
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
cv50 <- perry(fit50, splits = folds, fit = "both", 
    cost = rtmspe, trim = 0.1)

# 75\% subsets
fit75 <- ltsReg(Y ~ ., data = coleman, alpha = 0.75)
cv75 <- perry(fit75, splits = folds, fit = "both", 
    cost = rtmspe, trim = 0.1)

# combine results into one object
cv <- perrySelect("0.5" = cv50, "0.75" = cv75)
cv

# "perry" object
npe(cv50)
nfits(cv50)
peNames(cv50)
peNames(cv50) <- c("improved", "initial")
fits(cv50)
cv50

# "perrySelect" object
npe(cv)
nfits(cv)
peNames(cv)
peNames(cv) <- c("improved", "initial")
fits(cv)
fits(cv) <- 1:2
cv
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}}, \code{\link{perrySelect}},
  \code{\link{perryTuning}}
}
\keyword{utilities}

