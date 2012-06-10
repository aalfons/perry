\name{randomSplits}
\alias{print.randomSplits}
\alias{randomSplits}
\title{Random data splits}
\usage{
  randomSplits(n, m, R = 1, grouping = NULL)
}
\arguments{
  \item{n}{an integer giving the number of observations to
  be split into training and test data.  This is ignored if
  \code{grouping} is supplied in order to split groups of
  observations into folds.}

  \item{m}{an integer giving the number of observations or
  groups of observations to be used as test data.}

  \item{R}{an integer giving the number of random data
  splits.}

  \item{grouping}{a factor specifying groups of
  observations.  If supplied, the data are split according
  to the groups rather than individual observations such
  that all observations within a group belong either to the
  training or test data.}
}
\value{
  An object of class \code{"randomSplits"} with the
  following components:

  \item{n}{an integer giving the number of observations or
  groups.}

  \item{m}{an integer giving the number of observations or
  groups in the test data.}

  \item{R}{an integer giving the number of random data
  splits.}

  \item{subsets}{an integer matrix in which each column
  contains the indices of the observations or groups in the
  test data of the corresponding random data split.}

  \item{grouping}{a list giving the indices of the
  observations belonging to each group.  This is only
  returned if a grouping factor has been supplied.}
}
\description{
  Split observations or groups of observations into
  training and test data to be used for (repeated) random
  splitting (also known as random subsampling or Monte
  Carlo cross-validation).
}
\note{
  This is a simple wrapper function for
  \code{\link{perrySplits}} with a control object generated
  by \code{\link{splitControl}}.
}
\examples{
set.seed(1234)  # set seed for reproducibility
randomSplits(20, m = 5)
randomSplits(20, m = 5, R = 10)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perrySplits}}, \code{\link{splitControl}},
  \code{\link{cvFolds}}
}
\keyword{utilities}

