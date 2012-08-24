\name{perry-package}
\alias{perry-package}
\docType{package}
\title{
Resampling-based prediction error estimation for regression
models
}
\description{
Tools that allow developers to write functions for prediction
    error estimation with minimal programming effort and assist users with
    model selection in regression problems.
}
\details{
\tabular{ll}{
Package: \tab perry\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2012-08-23\cr
Depends: \tab R (>= 2.14.1), ggplot2, parallel, robustbase\cr
Imports: \tab stats\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

Index:
\preformatted{
accessors               Access or set information on resampling-based
                        prediction error results
aggregate.perry         Aggregate resampling-based prediction error
                        results
bootControl             Control object for bootstrap samples
bootPE                  Bootstrap prediction error estimation for
                        fitted models
bootSamples             Bootstrap samples
cost                    Prediction loss
cvFolds                 Cross-validation folds
foldControl             Control object for cross-validation folds
fortify.perry           Convert resampling-based prediction error
                        results into a data frame for plotting
perry                   Resampling-based prediction error for fitted
                        models
perry-package           Resampling-based prediction error estimation
                        for regression models
perryFit                Resampling-based prediction error for model
                        evaluation
perryPlot               Plot resampling-based prediction error results
perryReshape            Reshape resampling-based prediction error
                        results
perrySelect             Model selection via resampling-based prediction
                        error measures
perrySplits             Data splits for resampling-based prediction
                        error measures
perryTuning             Resampling-based prediction error for tuning
                        parameter selection
randomSplits            Random data splits
repCV                   (Repeated) cross-validation for fitted models
repRS                   (Repeated) random splitting for fitted models
reperry                 Recompute resampling-based prediction error
                        measures
splitControl            Control object for random data splits
subset.perry            Subsetting resampling-based prediction error
                        results
summary.perry           Summarize resampling-based prediction error
                        results
}
}
\author{
Andreas Alfons [aut, cre]

Maintainer: Andreas Alfons <andreas.alfons@kuleuven.be>
}
\keyword{package}
