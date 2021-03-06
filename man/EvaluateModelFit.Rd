% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EvaluateModelFit.R
\name{EvaluateModelFit}
\alias{EvaluateModelFit}
\title{Evaluate model fit}
\usage{
EvaluateModelFit(l, k, n)
}
\arguments{
\item{l}{the value of the log-likelihood function evaluated at MLEs of the parameters.}

\item{k}{the number of parameters in the model to be estimated.}

\item{n}{the number of observations used to estimate the parameters.}
}
\value{
A numeric vector containing (in order) AIC, AICc, BIC, and maximum log-likelihood.
}
\description{
Computes the AIC, corrected AIC (AICc), and BIC for a model.
}
