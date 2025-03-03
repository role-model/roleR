% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roleExperiment.R
\docType{methods}
\name{[,roleExperiment,ANY,ANY,ANY-method}
\alias{[,roleExperiment,ANY,ANY,ANY-method}
\alias{[,roleExperiment-method}
\alias{$}
\alias{$,roleExperiment-method}
\title{Extract parts of `roleExperiment` object}
\usage{
\S4method{[}{roleExperiment,ANY,ANY,ANY}(x, i, j, ..., drop = FALSE)

\S4method{$}{roleExperiment}(x, name)
}
\arguments{
\item{x}{a `roleExperiment` object}

\item{i}{row index}

\item{j}{column index}

\item{...}{additional arguemnts}

\item{drop}{boolean, whether to drop dimensions}

\item{name}{column name to extract}
}
\description{
Extract parts of `roleExperiment` object

Extract columns from `info` data.frame
}
