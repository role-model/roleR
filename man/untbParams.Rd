% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roleParams.R
\name{untbParams}
\alias{untbParams}
\title{Wrapper around roleParams to create a "untb-flavored" (Unified Neutral Theory of Biodiversity) roleModel}
\usage{
untbParams(
  individuals_local,
  individuals_meta,
  species_meta,
  speciation,
  dispersal_prob,
  init_type,
  niter,
  niterTimestep
)
}
\arguments{
\item{individuals_local}{individuals_local}

\item{individuals_meta}{individuals_meta}

\item{species_meta}{species_meta}

\item{speciation}{speciation}

\item{dispersal_prob}{dispersal_prob}

\item{init_type}{init_type}

\item{niter}{niter}

\item{niterTimestep}{niterTimestep}
}
\value{
a `roleParams` object
}
\description{
Only arguments relevant to a UNTB neutral model is included
}
