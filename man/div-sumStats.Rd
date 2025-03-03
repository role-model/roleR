% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roleData-divStats.R
\docType{methods}
\name{hillAbund}
\alias{hillAbund}
\alias{hillAbund,roleData-method}
\alias{hillGenetic}
\alias{hillGenetic,roleData-method}
\alias{hillTrait}
\alias{hillTrait,roleData-method}
\alias{hillPhylo}
\alias{hillPhylo,roleData-method}
\alias{richness}
\alias{richness,roleData-method}
\title{Hill number summary statistics for `roleData` objects.}
\usage{
hillAbund(x, q = 1:4, ...)

\S4method{hillAbund}{roleData}(x, q = 1:4)

hillGenetic(x, q = 1:4, ...)

\S4method{hillGenetic}{roleData}(x, q = 1:4)

hillTrait(x, q = 1:4, ...)

\S4method{hillTrait}{roleData}(x, q = 1:4)

hillPhylo(x, q = 1:4, ...)

\S4method{hillPhylo}{roleData}(x, q = 1:4)

richness(x, ...)

\S4method{richness}{roleData}(x)
}
\arguments{
\item{x}{The `roleData` object to calculate Hill numbers or richness from}

\item{q}{The Hill number exponents. Can be a single value or a vector of integers.}

\item{...}{additional parameters, ignored}
}
\description{
These functions calculate Hill numbers. The naming convention 
    indicates what type of Hill number each function calculates; e.g. 
    `hillAbund` calculates abundance-based Hill numbers, while `hillTrait` 
    calculates trait-based Hill numbers.
}
\details{
See Gaggiotti et al. 2018 for an overview of Hill numbers in the context of ecology and evolution.
}
\references{
Gaggiotti, Oscar E, Anne Chao, Pedro Peres-Neto, Chun-Huo Chiu, Christine Edwards, Marie-Josée Fortin, Lou Jost, Christopher M Richards, and Kimberly A Selkoe. “Diversity from Genes to Ecosystems: A Unifying Framework to Study Variation across Biological Metrics and Scales.” Evol. Appl. 11, no. 7 (2018): 1176–93. https://doi.org/10.1111/eva.12593.
}
