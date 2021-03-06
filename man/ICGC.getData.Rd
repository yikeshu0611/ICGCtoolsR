% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICGC.getData.R
\name{ICGC.getData}
\alias{ICGC.getData}
\title{To download and import ICGC data}
\usage{
ICGC.getData(Release, Project, donor = FALSE, specimen = FALSE,
  exp_array = FALSE, exp_seq = FALSE)
}
\arguments{
\item{Release}{the release verson, EX: release_28 or current}

\item{Project}{the project of data, EX: ALL-US}

\item{donor}{a logical param. TRUE: to download this data. FALSE, as defaulted}

\item{specimen}{a logical param. TRUE: to download this data. FALSE, as defaulted}

\item{exp_array}{a logical param. TRUE: to download this data. FALSE, as defaulted}

\item{exp_seq}{a logical param. TRUE: to download this data. FALSE, as defaulted}
}
\value{
a list
}
\description{
This function can be used to download and import ICGC data. Release and Project params can be finded in ICGC web. The other params are all logical words, which False as defaulted.
}
\examples{
ICGC.getData(Release="current",Project="ALL-US",donor=TRUE)
}
\author{
Jing Zhang
}
