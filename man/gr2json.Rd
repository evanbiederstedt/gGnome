% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gGnome.R
\name{gr2json}
\alias{gr2json}
\title{gr2json}
\usage{
gr2json(intervals, filename, y = rep("null", length(intervals)),
  labels = "", maxcn = 100, maxweight = 100)
}
\arguments{
\item{filename}{output json file}

\item{GRange}{input jab object}
}
\description{
Dumps GRanges into JSON with metadata features as data points in  "intervals"
}
\author{
Marcin Imielinski
}
