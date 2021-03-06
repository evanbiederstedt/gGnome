% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gGnome.R
\name{e2j}
\alias{e2j}
\title{given node and edges set, return GRangesList of the junction locations}
\usage{
e2j(segs, es, etype = "aberrant")
}
\arguments{
\item{segs}{the GRanges of nodes}

\item{es}{the data.table of edges}

\item{etype}{\code{character} scalar, specifying the type of edges to convert}
}
\value{
\code{junctions} object specifying the breakpoint pair and orientations of selected edges
}
\description{
Any adjacency between DNA segments can be defined as a pair of locations on the reference genome,
each with the orientation specifying the fused side of the molecule. This data is easily represented
by a \code{GRangesList}.
}
\details{
\code{etype} can be any of the following values, "all", "aberrant", "reference", "loose". It can also be
regular expressions that match these values.
}
