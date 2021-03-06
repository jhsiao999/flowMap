% Generated by roxygen2 (4.0.1): do not edit by hand
\name{statCrossLists}
\alias{statCrossLists}
\title{summary statistic across elements of lists}
\usage{
statCrossLists(obj, STAT)
}
\arguments{
\item{obj}{an object of one or more lists. List elements must be matrices of sample dimensions}

\item{STAT}{the statistic to be computed for each entry in the matrix across lists. Default value: median.}
}
\value{
est a matrix of the estimated statistic across random draws.
}
\description{
Compute an estimate of F-R statistics across random draws of cell population comparison.
}

