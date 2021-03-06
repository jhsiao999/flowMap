% Generated by roxygen2 (4.0.1): do not edit by hand
\name{getFRmat}
\alias{getFRmat}
\title{FR tests to compare two flow cytometry samples}
\usage{
getFRmat(XX1, XX2, sampleMethod, sampleSize, i = NULL)
}
\arguments{
\item{XX1}{a flow cytometry sample of cell populations, organized in a matrix or a data.frame of events (rows)
by features (columns) where cell populaiton memberships are indexed in the last column by
a variable named id.}

\item{XX2}{a flow cytometry sample of cell populations, organized in a matrix or a data.frame of events (rows)
by features (columns) where cell populaiton memberships are indexed in the last column by
a variable named id.}

\item{sampleMethod}{downsampling method, options include \emph{equalSize} or \emph{proportional}. Both methods
sample events without replacement from the combined events in a single cell population pair
comparison. Using \emph{equalSize}, each sample includes an equal number of events from the
two cell populations being compared. Using \emph{proportional}, the ratio of the event
membership is same as the ratio of event membershiop prior to sampling.
(default: \emph{proportional})}

\item{sampleSize}{specifies \emph{S},the number of events to be included in each sample. For \emph{equalSize} sampling
, S/2 is sampled from each population. For \emph{proportional} sampling, the ratio of event membership
is the same as the ratio of event membership prior to sampling.}

\item{i}{dummy variable to initialize parallel computing (\pkg{doParallel}).}
}
\value{
wmat a list of matrices containing sample FR statistics for each XX1 by XX2 population comparisons.

runsmat a list of matrices containing sample runs for each XX1 by XX2 population comparisons.

mumat a list of matrices containing sample expected number of runs for each XX1 by XX2 population comparisons.

sigma2mat a list of matrices containing sample estimated variance of runs for each XX1 by XX2 population comparisons.

pNormat a list of matrices containing sample one-sided p-value associated with the FR statistics of
               XX1 by XX2 population comparisons under the asymptotic normality assumption of the FR statistic.
}
\description{
Perform FR tests to compare cell populations across two flow cytometry samples.
The FR statistics are estimated based on a downsampling scheme designed to optimize runtime as well as
precision and accuracy of the statistics. The downsampling scheme samples from the pooled data of
events across a single cell population comparison. The events in the sample maintain the same cell
population membership ratio as the events in the pooled data (\emph{proportional}) or include an
equal number of events from the two cell populations in the comparison.
}
\author{
Chiaowen Joyce Hsiao \email{joyce.hsiao1@gmail.com}
}

