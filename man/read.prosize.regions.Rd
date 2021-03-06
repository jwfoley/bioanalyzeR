% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prosize.R
\name{read.prosize.regions}
\alias{read.prosize.regions}
\title{Read a ProSize smear analysis}
\usage{
read.prosize.regions(smear.csv)
}
\arguments{
\item{smear.csv}{A smear analysis CSV exported by ProSize.}
}
\value{
A list containing a list of assay metadata and a data frame of regions (a partial \code{\link{electrophoresis}} object).
}
\description{
This function reads smear analysis table from the ProSize software saved in CSV format. The smear analysis must have been exported in the "alternate" format.
}
\seealso{
\code{\link{read.prosize}}, \code{\link{read.prosize.electropherogram}}, \code{\link{read.prosize.peaks}}
}
