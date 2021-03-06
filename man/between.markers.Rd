% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/electrophoresis.R
\name{between.markers}
\alias{between.markers}
\title{Check whether data points are between markers}
\usage{
between.markers(electrophoresis, lower.marker.spread = 10)
}
\arguments{
\item{electrophoresis}{An \code{\link{electrophoresis}} object.}

\item{lower.marker.spread}{Proportion to scale the width of the lower marker peak, along the computed length scale, to compensate for underreporting in the Agilent software. Set to 1 to use the reported peak boundary, which works poorly.}
}
\value{
A vector of logicals with length \code{nrow{electrophoresis$data}}.
}
\description{
This function takes an electrophoresis object and reports whether each data point is between the lower and upper length markers.
}
\details{
Observations are considered to be between the markers if they are above the upper boundary of the lower marker and the lower boundary of the upper marker, as reported in the Agilent software's peak detection. If there is no upper marker, all points above the upper boundary of the lower marker are considered to be between the markers unless they have no estimated length (i.e. are beyond the last ladder peak).

The peak boundaries reported by Agilent tend to be too narrow for the lower marker and leave some residual fluorescence that can greatly distort some calculations, so by default this function overrides the 

Note: Data exported from the ProSize software are missing the peak boundaries, so in that situation only the precise lengths of the marker peaks are set as the boundaries. The inner half of each marker peak will still be included in the result.
}
\seealso{
\code{\link{in.peaks}}, \code{\link{in.regions}}
}
