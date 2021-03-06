% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot.R
\name{stdcrv.mobility}
\alias{stdcrv.mobility}
\title{Plot mobility standard curves}
\usage{
stdcrv.mobility(electrophoresis, n.simulate = 100, line.color = "red")
}
\arguments{
\item{electrophoresis}{An \code{\link{electrophoresis}} object.}

\item{n.simulate}{Number of data points to simulate for drawing the standard curve.}

\item{line.color}{Color of the standard curve and data points.}
}
\value{
A ggplot object containing several layers. You can draw it directly or customize it like any other ggplot object by adding more layers.
}
\description{
This function is a shortcut to plot the standard curve(s) of molecule length vs. migration speed from an \code{\link{electrophoresis}} object, wrapping \code{\link{ggplot}}. This allows you to check the quality of the model.
}
\details{
The positions of the ladder peaks reported by the Agilent software are shown in the selected color, and the fluorescence intensites within the peak boundaries are also plotted. If there are multiple ladders, each is shown as a separate facet.
}
\seealso{
\code{\link{qc.electrophoresis}}
}
