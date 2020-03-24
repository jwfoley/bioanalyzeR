#' Integrate a variable under each peak
#'
#' Compute the sum of some electrophoresis variable between the boundaries of each reported peak in an \code{electrophoresis} object.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param sum.variable Which variable to sum in each peak.
#'
#' @seealso \code{\link{integrate.regions}}, \code{\link{integrate.custom}}
#'
#' @export
integrate.peaks <- function(
	electrophoresis,
	sum.variable = "molarity"
) sapply(1:nrow(electrophoresis$peaks), function(peak) sum(electrophoresis$data[[sum.variable]][which(in.peak(electrophoresis, peak))]))

#' Integrate a variable in each region
#'
#' Compute the sum of some electrophoresis variable between the boundaries of each region of interest in an \code{electrophoresis} object.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param sum.variable Which variable to sum in each region.
#'
#' @seealso \code{\link{integrate.peaks}}, \code{\link{integrate.custom}}
#'
#' @export
integrate.regions <- function(
	electrophoresis,
	sum.variable = "molarity"
) sapply(1:nrow(electrophoresis$regions), function(region) sum(electrophoresis$data[[sum.variable]][which(in.region(electrophoresis, region))]))

#' Integrate a variable in a custom region
#'
#' Compute the sum of some electrophoresis variable in an \code{electrophoresis} object between specified boundaries. The variable is summed individually for each sample.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param lower.bound Lower boundary of the region to integrate.
#' @param upper.bound Upper boundary of the region to integrate.
#' @param bound.variable Which variable the boundaries refer to.
#' @param sum.variable Which variable to sum in each region.
#'
#' @seealso \code{\link{integrate.peaks}}, \code{\link{integrate.regions}}, \code{\link{region.ratio}}
#'
#' @export
integrate.custom <- function(
	electrophoresis,
	lower.bound = -Inf,
	upper.bound = Inf,
	bound.variable = "length",
	sum.variable = "molarity"
) {
	in.this.region <- in.custom.region(electrophoresis$data, lower.bound, upper.bound, bound.variable)
	sapply(1:nrow(electrophoresis$samples), function(sample) sum(electrophoresis$data[[sum.variable]][in.this.region & electrophoresis$data$sample.index == sample]))
}

#' Compare sums within regions
#'
#' Given two or more regions (pairs of lower and upper bounds), calculate the ratio of the integrated sum of each additional region relative to the integrated sum of the first region.  The ratio is compuate individually for each sample.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param bounds A list of two or more pairs (vectors) of boundaries, e.g. \code{list(c(100, 200), c(200, 500), c(500, 700)}.
#' @param bound.variable Which variable the boundaries refer to.
#' @param sum.variable Which variable to sum in each region.
#'
#' @return A matrix of ratios of sums within the regions, each region relative to the first region, for each sample.
#'
#' @seealso \code{\link{dv200}}, \code{\link{illumina.library.ratio}}
#'
#' @export
region.ratio <- function(
	electrophoresis,
	bounds,
	bound.variable = "length",
	sum.variable = "molarity"
) {
	stopifnot(length(bounds) >= 1)
	sum.matrix <- sapply(bounds, function(bound.pair) integrate.custom(electrophoresis, lower.bound = bound.pair[1], upper.bound = bound.pair[2], bound.variable = bound.variable, sum.variable = sum.variable))
	matrix(sum.matrix[,-1] / sum.matrix[,1], dimnames = list(NULL, sapply(bounds[-1], function(bound.pair) paste0(sum.variable, " ratio in ", bound.variable, " ", bound.pair[1], "-", bound.pair[2], "/", bounds[[1]][1], "-", bounds[[1]][2]))))
}

#' DV200 analysis
#'
#' Calculate the proportion of fragments above 200 bases (DV200). Only fragments between the lower and upper markers are considered.
#'
#' Note: Despite unclear wording, Agilent calculates DV200 as a proportion of total mass (concentration), rather than a proportion of molecules (molarity). For some purposes molar DV200 may be more relevant, but current protocols refer to the default DV200 based on mass.
#'
#' By default \code{lower.marker.spread = 1} (see \code{\link{between.markers}}) for consistency with the Agilent software.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param prop.variable Which variable to use for the proportion.
#' @param lower.marker.spread Amount to scale the width of the lower marker peak.
#'
#' @seealso \code{\link{region.ratio}}, \code{\link{illumina.library.ratio}}
#'
#' @export
dv200 <- function(electrophoresis, prop.variable = "concentration", lower.marker.spread = 1) region.ratio(electrophoresis$data[which(between.markers(electrophoresis, lower.marker.spread)),], bounds = list(c(-Inf, Inf), c(200, Inf)), sum.variable = prop.variable)


#' Ratio of good inserts to adapter dimers
#'
#' For Illumina sequencing libraries, compute the molar ratio of molecules with desirably long inserts to undesirable adapter dimers in each sample.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param min.sequenceable The shortest length of molecules that are likely to produce data on the sequencer (ignore shorter molecules that are probably just unclusterable free primers or adapters).
#' @param min.good.insert The shortest length of a desirable molecule (adapter length plus insert length).
#' @param max.sequenceable The longest length of molecules that are likely to produce data on the sequencer (ignore molecules that are too long for clustering).
#'
#' @seealso \code{\link{region.ratio}}, \code{\link{dv200}}
#'
#' @export
illumina.library.ratio <- function(
	electrophoresis,
	min.sequenceable =  100,
	min.good.insert =   200,
	max.sequenceable =  700
) region.ratio(electrophoresis, bounds = list(c(min.sequenceable, min.good.insert), c(min.good.insert, max.sequenceable)))


#' Normalize data to proportions
#'
#' For a given variable, normalize the data values in an \code{electrophoresis} object to proportions of the total. After normalization the sum of the variable for each sample is 1.
#'
#' Only data points between the markers are considered in the total. If there is no upper marker, all data points above the lower marker are considered.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param variable The name of a variable in \code{electrophoresis$data}.
#'
#' @return A vector whose length equals \code{nrow(electrophoresis$data)}, in which each value is an observation of the given variable normalized by the total for that sample, or NA if the observation is not between the markers.
#'
#' @seealso \code{\link{scale.by.differential}}
#'
#' @export
normalize.proportion <- function(electrophoresis, variable) {
	which.usable <- which(between.markers(electrophoresis) & ! is.na(electrophoresis$data[[variable]]))
	subset.usable <- electrophoresis$data[which.usable,]
	sample.sums <- sapply(1:nrow(electrophoresis$samples), function(sample) sum(subset.usable[[variable]][subset.usable$sample.index == sample]))
	result <- rep(NA, nrow(electrophoresis$data))
	result[which.usable] <- subset.usable[[variable]] / sample.sums[subset.usable$sample.index]
	result
}

