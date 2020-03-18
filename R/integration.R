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
	sapply(1:nrow(electrophoresis$samples), function(sample) sum(electrophoresis$data[[sum.variable]][in.this.region & from.sample(electrophoresis, sample)]))
}

#' Compare sums within regions
#'
#' Given two or more regions (pairs of lower and upper bounds), calculate the ratio of the integrated sum of each additional region relative to the integrated sum of the first region.  The ratio is compuate individually for each sample.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param bounds A list of two or more pairs (vectors) of boundaries, e.g. \code{list(c(100, 200), c(200, 500), c(500, 700)}.
#' @param ... Additional arguments passed to \code{\link{integrate.rawdata}}. In particular use \code{bound.variable} to set the boundary variable and \code{sum.variable} for the variable to sum.
#'
#' @return A matrix of ratios of sums within the regions, each region relative to the first region, for each sample.
#'
#' @seealso \code{\link{illumina.library.ratio}}
#'
#' @export
region.ratio <- function(
	electrophoresis,
	bounds,
	...
) {
	stopifnot(length(bounds) >= 1)
	sum.matrix <- sapply(bounds, function(bound.pair) integrate.custom(electrophoresis, lower.bound = bound.pair[1], upper.bound = bound.pair[2], ...))
	result <- as.matrix(sum.matrix[,-1] / sum.matrix[,1])
	colnames(result) <- sapply(bounds[-1], function(bound.pair) paste0(bound.pair[1], "-", bound.pair[2], " / ", bounds[[1]][1], "-", bounds[[1]][2]))
	result
}

#' Ratio of good inserts to adapter dimers
#'
#' For Illumina sequencing libraries, compute the molar ratio of molecules with desirably long inserts to undesirable adapter dimers in each sample.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param min.sequenceable The shortest length of molecules that are likely to produce data on the sequencer (ignore shorter molecules that are probably just unclusterable free primers or adapters).
#' @param min.good.insert The shortest length of a desirable molecule (adapter length plus insert length).
#' @param max.sequenceable The longest length of molecules that are likely to produce data on the sequencer (ignore molecules that are too long for clustering).
#'
#' @export
illumina.library.ratio <- function(
	electrophoresis,
	min.sequenceable =  100,
	min.good.insert =   200,
	max.sequenceable =  700
) region.ratio(electrophoresis, bounds = list(c(min.sequenceable, min.good.insert), c(min.good.insert, max.sequenceable)))

