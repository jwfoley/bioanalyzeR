#' Integrate a variable under the electrophoresis curve
#'
#' Compute the sum of some electrophoresis variable between some specified boundaries, or within an annotated peak, or both (intersection).
#'
#' This is an all-purpose base function that accepts any combination of sum variable, boundaries, boundary variable, and peak. Crucially, it operates directly on the \code{$data} member of an \code{electrophoresis} object rather than the entire object, and does not separate the results by sample. You can achieve this manually by taking a \code{\link{subset}} of the data but it is probably easier to use one of the functions that wrap this one.
#'
#' @param data A data frame of electrophoresis data, from the \code{$data} member of an electrophoresis object (not the whole object itself).
#' @param lower.bound Lower boundary of the region to integrate, or \code{NULL} to extend to negative infinity.
#' @param upper.bound Upper boundary of the region to integrate, or \code{NULL} to extend to positive infinity.
#' @param bound.variable Which variable the boundaries refer to, e.g. \code{"length"}.
#' @param peak The index of the peak to integrate in, or \code{NULL} to ignore peaks.
#' @param sum.variable Which variable to sum in the target region.
#'
#' @seealso \code{\link{integrate.peaks}}, \code{\link{integrate.regions}}, \code{\link{integrate.custom}}
#'
#' @export
integrate.rawdata <- function(
	data,
	lower.bound = NULL,
	upper.bound = NULL,
	bound.variable = "length",
	peak = NULL,
	sum.variable = "molarity"
) {
	above.lower <- if (is.null(lower.bound)) rep(TRUE, nrow(data)) else data[[bound.variable]] >= lower.bound
	below.upper <- if (is.null(upper.bound)) rep(TRUE, nrow(data)) else data[[bound.variable]] <= upper.bound
	in.peak <- if (is.null(peak)) rep(TRUE, nrow(data)) else ! is.na(data$peak) & data$peak == peak
	sum(data[[sum.variable]][above.lower & below.upper & in.peak])
}

#' Integrate a variable under each peak
#'
#' Compute the sum of some electrophoresis variable between the boundaries of each reported peak in an \code{electrophoresis} object.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param sum.variable Which variable to sum in each peak.
#'
#' @seealso \code{\link{integrate.rawdata}}, \code{\link{integrate.regions}}, \code{\link{integrate.custom}}
#'
#' @export
integrate.peaks <- function(
	electrophoresis,
	sum.variable = "molarity"
) sapply(1:nrow(electrophoresis$peaks), function(peak) integrate.rawdata(electrophoresis$data, peak = peak, sum.variable = sum.variable))

#' Integrate a variable in each region
#'
#' Compute the sum of some electrophoresis variable between the boundaries of each region of interest in an \code{electrophoresis} object.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param sum.variable Which variable to sum in each region.
#'
#' @seealso \code{\link{integrate.rawdata}}, \code{\link{integrate.peaks}}, \code{\link{integrate.custom}}
#'
#' @export
integrate.regions <- function(
	electrophoresis,
	sum.variable = "molarity"
) sapply(1:nrow(electrophoresis$region), function(region.index) integrate.rawdata(subset(electrophoresis$data, well.number == as.character(electrophoresis$regions$well.number[region.index])), electrophoresis$regions$lower.length[region.index], electrophoresis$regions$upper.length[region.index], sum.variable = sum.variable))

#' Integrate a variable in a custom region
#'
#' Compute the sum of some electrophoresis variable in an \code{electrophoresis} object between specified boundaries. The variable is summed individually for each sample.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param ... Arguments passed to \code{\link{integrate.rawdata}}. In particular use \code{lower.bound} and \code{upper.bound} to set the boundaries, \code{bound.variable} for the boundary variable, and \code{sum.variable} for the variable to sum.
#'
#' @seealso \code{\link{integrate.rawdata}}, \code{\link{integrate.peaks}}, \code{\link{integrate.regions}}, \code{\link{region.ratio}}
#'
#' @export
integrate.custom <- function(
	electrophoresis,
	...
) sapply(as.character(electrophoresis$samples$well.number), function(well) integrate.rawdata(subset(electrophoresis$data, well.number == well), ...))

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

