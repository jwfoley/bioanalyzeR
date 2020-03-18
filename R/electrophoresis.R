#' Combine multiple electrophoresis objects
#'
#' This function combines multiple \code{electrophoresis} objects into one so you can analyze and graph multiple batches together.
#'
#' All data frames are combined by \code{\link{rbind}} and lists are combined by \code{\link{c}}. Factor levels are expanded to the union of all inputs.
#'
#' @param ... Two or more objects of class \code{electrophoresis}.
#'
#' @return A new \code{electrophoresis} object containing all the data from the previous ones in the provided order.
#'
#' @export
rbind.electrophoresis <- function(...) {
	arg.list <- list(...)
	structure(list(
		data = do.call(rbind, lapply(arg.list, function(x) x$data)),
		assay.info = do.call(c, lapply(arg.list, function(x) x$assay.info)),
		samples = do.call(rbind, lapply(arg.list, function(x) x$samples)),
		wells.by.ladder = do.call(c, lapply(arg.list, function(x) x$wells.by.ladder)),
		peaks = do.call(rbind, lapply(arg.list, function(x) x$peaks)),
		regions = do.call(rbind, lapply(arg.list, function(x) x$regions)),
		mobility.functions = do.call(c, lapply(arg.list, function(x) x$mobility.functions)),
		mass.coefficients = do.call(c, lapply(arg.list, function(x) x$mass.coefficients))
	), class = "electrophoresis")
}

#' Check whether data points are from a certain sample
#'
#' This function takes an electrophoresis object and checks whether each data point is from a specified sample.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param which.sample The integer index of a sample in \code{electrophoresis$samples}.
#'
#' @return A vector of logicals with length \code{nrow(data)}.
#'
#' @seealso \code{\link{from.samples}}, \code{\link{in.peak}}, \code{\link{in.region}}
#'
#' @export
from.sample <- function(electrophoresis, which.sample) {
	electrophoresis$data$batch == electrophoresis$samples$batch[which.sample] &
	electrophoresis$data$well.number == electrophoresis$samples$well.number[which.sample]
}

#' Match data points to the samples they are from
#'
#' This function takes an electrophoresis object and reports which of the samples each data point belongs to.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#'
#' @return A vector of integers with length \code{nrow(data)}. Each element is the integer index of the sample in \code{electrophoresis$samples} that the data point belongs to.
#'
#' @seealso \code{\link{from.sample}}, \code{\link{in.peaks}}, \code{\link{in.regions}}
#'
#' @export
from.samples <- function(electrophoresis) {
	result <- rep(NA, nrow(electrophoresis$data))
	for (sample in 1:nrow(electrophoresis$samples)) result[which(from.sample(electrophoresis, sample))] <- sample
	result
}

#' Check whether data points are within a custom region
#'
#' This function takes the \code{$data} element of an \code{electrophoresis} object and the boundaries of a region, for a desired variable, and checks whether each data point is within that region.
#'
#' @param data A data frame of electrophoresis data, from the \code{$data} member of an \code{electrophoresis} object (not the whole object itself).
#' @param lower.bound Lower boundary of the region.
#' @param upper.bound Upper boundary of the region.
#' @param bound.variable Which variable the boundaries refer to.
#'
#' @return A vector of logicals with length \code{nrow(data)}.
#'
#' @seealso \code{\link{in.region}}, \code{\link{in.peak}}
#'
#' @export
in.custom.region <- function(
	data,
	lower.bound = -Inf,
	upper.bound = Inf,
	bound.variable = "length"
) data[[bound.variable]] >= lower.bound & data[[bound.variable]] <= upper.bound

#' Check whether data points are within a reported peak
#'
#' This function takes an electrophoresis object and checks whether each data point is within a specified peak.
#'
#' Because each peak in the table belongs to a specific sample, only data points from that sample are \code{TRUE}.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param which.peak The integer index of a peak in \code{electrophoresis$peaks}.
#'
#' @return A vector of logicals with length \code{nrow(data)}.
#'
#' @seealso \code{\link{in.region}}, \code{\link{in.custom.region}}, \code{\link{in.peaks}}
#'
#' @export
in.peak <- function(electrophoresis, which.peak) {
	electrophoresis$data$batch == electrophoresis$peaks$batch[which.peak] &
	electrophoresis$data$well.number == electrophoresis$peaks$well.number[which.peak] &
	! is.na(electrophoresis$data$length) &
	electrophoresis$data$length >= electrophoresis$peaks$lower.length[which.peak] &
	electrophoresis$data$length <= electrophoresis$peaks$upper.length[which.peak]
}

#' Check whether data points are within a reported region
#'
#' This function takes an electrophoresis object and checks whether each data point is within a specified region.
#'
#' Because each region in the table belongs to a specific sample, only data points from that sample are \code{TRUE}.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param which.region The integer index of a region in \code{electrophoresis$regions}.
#'
#' @return A vector of logicals with length \code{nrow(data)}.
#'
#' @seealso \code{\link{in.peak}}, \code{\link{in.custom.region}}, \code{\link{in.regions}}
#'
#' @export
in.region <- function(electrophoresis, which.region) {
	electrophoresis$data$batch == electrophoresis$regions$batch[which.region] &
	electrophoresis$data$well.number == electrophoresis$regions$well.number[which.region] &
	! is.na(electrophoresis$data$length) &
	electrophoresis$data$length >= electrophoresis$regions$lower.length[which.region] &
	electrophoresis$data$length <= electrophoresis$regions$upper.length[which.region]
}

#' Match data points to the peaks they are in
#'
#' This function takes an electrophoresis object and reports which of the peaks each data point belongs to, if any.
#'
#' Warning: If peaks in the reported table overlap, any data point in more than one peak will only be matched to the last peak it belongs to.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#'
#' @return A vector of integers with length \code{nrow(data)}. Each element is either the integer index of the peak in \code{electrophoresis$peaks} that the data point belongs to, or NA if it is not in any of the annotated peaks.
#'
#' @seealso \code{\link{in.regions}}, \code{\link{in.peak}}
#'
#' @export
in.peaks <- function(electrophoresis) {
	result <- rep(NA, nrow(electrophoresis$data))
	if (! is.null(electrophoresis$peaks)) for (i in 1:nrow(electrophoresis$peaks)) result[which(in.peak(electrophoresis, i))] <- i
	result
}

#' Match data points to the regions they are in
#'
#' This function takes an electrophoresis object and reports which of the regions each data point belongs to, if any.
#'
#' Warning: If regions in the reported table overlap, any data point in more than one region will only be matched to the last region it belongs to.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#'
#' @return A vector of integers with length \code{nrow(data)}. Each element is either the integer index of the region in \code{electrophoresis$regions} that the data point belongs to, or NA if it is not in any of the annotated regions.
#'
#' @seealso \code{\link{in.peaks}}, \code{\link{in.region}}
#'
#' @export
in.regions <- function(electrophoresis) {
	result <- rep(NA, nrow(electrophoresis$data))
	if (! is.null(electrophoresis$regions)) for (i in 1:nrow(electrophoresis$regions)) result[which(in.region(electrophoresis, i))] <- i
	result
}

