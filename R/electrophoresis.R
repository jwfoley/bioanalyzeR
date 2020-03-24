bioanalyzer.first.char <- "<" # XML opening bracket that distinguishes Bioanalyzer XML exports
tapestation.first.char <- rawToChar(as.raw(c(239, 187, 191))) # byte order mark that distinguishes TapeStation XML exports

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
	if (length(arg.list) == 1) return(arg.list[[1]]) # shortcut if only one input
	
 	# increment the sample indexes in the data so they'll match the new table 
	for (i in 1:(length(arg.list) - 1)) for (j in (i + 1):length(arg.list)) {
		arg.list[[j]]$data$sample.index <- arg.list[[j]]$data$sample.index + nrow(arg.list[[i]]$samples)
		if (! is.null(arg.list[[j]]$peaks)) arg.list[[j]]$peaks$sample.index <- arg.list[[j]]$peaks$sample.index + nrow(arg.list[[i]]$samples)
		if (! is.null(arg.list[[j]]$regions)) arg.list[[j]]$regions$sample.index <- arg.list[[j]]$regions$sample.index + nrow(arg.list[[i]]$samples)
	}
	
	structure(list(
		data = do.call(rbind, lapply(arg.list, function(x) x$data)),
		assay.info = do.call(c, lapply(arg.list, function(x) x$assay.info)),
		samples = do.call(rbind, lapply(arg.list, function(x) x$samples)),
		peaks = do.call(rbind, lapply(arg.list, function(x) x$peaks)),
		regions = do.call(rbind, lapply(arg.list, function(x) x$regions)),
		mobility.functions = do.call(c, lapply(arg.list, function(x) x$mobility.functions)),
		mass.coefficients = do.call(c, lapply(arg.list, function(x) x$mass.coefficients))
	), class = "electrophoresis")
}

#' Read files into an electrophoresis object
#'
#' This function reads one or more XML files exported from the Agilent software (and accompanying PNG files if from a TapeStation) and calls the appropriate function to read them into an \code{electrophoresis} object.
#'
#' Multiple input files are combined with \code{\link{rbind.electrophoresis}}.
#'
#' Warning: The filename of each TapeStation gel image must be the same as the corresponding XML file, e.g. \code{run1.png} goes with \code{run1.xml}.
#'
#' @param ... One or more XML files.
#' @param fit The method used to fit the mobility model of molecule length vs. migration distance, one of \code{"interpolation"} (linear interpolation via \code{\link{approxfun}}), \code{"spline"} (splines via \code{\link{splinefun}}), or \code{"regression"} (log-linear regression via \code{\link{lm}} with the model \code{relative.distance ~ log(length)}).
#'
#' @seealso \code{\link{read.bioanalyzer}}, \code{\link{read.tapestation}}
#'
#' @export
read.electrophoresis <- function(..., fit = "spline") do.call(rbind, lapply(list(...), function(xml.file) {
	first.char <- readChar(xml.file, 1)
	if (first.char == bioanalyzer.first.char)
		read.bioanalyzer(xml.file, fit = fit)
	else if (first.char == tapestation.first.char)
		read.tapestation(xml.file, fit = fit)
	else
		stop("unrecognized XML file format")
}))

#' Subset samples an electrophoresis object
#'
#' This function takes a subset of the samples in an \code{electrophoresis} object.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param ... A logical expression indicating samples to keep, and any other arguments passed to \code{\link{subset}}.
#'
#' @return A new \code{electrophoresis} object containing only the data from the subset of samples that match the given expression, with sample indices renumbered and factors releveled.
#'
#' @export
subset.electrophoresis <- function(electrophoresis, ...) {
	nrow.initial <- nrow(electrophoresis$samples)
	electrophoresis$samples <- subset(electrophoresis$samples, ...)
	if (nrow(electrophoresis$samples) == nrow.initial) { # shortcut if all samples are kept
		return(electrophoresis)
	} else if (nrow(electrophoresis$samples) == 0) { # shortcut if no samples are kept
		stop("empty subset")	
	}
	remaining.samples <- as.integer(rownames(electrophoresis$samples))
	
	# remove unwanted data
	electrophoresis$data <- subset(electrophoresis$data, sample.index %in% remaining.samples)
	electrophoresis$assay.info <- electrophoresis$assay.info[names(electrophoresis$assay.info) %in% as.character(electrophoresis$samples$batch)]
	if (! is.null(electrophoresis$peaks)) {
		electrophoresis$peaks <- subset(electrophoresis$peaks, sample.index %in% remaining.samples)
		if (nrow(electrophoresis$peaks) == 0) electrophoresis$peaks <- NULL
	}
	if (! is.null(electrophoresis$regions)) {
		electrophoresis$regions <- subset(electrophoresis$regions, sample.index %in% remaining.samples)
		if (nrow(electrophoresis$regions) == 0) electrophoresis$regions <- NULL
	}
#	electrophoresis$mobility.functions <- 
	if (! is.null(electrophoresis$mass.coefficients)) electrophoresis$mass.coefficients <- electrophoresis$mass.coefficients[remaining.samples]
	
	# renumber sample indices
	new.indices <- 1:nrow(electrophoresis$samples)
	names(new.indices) <- rownames(electrophoresis$samples)
	electrophoresis$data$sample.index <- new.indices[as.character(electrophoresis$data$sample.index)]
	if (! is.null(electrophoresis$peaks)) electrophoresis$peaks$sample.index <- new.indices[as.character(electrophoresis$peaks$sample.index)]
	if (! is.null(electrophoresis$regions)) electrophoresis$regions$sample.index <- new.indices[as.character(electrophoresis$regions$sample.index)]
#	electrophoresis$mobility.functions <-
	
	# rename rows
	rownames(electrophoresis$data) <- 1:nrow(electrophoresis$data)
	rownames(electrophoresis$samples) <- new.indices
	if (! is.null(electrophoresis$peaks)) rownames(electrophoresis$peaks) <- 1:nrow(electrophoresis$peaks)
	if (! is.null(electrophoresis$regions)) rownames(electrophoresis$regions) <- 1:nrow(electrophoresis$regions)
	
	electrophoresis
}

#' Get the original x-variable
#'
#' This function takes an \code{electrophoresis} object and returns the name of the x-value that was used to fit the mobility model.
#'
#' The result should only be either \code{"aligned time"} for Bioanalyzer data or \code{"relative.distance"} for TapeStation data.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#'
#' @return A character giving the name of the x-value.
#'
#' @export
get.x.name <- function(electrophoresis) {
	possible.x.names <- c("aligned.time", "relative.distance")
	result <- possible.x.names[possible.x.names %in% colnames(electrophoresis$data)]
	stopifnot(length(result) == 1)
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
) ! is.na(data[[bound.variable]]) & data[[bound.variable]] >= lower.bound & data[[bound.variable]] <= upper.bound

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
	electrophoresis$data$sample.index == electrophoresis$peaks$sample.index[which.peak] &
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
#' @return A vector of logicals with length \code{nrow(electrophoresis$data)}.
#'
#' @seealso \code{\link{in.peak}}, \code{\link{in.custom.region}}, \code{\link{in.regions}}
#'
#' @export
in.region <- function(electrophoresis, which.region) {
	electrophoresis$data$sample.index == electrophoresis$regions$sample.index[which.region] &
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
#' @seealso \code{\link{in.regions}}, \code{\link{in.peak}}, \code{\link{between.markers}}
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

#' Check whether data points are between markers
#'
#' This function takes an electrophoresis object and reports whether each data point is between the lower and upper length markers.
#'
#' Observations are considered to be between the markers if they are above the upper boundary of the lower marker and the lower boundary of the upper marker, as reported in the Agilent software's peak detection. If there is no upper marker, all points above the upper boundary of the lower marker are considered to be between the markers.
#'
#' @param electrophoresis
#'
#' @return A vector of logicals with length \code{nrow{electrophoresis$data}}.
#'
#' @seealso \code{\link{\in.peaks}}
#'
#' @export
between.markers <- function(electrophoresis) {
	x.name <- get.x.name(electrophoresis)
	result <- rep(F, nrow(electrophoresis$data))
	# first set all points above the lower marker to TRUE
	for (lower.marker in which(electrophoresis$peaks$peak.observations %in% c("Lower Marker", "edited Lower Marker"))) result[electrophoresis$data$sample.index == electrophoresis$peaks$sample.index[lower.marker] & (
		if (x.name == "relative.distance") electrophoresis$data[[x.name]] < electrophoresis$peaks[[paste0("lower.", x.name)]][lower.marker]
		else electrophoresis$data[[x.name]] > electrophoresis$peaks[[paste0("upper.", x.name)]][lower.marker]
	)] <- T
	# then set all points in or above the upper marker, if there is one, to FALSE
	for (upper.marker in which(electrophoresis$peaks$peak.observations %in% c("Upper Marker", "edited Upper Marker"))) result[electrophoresis$data$sample.index == electrophoresis$peaks$sample.index[upper.marker] & (
		if (x.name == "relative.distance") electrophoresis$data[[x.name]] < electrophoresis$peaks[[paste0("upper.", x.name)]][upper.marker]
		else electrophoresis$data[[x.name]] > electrophoresis$peaks[[paste0("lower.", x.name)]][upper.marker]
	)] <- F
	result
}

#' Scale data by a differential
#'
#' Given an x-variable and a y-variable, this function scales the y-values from the observed data points by the differentials of the x-values. The resulting values of y/dx can then be used to make visually accurate graphs.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param x The name of the x-variable in \code{electrophoresis$data}.
#' @param y The name of the y-variable in \code{electrophoresis$data}.
#'
#' @return A vector of the y-values, one for each row of \code{electrophoresis$data}, divided by the differentials of the corresponding x-values.
#'
#' @seealso \code{\link{normalize.proportion}}
#'
#' @export
scale.by.differential <- function(electrophoresis, x, y) {
	stopifnot(all(diff(electrophoresis$data$sample.index) %in% c(0, 1))) # assume data points from each sample are contiguous and ordered by sample
	delta.x <- do.call(c, lapply(unique(electrophoresis$data$sample.index), function(i) c(NA, diff(electrophoresis$data[electrophoresis$data$sample.index == i,x])))) # apply by sample to make sure we don't get a weird delta at the sample boundary
	if (all(delta.x < 0, na.rm = T)) delta.x <- -delta.x else stopifnot(all(delta.x > 0, na.rm = T)) # assume data points are monotonic; if negative (like migration distance) make them positive so the math comes out clean
	
	electrophoresis$data[[y]] / delta.x
}

