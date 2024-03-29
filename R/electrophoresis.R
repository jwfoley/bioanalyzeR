# file type identifiers
BIOANALYZER.FIRST.CHAR <- "<" # XML opening bracket that distinguishes Bioanalyzer XML exports
TAPESTATION.FIRST.CHAR <- rawToChar(as.raw(239)) # first byte of the byte order mark that distinguishes TapeStation XML exports
GZIP.FIRST.CHAR <- rawToChar(as.raw(31)) # first byte of the gzip magic number

# marker peak identifiers
LOWER.MARKER.NAMES <- c("Lower Marker", "edited Lower Marker")
UPPER.MARKER.NAMES <- c("Upper Marker", "edited Upper Marker")

#' Electrophoresis class
#'
#' This S3 class is a generic container for electrophoresis data and metadata. The constructor simply assembles the provided members into a class structure without checking them: if they are formatted incorrectly there will be no warnings at this stage.
#'
#' @param data A tall data frame of the run data, specifically:
#' * `time` (Bioanalyzer, ProSize) - time when this data point was measured
#' * `aligned.time` (Bioanalyzer, ProSize) - measurement time aligned between the expected times of the marker peaks
#' * `distance` (TapeStation) - migration distance of the measurement from the top of the gel area
#' * `relative.distance` (TapeStation) - migration distance normalized relative to the marker peaks
#' * `fluorescence` - fluorescence reading at this point
#' * `length` - estimated molecule length at this point
#' * `concentration` - estimated concentration of the area under the curve between this point and the previous one
#' * `molarity` - estimated molarity of the area under the curve between this point and the previous one
#' @param assay.info A list of metadata about each batch and the assay kit used.
#' @param samples A data frame of metadata for each sample (also annotated to `data`, `peaks`, and `regions` with the same factor levels), specifically:
#' * `batch` - the batch (instrument run) of the sample, from the file name
#' * `well.number` - the well number in which the sample was loaded
#' * `sample.name` - the name of the sample
#' * `sample.observations` - notes about this sample supplied by the user or the Agilent software
#' * `sample.comment` - notes about this sample supplied by the user
#' * `reagent.id` (TapeStation) - the name of the ScreenTape used for this sample
#' * `ladder.well` - which well contains the ladder that calibrates this sample
#' * `RIN`, `DIN`, `RQN`, `28S/18S`, etc. - sample quality metrics reported by the Agilent software
#' @param peaks A data frame of peaks reported by the Agilent software, annotated with their lower and upper boundaries in various scales.
#' @param regions A data frame of regions of interest reported by the Agilent software, annotated with their lower and upper boundaries in varous scales.
#' @param calibration A list of mobility functions to convert migration speed measurements (aligned time or relative distance) into estimated molecule lengths. For each ladder in each batch, there is a list containing:
#' * `mobility.function` - the function that estimates molecule lengths from migration speed measurements
#' * `mobility.inverse` - a function fit by the same method for estimating the migration speed measurement of a given molecule length
#' * `ladder.peaks` - the known lengths and migration speed measurements of this ladder's peaks used for calibration
#'
#' @export electrophoresis
#' @exportClass electrophoresis
#' @md
electrophoresis <- function(
	data = NULL,
	assay.info = NULL,
	samples = NULL,
	peaks = NULL,
	regions = NULL,
	calibration = NULL
) structure(list(
	data = data,
	assay.info = assay.info,
	samples = samples,
	peaks = peaks,
	regions = regions,
	calibration = calibration
), class = "electrophoresis") 


#' Combine multiple electrophoresis objects
#'
#' This function combines multiple \code{\link{electrophoresis}} objects into one so you can analyze and graph multiple batches together.
#'
#' All data frames are combined by \code{\link{rbind}} and lists are combined by \code{\link{c}}. Factor levels are expanded to the union of all inputs.
#'
#' @param ... Two or more objects of class \code{\link{electrophoresis}}.
#'
#' @return A new \code{\link{electrophoresis}} object containing all the data from the previous ones in the provided order.
#'
#' @export
#' @importFrom plyr rbind.fill
rbind.electrophoresis <- function(...) {
	arg.list <- list(...)
	if (length(arg.list) == 1) return(arg.list[[1]]) # shortcut if only one input
	
 	# increment the sample indexes in the data so they'll match the new table 
	for (i in 1:(length(arg.list) - 1)) for (j in (i + 1):length(arg.list)) {
		arg.list[[j]]$data$sample.index <- arg.list[[j]]$data$sample.index + nrow(arg.list[[i]]$samples)
		if (! is.null(arg.list[[j]]$peaks)) arg.list[[j]]$peaks$sample.index <- arg.list[[j]]$peaks$sample.index + nrow(arg.list[[i]]$samples)
		if (! is.null(arg.list[[j]]$regions)) arg.list[[j]]$regions$sample.index <- arg.list[[j]]$regions$sample.index + nrow(arg.list[[i]]$samples)
	}
	
	electrophoresis(
		data = do.call(rbind.fill, lapply(arg.list, function(x) x$data)),
		assay.info = do.call(c, lapply(arg.list, function(x) x$assay.info)),
		samples = do.call(rbind.fill, lapply(arg.list, function(x) x$samples)),
		peaks = do.call(rbind.fill, lapply(arg.list, function(x) x$peaks)),
		regions = do.call(rbind.fill, lapply(arg.list, function(x) x$regions)),
		calibration = do.call(c, lapply(arg.list, function(x) x$calibration))
	)
}


#' Read files into an electrophoresis object
#'
#' These functions read one or more XML, CSV, or ZIP files exported from the Agilent software (and accompanying unaligned electropherogram CSV files if from a TapeStation) and calls the appropriate function to read them into an \code{\link{electrophoresis}} object, which is filled out with estimates of molecule length, concentration, and molarity (see \code{\link{calibrate.electrophoresis}}). \code{read.electrophoresis} is the easiest to use as it automatically infers the correct file type.
#'
#' @param xml.file The filename of an XML file exported from the Bioanalyzer or TapeStation software. The XML file may be compressed with `gzip` and the filename can be a remote URL. The filename is expected to end in \code{.xml} or \code{.xml.gz} and the name before that extension is used as the name of the batch. A TapeStation XML file must be accompanied by a CSV file, whose path is inferred from the XML file's if its name is in the expected format.
#' @param ... One or more XML files exported from the Bioanalyzer or TapeStation software, or CSV or ZIP files exported from the ProSize software.
#' @param mc.cores Maximum number of CPU cores to use (passed to \code{\link[parallel]{mclapply}}). Only one core is used per input file.
#'
#' @inheritParams read.bioanalyzer
#' @inheritParams read.tapestation
#' @inheritParams read.prosize
#' @inheritParams read.prosize.zip
#'
#' @name read.electrophoresis
#' 
#' @export
#' @importFrom parallel mclapply detectCores
read.electrophoresis <- function(
	...,
	method = "hyman",
	extrapolate = FALSE,
	mc.cores = if (.Platform$OS.type == "windows") 1 else detectCores()
) do.call(rbind, mclapply(list(...), function(file.path) {
	if (endsWith(file.path, ".csv")) {
		read.prosize(file.path, method = method, extrapolate = extrapolate)
	} else if (endsWith(file.path, ".zip")) {
		read.prosize.zip(file.path, method = method, extrapolate = extrapolate)
	} else {
		xml.con <- file(file.path)
		first.char <- readChar(xml.con, 1, useBytes = T)
		if (first.char == GZIP.FIRST.CHAR) first.char <- readChar(gzcon(xml.con), 1, useBytes = T) # if gzipped, uncompress and try again
		close(xml.con) # if not explicitly closed, R gives a warning
		if (first.char == BIOANALYZER.FIRST.CHAR)
			read.bioanalyzer(file.path, method = method, extrapolate = extrapolate)
		else if (first.char == TAPESTATION.FIRST.CHAR)
			read.tapestation(file.path, method = method, extrapolate = extrapolate)
		else
			stop("unrecognized XML file format")
	}
}, mc.cores = mc.cores))


#' Add annotations to an electrophoresis object
#'
#' This function adds columns of annotations to the \code{$samples} table of an \code{\link{electrophoresis}} object. These new annotations can be used for subsetting, plot faceting, etc.
#'
#' The input annotations can be a data frame or either a \code{\link{connection}} or a \code{\link{character}} containing the path of a file, which is then read into a data frame by \code{\link{read.table}}. Therefore a file will probably need to be in CSV format or similar.
#'
#' The first column of \code{annotations} is used to identify samples in the matching column of \code{electrophoresis$samples}, so it needs a matching column name, usually \code{sample.name} or \code{well.number}. But you can use any column that exists in the table, such as \code{batch} or \code{reagent.id} or even a column added by a previous call to this function. The elements of this first column are interpreted as class \code{\link{character}} and cannot contain any duplicates. However, this column does not need to be in the same order as \code{electrophoresis$samples} and not all identifiers in the sample table (e.g. sample names) need to appear in the annotation table, nor vice versa.
#'
#' Each additional column of \code{annotations} produces a new column of the same name in \code{electrophoresis$samples}. If there is already a column of that name, it is overwritten. Samples that do not appear in the annotation table receive \code{NA} for their new annotations. Samples with the same identifier (e.g. the same \code{sample.name} because they are replicates) receive the same annotations.
#'
#' If \code{stringsAsFactors} is true then variables parsed as strings will be cast as factors, but the factor levels will be kept in the order they appear in the table, unlike the behavior of \code{\link{read.table}}, which alphabetizes them. Thus the order in the annotation file determines the order of factor levels.
#'
#' @param electrophoresis An \code{\link{electrophoresis}} object.
#' @param annotation Either a data frame or the path of a file that can be read by \code{\link{read.table}}.
#' @param ... Additional arguments passed to \code{\link{read.table}}.
#'
#' @return A new \code{\link{electrophoresis}} object with additional columns added to its \code{$samples} element.
#'
#' @export
annotate.electrophoresis <- function(
	electrophoresis,
	annotations,
	header = TRUE,
	row.names = NULL,
	sep = "\t",
	stringsAsFactors = TRUE,
	...
) {
	if (any(class(annotations) %in% c("character", "connection"))) annotations <- read.table(
		annotations,
		header = header,
		row.names = row.names,
		sep = sep,
		stringsAsFactors = F, # cast them as factors later
		...
	)
	stopifnot(
		"empty annotations" = ncol(annotations) > 1,
		"no label recognized in annotations" = colnames(annotations)[1] %in% colnames(electrophoresis$samples),
		"duplicate annotations" = anyDuplicated(annotations[,1]) == 0
	)
	
	identifiers <- as.character(electrophoresis$samples[,colnames(annotations)[1]])
	for (col in 2:ncol(annotations)) {
		annotation.lookup <- annotations[,col]
		if (stringsAsFactors && class(annotation.lookup) == "character") annotation.lookup <- factor(annotation.lookup, levels = unique(annotation.lookup))
		names(annotation.lookup) <- annotations[,1]
		electrophoresis$samples[,colnames(annotations)[col]] <- annotation.lookup[identifiers]
	}
	
	electrophoresis
}


#' Subset samples an electrophoresis object
#'
#' This function takes a subset of the samples in an \code{\link{electrophoresis}} object.
#'
#' @param electrophoresis An \code{\link{electrophoresis}} object.
#' @param ... A logical expression indicating samples to keep, and any other arguments passed to \code{\link{subset}}.
#'
#' @return A new \code{\link{electrophoresis}} object containing only the data from the subset of samples that match the given expression, with sample indices renumbered and factors releveled. Note: mobility calibrations from unused samples are kept.
#'
#' @export
subset.electrophoresis <- function(electrophoresis, ...) {
	nrow.initial <- nrow(electrophoresis$samples)
	electrophoresis$samples <- subset(electrophoresis$samples, ...)
	if (nrow(electrophoresis$samples) == nrow.initial) { # shortcut if all samples are kept
		return(electrophoresis)
	} else if (nrow(electrophoresis$samples) == 0) { # shortcut if no samples are kept
		for (member in names(electrophoresis)) electrophoresis[[member]] <- NULL
		return(electrophoresis)
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
	
	# renumber sample indices
	new.indices <- setNames(seq(nrow(electrophoresis$samples)), rownames(electrophoresis$samples))
	electrophoresis$data$sample.index <- new.indices[as.character(electrophoresis$data$sample.index)]
	if (! is.null(electrophoresis$peaks)) electrophoresis$peaks$sample.index <- new.indices[as.character(electrophoresis$peaks$sample.index)]
	if (! is.null(electrophoresis$regions)) electrophoresis$regions$sample.index <- new.indices[as.character(electrophoresis$regions$sample.index)]
	
	# rename rows
	rownames(electrophoresis$data) <- seq(nrow(electrophoresis$data))
	rownames(electrophoresis$samples) <- new.indices
	if (! is.null(electrophoresis$peaks)) rownames(electrophoresis$peaks) <- seq(nrow(electrophoresis$peaks))
	if (! is.null(electrophoresis$regions)) rownames(electrophoresis$regions) <- seq(nrow(electrophoresis$regions))
	
	electrophoresis
}


#' Get the original x-variable
#'
#' This function takes an \code{\link{electrophoresis}} object and returns the name of the x-variable that was used to fit the mobility model.
#'
#' If `raw == FALSE` the result should only be either \code{"aligned time"} for Bioanalyzer data or \code{"relative.distance"} for TapeStation data. If `raw == TRUE` the result should be \code{"time"} or \code{"distance"}.
#'
#' @param electrophoresis An \code{\link{electrophoresis}} object.
#' @param raw Whether to return the name of the raw variable instead of the aligned variable.
#' @param allow.multiple Whether to allow multiple raw variables (data combined from multiple platforms).
#'
#' @return A character giving the name of the x-variable.
#'
#' @export
get.x.name <- function(electrophoresis, raw = FALSE, allow.multiple = FALSE) {
	possible.x.names <- if (raw) c("time", "distance") else c("aligned.time", "relative.distance")
	result <- intersect(possible.x.names, colnames(electrophoresis$data))
	stopifnot(
		"no x-variable found" = length(result) > 0,
		"multiple x-variables" = allow.multiple || length(result) == 1
	)
	result
}


#' Check whether data points are within a custom region
#'
#' This function takes the \code{$data} element of an \code{\link{electrophoresis}} object and the boundaries of a region, for a desired variable, and checks whether each data point is within that region.
#'
#' @param data A data frame of electrophoresis data, from the \code{$data} member of an \code{\link{electrophoresis}} object (not the whole object itself).
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


#' Check whether data points are within a reported peak or region
#' 
#' These functions take an electrophoresis object and check whether each data point is within a specified peak or region.
#'
#' Data points are considered to be within the boundaries of a peak if their original x-value (`aligned.time` for Bioanalyzer, `relative.distance` for TapeStation) is within the boundaries reported by the Agilent software, but for regions the length boundaries are used.
#'
#' @param electrophoresis An \code{\link{electrophoresis}} object.
#' @param which.peak The integer index of a peak in \code{electrophoresis$peaks}.
#' @param which.region The integer index of a region in \code{electrophoresis$regions}.
#' 
#' @return A vector of logicals with length \code{nrow(electrophoresis$data)}.
#'
#' @seealso \code{\link{in.custom.region}}, \code{\link{in.peaks}}, \code{\link{in.regions}}
#'
#' @name in.peak.region
NULL


#' @rdname in.peak.region
#' @export
in.peak <- function(electrophoresis, which.peak) {
	x.names <- get.x.name(electrophoresis, allow.multiple = T)
	peak.x.name <- x.names[which(! is.na(electrophoresis$peaks[which.peak, x.names]))]
	if (length(peak.x.name) == 0) { # no non-NA peak x-value found so nothing is in it
		rep(NA, nrow(electrophoresis$data))
	} else {
		stopifnot("multiple peak x-values" = length(peak.x.name) == 1)
		electrophoresis$data$sample.index == electrophoresis$peaks$sample.index[which.peak] &
			electrophoresis$data[[peak.x.name]] >= electrophoresis$peaks[[paste0("lower.", peak.x.name)]][which.peak] &
			electrophoresis$data[[peak.x.name]] <= electrophoresis$peaks[[paste0("upper.", peak.x.name)]][which.peak]
	}
}


#' @rdname in.peak.region
#' @export
in.region <- function(electrophoresis, which.region) {
	electrophoresis$data$sample.index == electrophoresis$regions$sample.index[which.region] &
	! is.na(electrophoresis$data$length) &
	electrophoresis$data$length >= electrophoresis$regions$lower.length[which.region] &
	electrophoresis$data$length <= electrophoresis$regions$upper.length[which.region]
}


#' Match data points to the peaks or regions they are in
#'
#' These functions take an electrophoresis object and report which of the peaks or regions each data point belongs to, if any.
#'
#' Data points are considered to be within the boundaries of a peak if their original x-value (`aligned.time` for Bioanalyzer, `relative.distance` for TapeStation) is within the boundaries reported by the Agilent software, but for regions the length boundaries are used.
#'
#' Warning: If peaks or regions in the reported table overlap, any data point in more than one peak or region will only be matched to the last peak or region it belongs to.
#'
#' @param electrophoresis An \code{\link{electrophoresis}} object.
#'
#' @return A vector of integers with length \code{nrow(data)}. Each element is either the integer index of the peak in \code{electrophoresis$peaks} or region in \code{electrophoresis$regions} that the data point belongs to, or NA if it is not in any of the annotated peaks or regions.
#'
#' @seealso \code{\link{in.peak}}, \code{\link{in.region}}, \code{\link{between.markers}}
#'
#' @name in.peaks.regions
NULL


#' @rdname in.peaks.regions
#' @export
in.peaks <- function(electrophoresis) {
	result <- rep(NA, nrow(electrophoresis$data))
	if (! is.null(electrophoresis$peaks)) for (i in seq(nrow(electrophoresis$peaks))) result[which(in.peak(electrophoresis, i))] <- i
	result
}


#' @rdname in.peaks.regions
#' @export
in.regions <- function(electrophoresis) {
	result <- rep(NA, nrow(electrophoresis$data))
	if (! is.null(electrophoresis$regions)) for (i in seq(nrow(electrophoresis$regions))) result[which(in.region(electrophoresis, i))] <- i
	result
}


#' Check whether data points are between markers
#'
#' This function takes an electrophoresis object and reports whether each data point is between the lower and upper length markers.
#'
#' Observations are considered to be between the markers if they are above the upper boundary of the lower marker and the lower boundary of the upper marker, as reported in the Agilent software's peak detection. If there is no upper marker, all points above the upper boundary of the lower marker are considered to be between the markers unless they have no estimated length (i.e. are beyond the last ladder peak).
#'
#' The peak boundaries reported by Agilent tend to be too narrow for the lower marker and leave some residual fluorescence that can greatly distort some calculations, so by default this function overrides the 
#'
#' Note: Data exported from the ProSize software are missing the peak boundaries, so in that situation only the precise lengths of the marker peaks are set as the boundaries. The inner half of each marker peak will still be included in the result.
#'
#' @param electrophoresis An \code{\link{electrophoresis}} object.
#' @param lower.marker.spread Proportion to scale the width of the lower marker peak, along the computed length scale, to compensate for underreporting in the Agilent software. Set to 1 to use the reported peak boundary, which works poorly.
#'
#' @return A vector of logicals with length \code{nrow{electrophoresis$data}}.
#'
#' @seealso \code{\link{in.peaks}}, \code{\link{in.regions}}
#'
#' @export
between.markers <- function(electrophoresis, lower.marker.spread = 10) {
	result <- rep(F, nrow(electrophoresis$data))
	# first set all points above the lower marker to TRUE
	for (lower.marker in which(electrophoresis$peaks$peak.observations %in% LOWER.MARKER.NAMES)) {
		lower.bound <- if (is.na(electrophoresis$peaks$upper.length[lower.marker])) electrophoresis$peaks$length[lower.marker] else lower.marker.spread * (electrophoresis$peaks$upper.length[lower.marker] - electrophoresis$peaks$length[lower.marker])
		result[
			electrophoresis$data$sample.index == electrophoresis$peaks$sample.index[lower.marker] & 
			electrophoresis$data$length > lower.bound
		] <- T
	}
	# then set all points in or above the upper marker, if there is one, to FALSE
	for (upper.marker in which(electrophoresis$peaks$peak.observations %in% UPPER.MARKER.NAMES)) {
		upper.bound <- if (is.na(electrophoresis$peaks$lower.length[upper.marker])) electrophoresis$peaks$length[upper.marker] else electrophoresis$peaks$lower.length[upper.marker]
		result[
			electrophoresis$data$sample.index == electrophoresis$peaks$sample.index[upper.marker] &
			electrophoresis$data$length >= upper.bound
		] <- F
	}
	result
}


#' Scale data by a differential
#'
#' Given an x-variable and a y-variable, this function scales the y-values from the observed data points by the differentials of the x-values. The resulting values of y/dx can then be used to make visually accurate graphs.
#'
#' @param electrophoresis An \code{\link{electrophoresis}} object.
#' @param x The name of the x-variable in \code{electrophoresis$data}.
#' @param y The name of the y-variable in \code{electrophoresis$data}.
#'
#' @return A vector of the y-values, one for each row of \code{electrophoresis$data}, divided by the differentials of the corresponding x-values.
#'
#' @seealso \code{\link{normalize.proportion}}
#'
#' @export
differential.scale <- function(electrophoresis, x, y) {
	stopifnot("sample indexes out of order" = all(diff(electrophoresis$data$sample.index) >= 0))
	delta.x <- unlist(by(electrophoresis$data, electrophoresis$data$sample.index, function(data.subset) c(NA, diff(data.subset[[x]])), simplify = F)) # apply by sample to make sure we don't get a weird delta at the sample boundary
	if (all(delta.x < 0, na.rm = T)) delta.x <- -delta.x else stopifnot("x-values out of order" = all(delta.x > 0, na.rm = T)) # assume data points are monotonic; if negative (like migration distance) make them positive so the math comes out clean
	
	electrophoresis$data[[y]] / delta.x
}

