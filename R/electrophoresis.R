# combine multiple electrophoresis objects
#' @export
rbind.electrophoresis <- function(...) {
	arg.list <- list(...)
	
	# increment the peak indexes in the data so they'll match the new table
	for (i in 1:(length(arg.list) - 1)) for (j in (i + 1):length(arg.list)) arg.list[[j]]$data$peak <- arg.list[[j]]$data$peak + nrow(arg.list[[i]]$peaks)
	
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


# compute the sum of some variable under the electrophoresis curve between some boundaries, or within an annotated peak (or both)
# to do it separately by sample you must provide this function a subset
#' @export
integrate.rawdata <- function(
	data,
	lower.bound = NULL,
	upper.bound = NULL,
	bound.unit = "length",
	peak = NULL,
	sum.unit = "molarity"
) {
	above.lower <- if (is.null(lower.bound)) rep(TRUE, nrow(data)) else data[[bound.unit]] >= lower.bound
	below.upper <- if (is.null(upper.bound)) rep(TRUE, nrow(data)) else data[[bound.unit]] <= upper.bound
	in.peak <- if (is.null(peak)) rep(TRUE, nrow(data)) else ! is.na(data$peak) & data$peak == peak
	sum(data[[sum.unit]][above.lower & below.upper & in.peak])
}

# compute the sum of some variable under the electrophoresis curve within each peak in the peak table
#' @export
integrate.peaks <- function(
	electrophoresis,
	sum.unit = "molarity"
) sapply(1:nrow(electrophoresis$peaks), function(peak) integrate.rawdata(electrophoresis$data, peak = peak, sum.unit = sum.unit))

# compute the sum of some variable under the electrophoresis curve within each region in the region table
#' @export
integrate.regions <- function(
	electrophoresis,
	sum.unit = "molarity"
) sapply(1:nrow(electrophoresis$region), function(region.index) integrate.rawdata(subset(electrophoresis$data, well.number == as.character(electrophoresis$regions$well.number[region.index])), electrophoresis$regions$lower.length[region.index], electrophoresis$regions$upper.length[region.index], sum.unit = sum.unit))

# compute the sum of some variable under the electrophoresis curve within a specified region, for each sample
#' @export
integrate.custom <- function(
	electrophoresis,
	...
) sapply(as.character(electrophoresis$samples$well.number), function(well) integrate.rawdata(subset(electrophoresis$data, well.number == well), ...))

# given two or more regions (pairs of lower and upper bounds), calculate the ratio of the integrated sum of each additional region relative to the integrated sum of the first region, for each sample
# bounds should be a list of pairs of lower and upper bounds, e.g. list(c(100, 200), c(200, 700))
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

# shortcut to do region.ratio for common Illumina library size vs. adapter dimer contamination
#' @export
illumina.library.ratio <- function(
	electrophoresis,
	bounds = list(c(100, 200), c(200, 700)),
	...
) region.ratio(electrophoresis, bounds, ...)

