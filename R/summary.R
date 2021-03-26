#' Summarize lengths in a subset of data
#'
#' Compute summary statistics of the molecule length distribution in a subset of an \code{\link{electrophoresis}} object. This helper function is called by other functions that split the object into subsets.
#'
#' @param sample.frame A subset of an \code{\link{electrophoresis}} object containing only contiguous data from one sample.
#'
#' @seealso \code{\link{summarize.peak.region}}, \code{\link{summarize.custom}}
#'
#' @export
summarize.subset <- function(sample.frame) {
	total.molarity <- sum(sample.frame$molarity)
	sample.median <- round(sample.frame$length[min(which(cumsum(sample.frame$molarity) >= total.molarity / 2))])
	sample.mean <- sum(sample.frame$molarity * sample.frame$length) / total.molarity
	length.residuals <- sample.frame$length - sample.mean
	sample.sd <- sqrt(sum(sample.frame$molarity * length.residuals^2) / total.molarity)
	sample.skewness <- sum(sample.frame$molarity * length.residuals^3) / total.molarity / sample.sd^3
	sample.kurtosis <- sum(sample.frame$molarity * length.residuals^4) / total.molarity / sample.sd^4
	
	c(
		Median = sample.median,
		Mean = sample.mean,
		SD = sample.sd,
		Skewness = sample.skewness,
		Kurtosis = sample.kurtosis
	)
}


#' Summarize lengths in peaks or regions
#'
#' Compute summary statistics of the molecule length distribution in the reported peaks or regions in an \code{\link{electrophoresis}} object.
#'
#' @param electrophoresis An \code{\link{electrophoresis}} object.
#' @param index The index, or a vector of indexes, of the peaks or regions to summarize (row numbers in \code{electrophoresis$peaks} or \code{electrophoresis$regions}).
#'
#' @seealso \code{\link{summarize.custom}}
#'
#' @name summarize.peak.region
NULL


#' @rdname summarize.peak.region
#' @export
summarize.peak <- function(
	electrophoresis,
	index = seq(nrow(electrophoresis$peaks))
) as.data.frame(t(sapply(index, function(i) summarize.subset(electrophoresis$data[in.peak(electrophoresis, i),]))))


#' @rdname summarize.peak.region
#' @export
summarize.region <- function(
	electrophoresis,
	index = seq(nrow(electrophoresis$regions))
) as.data.frame(t(sapply(index, function(i) summarize.subset(electrophoresis$data[in.region(electrophoresis, i),]))))


#' Summarize lengths in a custom region
#'
#' Compute summary statistics of the molecule length distribution between specified boundaries. The summary is computed individually for each sample.
#'
#' @param electrophoresis An \code{\link{electrophoresis}} object.
#' @param lower.bound Lower boundary of the region to summarize.
#' @param upper.bound Upper boundary of the region to summarize.
#'
#' @seealso \code{\link{summarize.peak}}, \code{\link{summarize.region}}
#'
#' @export
summarize.custom <- function(
	electrophoresis,
	lower.bound = -Inf,
	upper.bound = Inf
) {
	stopifnot("upper bound must be greater than lower bound" = upper.bound > lower.bound)
	in.this.region <- in.custom.region(electrophoresis$data, lower.bound, upper.bound, "length")
	result <- as.data.frame(t(simplify2array(by(electrophoresis$data[in.this.region,], electrophoresis$data$sample.index[in.this.region], summarize.subset))))
	
	if (lower.bound == -Inf) {
		if (upper.bound != Inf) { # bounded only on right
			colnames(result) <- paste(colnames(result), "below", upper.bound)
		}
	} else if (upper.bound == Inf) { # bounded only on left
		colnames(result) <- paste(colnames(result), "above", lower.bound)
	} else { # bounded on both sides
		colnames(result) <- paste0(colnames(result), " in ", lower.bound, "-", upper.bound)
	}
	
	result
}


