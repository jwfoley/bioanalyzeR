# compute the sum of some variable under the electrophoresis curve between some boundaries, or within an annotated peak (or both)
# to do it separately by sample you must provide this function a subset
integrate.electrophoresis <- function(
	data,
	lower.bound = NULL,
	upper.bound = NULL,
	peak = NULL,
	bound.unit = "length",
	sum.unit = "delta.molarity"
) {
	above.lower <- if (is.null(lower.bound)) rep(TRUE, nrow(data)) else data[[bound.unit]] >= lower.bound
	below.upper <- if (is.null(upper.bound)) rep(TRUE, nrow(data)) else data[[bound.unit]] <= upper.bound
	in.peak <- if (is.null(peak)) rep(TRUE, nrow(data)) else ! is.na(data$peak) & data$peak == peak
	sum(data[[sum.unit]][above.lower & below.upper & in.peak])
}

