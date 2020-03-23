#' Variable labels for electrophoresis data
#'
#' This function generates a customized descriptive label for a variable in an \code{electrophoresis} object.
#'
#' If your \code{electrophoresis} object contains multiple batches with different units for some reason, this function gives a warning and does not include a unit in the label.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param variable The name of the variable to convert into a label. One of \code{"time"}, \code{"aligned.time"}, \code{"distance"}, \code{"relative.distance"}, \code{"fluorescence"}, \code{"length"}, \code{"concentration"}, \code{"molarity"}.
#'
#' @param variable2 Optionally, the name of a second variable that is the denominator of the first (e.g. molarity per length).
#'
#' @return A \code{"character"} object containing a descriptive, human-readable label including the correct units (e.g. ng/µl, nM) specified in the original metadata.
#'
#' @export
variable.label <- function(electrophoresis, variable, variable2 = NULL) if (is.null(variable2)) switch(variable,
	time =               "time (s)",
	aligned.time =       "aligned time relative to markers (s)",
	distance =           "distance migrated",
	relative.distance =  "distance migrated relative to markers",
	fluorescence =       "fluorescence",
	length = {
		length.units <- unique(sapply(electrophoresis$assay.info, function(x) x$length.unit))
		if (length(length.units) == 1) paste0("length (", length.units, ")") else {
			warning("incompatible length units")
			"length"
		}
	},
	concentration = {
		concentration.units <- unique(sapply(electrophoresis$assay.info, function(x) x$concentration.unit))
		if (length(concentration.units) == 1) paste0("concentration (", concentration.units, ")") else {
			warning("incompatible concentration units")
			"concentration"
		}
	},
	molarity = {
		molarity.units <- unique(sapply(electrophoresis$assay.info, function(x) x$molarity.unit))
		if (length(molarity.units) == 1) paste0("molarity (", molarity.units, ")") else {
			warning("incompatible molarity units")
			"molarity"
		}
	},	
) else paste(variable.label(electrophoresis, variable), "per", variable.label(electrophoresis, variable2))

#' Labeller for electrophoresis samples
#'
#' This is a labeller function compatible with \code{\link{facet_wrap}} and \code{\link{facet_grid}}. It allows you to facet the data from an \code{electrophoresis} object on \code{sample.index}, to keep the samples in the observed order, but replaces the facet labels with the annotated sample names.
#'
#' @param electrophoresis An \code{electrophoresis} object whose \code{samples} member contains all the samples in your plot.
#'
#' @export
labeller.electrophoresis <- function(electrophoresis) function(factor.frame) {
	stopifnot(ncol(factor.frame) == 1)
	stopifnot(class(factor.frame[,1]) == "integer")
	list(as.character(electrophoresis$samples$sample.name[factor.frame[,1]]))
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
#' @export
scale.by.differential <- function(electrophoresis, x, y) {
	stopifnot(all(diff(electrophoresis$data$sample.index) %in% c(0, 1))) # assume data points from each sample are contiguous and ordered by sample
	delta.x <- do.call(c, lapply(unique(electrophoresis$data$sample.index), function(i) c(NA, diff(electrophoresis$data[electrophoresis$data$sample.index == i,x])))) # apply by sample to make sure we don't get a weird delta at the sample boundary
	if (all(delta.x < 0, na.rm = T)) delta.x <- -delta.x else stopifnot(all(delta.x > 0, na.rm)) # assume data points are monotonic; if negative (like migration distance) make them positive so the math comes out clean
	
	electrophoresis$data[[y]] / delta.x
}

#' Plot electrophoresis data
#'
#' This function is a shortcut to plot the data from an \code{electrophoresis} object, wrapping \code{\link{ggplot}} similarly to \code{\link{qplot}}. The result is analogous to electropherograms produced by the Agilent software.
#'
#' Before plotting, the y-variable is scaled by the differentials in the x-value. Thus the units of the y-axis are divided by the units of the x-axis, e.g. molarity per length. This ensures that the area under the curve (width times height) represents the desired variable in the correct units. For example, if the x-variable is length in bp, the graph will be equivalent to a histogram with one bar for each possible molecule length in base pairs.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param x The variable to use as the x-value of each point in the graph. Can be one of \code{"time"}, \code{"aligned.time"}, \code{"distance"}, \code{"relative.distance"}, or \code{"length"}.
#' @param y The variable to use as the y-value of each point in the graph. Can be one of \code{"fluorescence"}, \code{"concentration"}, or \code{"molarity"}.
#' @param log Which variables to log-transform (\code{"x"}, \code{"y"}, or \code{"xy"}).
#' @param facets Faceting formula to use. Picks \code{\link{facet_wrap}} or \code{\link{facet_grid}} depending on whether the formula is one- or two-sided. If \code{NULL}, overlay all samples in one color-coded graph.
#' @param margins Display marginal facets (via \code{\link{facet_grid}}) if using a two-side faceting formula.
#' @param scales Scaling rules for the facets, passed to \code{\link{facet_wrap}}.
#' @param geom Name of the geom to draw. Currently only \code{"line"} (\code{\link{geom_line}}, to get continuous lines) and \code{"area"} (\code{\link{geom_area}}, to fill the area under the curves) are supported.
#' @param include.ladder If \code{FALSE}, graph only the actual samples and not the ladder(s) wells.
#' @param between.markers If \code{TRUE}, graph only data between the marker peaks.
#' @param xlim, ylim Limits of x- and y-axes 
#' @param peak.fill Color to fill the area under reported peaks. Set to \code{NA} to skip plotting the peaks.
#' @param region.alpha Alpha-transparency of the highlight in the reported regions of interested. Set to \code{NA} to skip plotting the regions.
#' @param area.alpha Alpha-transparency of the filled areas under the curves, if they are overlaid in one graph (\code{facet = FALSE} and \code{geom = "area"}), to make them visible through one another.
#' @param title, xlab, ylab Plot title, x-axis label, and y-axis label.
#'
#' @return A ggplot object containing several layers. You can draw it directly or customize it like any other ggplot object by adding more layers.
#'
#' @export
#' @import ggplot2
qplot.electrophoresis <- function(electrophoresis,
	x = "length",
	y = "molarity",
	log = "",
	facets = ~ sample.index,
	margins = FALSE,
	scales = "fixed",
	geom = "line",
	include.ladder = FALSE,
	between.markers = TRUE,
	xlim = c(NA, NA),
	ylim = c(NA, NA),
	peak.fill = "darkred",
	region.alpha = 0.2,
	area.alpha = 0.2,
	title = NULL,
	xlab = NULL,
	ylab = NULL
) {

	# remove ladders
	if (! include.ladder) electrophoresis <- subset(electrophoresis, well.number != ladder.well)
	
	# remove data in unusable ranges
	electrophoresis$data <- electrophoresis$data[! is.na(electrophoresis$data[[x]]) & ! is.na(electrophoresis$data[[y]]),]
	if (log %in% c("x", "xy")) electrophoresis$data <- electrophoresis$data[electrophoresis$data[[x]] > 0,]
	if (log %in% c("y", "xy")) electrophoresis$data <- electrophoresis$data[electrophoresis$data[[y]] > 0,]
	
	# remove data outside the space between markers
	if (between.markers) for (i in 1:nrow(electrophoresis$peaks)) {
		if (electrophoresis$peaks$peak.observations[i] == "Lower Marker") {
			electrophoresis$data <- subset(electrophoresis$data, ! (
				sample.index == electrophoresis$peaks$sample.index[i] &
				length <= electrophoresis$peaks$upper.length[i]
			))
		} else if (electrophoresis$peaks$peak.observations[i] == "Upper Marker") {
			electrophoresis$data <- subset(electrophoresis$data, ! (
				sample.index == electrophoresis$peaks$sample.index[i] &
				length >= electrophoresis$peaks$lower.length[i]
			))
		}
	}
	
	# annotate data with metadata
	electrophoresis$data <- cbind(electrophoresis$data, electrophoresis$samples[electrophoresis$data$sample.index,])
	if (! is.null(electrophoresis$peaks)) electrophoresis$peaks <- cbind(electrophoresis$peaks, electrophoresis$samples[electrophoresis$peaks$sample.index,])
	if (! is.null(electrophoresis$regions)) electrophoresis$regions <- cbind(electrophoresis$regions, electrophoresis$samples[electrophoresis$regions$sample.index,])
	
	# scale y-values to dx
	electrophoresis$data$y.scaled <- scale.by.differential(electrophoresis, x, y)
	
	# create plot but don't add the geom yet
	this.plot <- ggplot(electrophoresis$data)
	
	# add regions
	if (! is.null(facets) & ! is.na(region.alpha) & ! is.null(electrophoresis$regions)) this.plot <- this.plot + geom_rect(aes_(xmin = as.name(paste0("lower.", x)), xmax = as.name(paste0("upper.", x)), ymin = -Inf, ymax = Inf), data = electrophoresis$regions, alpha = region.alpha)
	
	# finally add the geom (after the regions so it's in front)
	this.plot <- this.plot + switch(geom,
		line = geom_line(if (! is.null(facets))
			aes_(x = as.name(x), y = as.name("y.scaled"))
		else 
			aes_(x = as.name(x), y = as.name("y.scaled"), color = as.name("sample.name"))
		),
		area = if (! is.null(facets))
			geom_area(aes_(x = as.name(x), y = as.name("y.scaled")))
		else
			geom_area(aes_(x = as.name(x), y = as.name("y.scaled"), fill = as.name("sample.name")), alpha = area.alpha)
	)
	
	# add peaks
	if (! is.null(facets) & ! is.na(peak.fill) & ! is.null(electrophoresis$peaks)) {
		peak.data <- subset(cbind(electrophoresis$data, peak = in.peaks(electrophoresis)), ! is.na(peak))
		this.plot <- this.plot + geom_area(aes_(x = as.name(x), y = as.name("y.scaled"), group = as.name("peak")), data = peak.data, fill = peak.fill)
	}
	
	# add faceting
	if (! is.null(facets)) this.plot <- this.plot +
		if (length(facets) == 2) # one-sided formula
			facet_wrap(facets, scales = scales, labeller = if (facets == ~ sample.index) labeller.electrophoresis(electrophoresis) else "label_value")
		else
			facet_grid(facets, margins = margins, scales = scales)
	
	# apply limits
	if (! all(is.na(xlim))) this.plot <- this.plot + lims(x = xlim)
	if (! all(is.na(ylim))) this.plot <- this.plot + lims(y = ylim)
	
	# apply log transformations
	if (log %in% c("x", "xy")) this.plot <- this.plot + scale_x_log10()
	if (log %in% c("y", "xy")) this.plot <- this.plot + scale_y_log10()
	
	# set labels and other settings for specific x & y variables
	this.plot <- this.plot + labs(
		x = if (! is.null(xlab)) xlab else variable.label(electrophoresis, x),
		y = if (! is.null(ylab)) ylab else variable.label(electrophoresis, y, x),
		title = title
	)
	if (x %in% c("distance", "relative.distance")) this.plot <- this.plot + scale_x_reverse()
	
	this.plot
}

#' Plot mobility standard curves
#'
#' This function is a shortcut to plot the standard curve(s) of molecule length vs. migration speed from an \code{electrophoresis} object, wrapping \code{\link{ggplot}}. This allows you to check the quality of the model.
#'
#' The positions of the ladder peaks reported by the Agilent software are shown in the selected color, and the fluorescence intensites within the peak boundaries are also plotted. If there are multiple ladders, each is shown as a separate facet.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param n.simulate Number of data points to simulate for drawing the standard curve.
#' @param line.color Color of the standard curve and data points.
#'
#' @return A ggplot object containing several layers. You can draw it directly or customize it like any other ggplot object by adding more layers.
#' 
#' @seealso \code{\link{qc.electrophoresis}}
#'
#' @import ggplot2
#'
#' @export
stdcrv.mobility <- function(electrophoresis, n.simulate = 100, line.color = "red") { # returns a ggplot object, which can be extended by adding more features
	which.ladders <- which(electrophoresis$samples$well.number == electrophoresis$samples$ladder.well)
	ladder.data <- subset(cbind(electrophoresis$data, peak = in.peaks(electrophoresis)), sample.index %in% which.ladders & ! is.na(peak))
	ladder.data$true.length <- electrophoresis$peaks$length[ladder.data$peak]
	ladder.peaks <- subset(electrophoresis$peaks, sample.index %in% which.ladders & ! is.na(length))
	
	# determine which kind of data we have
	possible.x.names <- c("aligned.time", "relative.distance")
	x.name <- possible.x.names[possible.x.names %in% colnames(electrophoresis$data)]
	stopifnot(length(x.name) == 1)	
	
	simulated.data <- do.call(rbind, lapply(which.ladders, function(ladder.index) {
		x.range <- range(ladder.peaks[[x.name]][ladder.peaks$sample.index == ladder.index])
		x.diff <- diff(x.range)
		batch <- as.character(electrophoresis$samples$batch[ladder.index])
		well.number <- as.character(electrophoresis$samples$well.number[ladder.index])
		result <- data.frame(batch, well.number, x = x.range[1] + x.diff * (0:(n.simulate - 1) / (n.simulate - 1)))
		result$estimated.length <- electrophoresis$mobility.functions[[batch]][[well.number]](result$x)
		result
	}))
	this.plot <- ggplot(ladder.data, aes_(x = as.name("true.length"), y = as.name(x.name), color = as.name("fluorescence"))) +
		geom_point() + 
		geom_point(aes_(x = as.name("length"), y = as.name(x.name)), data = ladder.peaks, color = line.color) + # overlay the reported peak positions
		geom_line(aes(x = estimated.length, y = x), data = simulated.data, col = line.color) + # overlay the simulated data from the mobility function
		xlab(paste("true", variable.label(electrophoresis, "length"))) +
		ylab(variable.label(electrophoresis, x.name)) +
		facet_wrap(~ batch * well.number)
	if (x.name == "relative.distance") this.plot <- this.plot + scale_y_reverse()
	
	this.plot
}

#' Compare estimates with Agilent software's output
#'
#' This function is a shortcut to plot the estimates of a variable from the \code{\link{bioanalyzeR}} package against the values reported by the Agilent software, at all reported peaks. This allows you to check the quality of the estimates.
#'
#' If \code{variable = "length"}, apply the mobility model(s) fit by \code{\link{read.bioanalyzer}} or \code{\link{read.tapestation}} to the peak centers' aligned times or relative distance and compare those estimated molecule lengths with the lengths reported by the Agilent software. If \code{variable = "concentration"} or \code{variable = "molarity"}, integrate the appropriate variable between the boundaries of each peak and compare that sum with the figure reported by the Agilent software. For the ladder Agilent's reported values are hardcoded to the known, correct properties of the ladder molecules.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param variable Which variable to compare. One of \code{"length"}, \code{"concentration"}, \code{"molarity"}.
#' @param log Whether to log-scale both axes. Typically the data are more evenly spaced in log scale.
#'
#' @return A ggplot object containing several layers. You can draw it directly or customize it like any other ggplot object by adding more layers.
#' 
#' @seealso \code{\link{stdcrv.mobility}}
#'
#' @export
#' @import ggplot2
qc.electrophoresis <- function(electrophoresis, variable, log = TRUE) {
	peaks <- switch(variable,
		length = {
			possible.x.names <- c("aligned.time", "relative.distance")
			x.name <- possible.x.names[possible.x.names %in% colnames(electrophoresis$data)]
			stopifnot(length(x.name) == 1)	
			
			result <- cbind(electrophoresis$peaks, estimated.variable = NA)
			for (i in 1:nrow(electrophoresis$samples)) {
				which.peaks <- result$sample.index == i
				result$estimated.variable[which.peaks] <- electrophoresis$mobility.functions[[as.character(electrophoresis$samples$batch[i])]][[as.character(electrophoresis$samples$ladder.well[i])]](result[[x.name]][which.peaks])
			}
			result
		},
		
		concentration = cbind(electrophoresis$peaks, estimated.variable = integrate.peaks(electrophoresis, "concentration")),
		
		molarity = cbind(electrophoresis$peaks, estimated.variable = integrate.peaks(electrophoresis, "molarity"))
	)
	
	peaks <- subset(peaks, ! is.na(estimated.variable)) # remove NA's so they don't affect the x-limits and throw a warning
	
	result <- ggplot(peaks, aes_(as.name(variable), as.name("estimated.variable"), color = as.name("peak.observations"))) +
		geom_abline() +
		geom_point() +
		geom_smooth(method = "lm") +
		xlab(paste("software-reported", variable.label(electrophoresis, variable))) +
		ylab(paste("estimated", variable.label(electrophoresis, variable))) +
		facet_wrap(~ sample.index, labeller = labeller.electrophoresis(electrophoresis))
	
	if (log) result <- result + scale_x_log10() + scale_y_log10()
	
	result
}

