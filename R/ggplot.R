#' Variable labels for electrophoresis data
#'
#' This function generates a customized descriptive label for a variable in an \code{electrophoresis} object.
#'
#' If your \code{electrophoresis} object contains multiple batches with different units for some reason, this function gives a warning and does not include a unit in the label.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param variable The name of the variable to convert into a label. One of \code{"time"}, \code{"aligned.time"}, \code{"distance"}, \code{"relative.distance"}, \code{"fluorescence"}, \code{"length"}, \code{"concentration"}, \code{"molarity"}.
#'
#' @return A \code{"character"} object containing a descriptive, human-readable label including the correct units (e.g. ng/µl, nM) specified in the original metadata.
#'
#' @export
variable.label <- function(electrophoresis, variable) switch(variable,
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
)

#' Labeller for electrophoresis samples
#'
#' This is a labeller function compatible with \code{\link{facet_wrap}} and \code{\link{facet_grid}}. It allows you to facet the data from an \code{electrophoresis} object on \code{batch * well.number}, to keep the samples in the observed order, but replaces the facet labels with the annotated sample names.
#'
#' @param electrophoresis An \code{electrophoresis} object whose \code{samples} member contains all the samples in your plot.
#'
#' @export
labeller.electrophoresis <- function(electrophoresis) function(factor.frame) list(
	apply(factor.frame, 1, function(labels) {
		which.sample <- which(electrophoresis$samples$batch == labels[1] & electrophoresis$samples$well.number == labels[2])
		stopifnot(length(which.sample) == 1)	
		as.character(electrophoresis$samples$sample.name[which.sample])
	})
)

#' Plot electrophoresis data
#'
#' This function is a shortcut to plot the data from an \code{electrophoresis} object, wrapping \code{\link{ggplot}} similarly to \code{\link{qplot}}. The result is analogous to electropherograms produced by the Agilent software.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param x The variable to use as the x-value of each point in the graph. Can be one of \code{"time"}, \code{"aligned.time"}, \code{"distance"}, \code{"relative.distance"}, or \code{"length"}.
#' @param y The variable to use as the y-value of each point in the graph. Can be one of \code{"fluorescence"}, \code{"concentration"}, or \code{"molarity"}.
#' @param log Which variables to log-transform (\code{"x"}, \code{"y"}, or \code{"xy"}).
#' @param facet If \code{TRUE}, display each sample in a separate facet via \code{\link{facet_wrap}}, as in the Agilent software's usual display mode. If \code{FALSE}, overlay all samples in one graph and color-coded, as in the Agilent software's comparison mode.
#' @param scales Scaling rules for the facets, passed to \code{\link{facet_wrap}}.
#' @param geom Name of the geom to draw. Currently only \code{"line"} (\code{\link{geom_line}}, to get continuous lines) and \code{"area"} (\code{\link{geom_area}}, to fill the area under the curves) are supported.
#' @param include.ladder If \code{FALSE}, graph only the actual samples and not the ladder(s) wells.
#' @param between.markers If \code{TRUE}, graph only data between the marker peaks.
#' @param peak.fill Color to fill the area under reported peaks. Set to \code{NA} to skip plotting the peaks.
#' @param region.alpha Alpha-transparency of the highlight in the reported regions of interested. Set to \code{NA} to skip plotting the regions.
#' @param area.alpha Alpha-transparency of the filled areas under the curves, if they are overlaid in one graph (\code{facet = FALSE} and \code{geom = "area"}), to make them visible through one another.
#'
#' @return A ggplot object containing several layers. You can draw it directly or customize it like any other ggplot object by adding more layers.
#'
#' @export
#' @import ggplot2
qplot.electrophoresis <- function(electrophoresis,
	x = "length",
	y = "molarity",
	log = "",
	facet = TRUE,
	scales = "fixed",
	geom = "line",
	include.ladder = FALSE,
	between.markers = TRUE,
	peak.fill = "darkred",
	region.alpha = 0.2,
	area.alpha = 0.2
) {

	graph.data <- electrophoresis$data
	
	# remove ladders
	if (! include.ladder) graph.data <- subset(graph.data, ! is.ladder)
	
	# remove data in unusable ranges
	graph.data <- graph.data[! is.na(graph.data[[x]]) & ! is.na(graph.data[[y]]),]
	if (log %in% c("x", "xy")) graph.data <- graph.data[graph.data[[x]] > 0,]
	if (log %in% c("y", "xy")) graph.data <- graph.data[graph.data[[y]] > 0,]
	
	# remove data outside the space between markers
	if (between.markers) for (i in 1:nrow(electrophoresis$peaks)) {
		if (electrophoresis$peaks$peak.observations[i] == "Lower Marker") {
			graph.data <- subset(graph.data, ! (
				batch == electrophoresis$peaks$batch[i] &
				well.number == electrophoresis$peaks$well.number[i] &
				length <= electrophoresis$peaks$upper.length[i]
			))
		} else if (electrophoresis$peaks$peak.observations[i] == "Upper Marker") {
			graph.data <- subset(graph.data, ! (
				batch == electrophoresis$peaks$batch[i] &
				well.number == electrophoresis$peaks$well.number[i] &
				length >= electrophoresis$peaks$lower.length[i]
			))
		}
	}
	
	# create plot but don't add the geom yet
	this.plot <- ggplot(graph.data)
	
	# add regions
	if (facet & ! is.na(region.alpha) & ! is.null(electrophoresis$regions)) this.plot <- this.plot + geom_rect(aes_(xmin = as.name(paste0("lower.", x)), xmax = as.name(paste0("upper.", x)), ymin = -Inf, ymax = Inf), data = electrophoresis$regions, alpha = region.alpha)
	
	# finally add the geom (after the regions so it's in front)
	this.plot <- this.plot + switch(geom,
		line = geom_line(if (facet)
			aes_(x = as.name(x), y = as.name(y)) 
		else 
			aes_(x = as.name(x), y = as.name(y), color = as.name("sample.name"))
		),
		area = if (facet)
			geom_area(aes_(x = as.name(x), y = as.name(y))) 
		else
			geom_area(aes_(x = as.name(x), y = as.name(y), fill = as.name("sample.name")), alpha = area.alpha)
	)
	
	# add peaks
	if (facet & ! is.na(peak.fill) & ! is.null(electrophoresis$peaks)) this.plot <- this.plot + geom_area(aes_(x = as.name(x), y = as.name(y), group = as.name("peak")), data = subset(graph.data, ! is.na(peak)), fill = peak.fill)
	
	# add faceting
	if (facet) this.plot <- this.plot + facet_wrap(~ batch * well.number, scales = scales, labeller = labeller.electrophoresis(electrophoresis))
	
	# apply log transformations
	if (log %in% c("x", "xy")) this.plot <- this.plot + scale_x_log10()
	if (log %in% c("y", "xy")) this.plot <- this.plot + scale_y_log10()
	
	# set labels and other settings for specific x & y variables
	this.plot <- this.plot + xlab(variable.label(electrophoresis, x)) + ylab(variable.label(electrophoresis, y))
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
stdcrv.mobility <- function(electrophoresis, n.simulate = 100, line.color = "red") { # returns a ggplot object, which can be extended by adding more features
	ladder.data <- subset(electrophoresis$data, is.ladder & ! is.na(peak))
	ladder.data$true.length <- electrophoresis$peaks$length[ladder.data$peak]
	good.peaks <- subset(electrophoresis$peaks, ! is.na(length))
	
	# determine which kind of data we have
	possible.x.names <- c("aligned.time", "relative.distance")
	x.name <- possible.x.names[possible.x.names %in% colnames(electrophoresis$data)]
	stopifnot(length(x.name) == 1)	
	
	simulated.data <- do.call(rbind, lapply(names(electrophoresis$wells.by.ladder), function(batch) do.call(rbind, lapply(names(electrophoresis$wells.by.ladder[[batch]]), function(ladder.well) {
		x.range <- range(subset(good.peaks, batch == batch & well.number == ladder.well)[[x.name]])
		x.diff <- diff(x.range)
		result <- data.frame(batch, well.number = ladder.well, x = x.range[1] + x.diff * (0:(n.simulate - 1) / (n.simulate - 1)))
		result$estimated.length <- electrophoresis$mobility.functions[[batch]][[ladder.well]](result$x)
		result
	}))))
	this.plot <- ggplot(ladder.data, aes_(x = as.name("true.length"), y = as.name(x.name), color = as.name("fluorescence"))) +
		geom_point() + 
		geom_point(aes_(x = as.name("length"), y = as.name(x.name)), data = subset(good.peaks, is.ladder), color = line.color) + # overlay the reported peak positions
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
			
			peaks <- cbind(electrophoresis$peaks, estimated.variable = NA)
			for (batch in names(electrophoresis$wells.by.ladder)) for (ladder.well in names(electrophoresis$wells.by.ladder[[batch]])) {
				which.peaks <- peaks$well.number %in% electrophoresis$wells.by.ladder[[batch]][[ladder.well]]
				peaks$estimated.variable[which.peaks] <- electrophoresis$mobility.functions[[batch]][[ladder.well]](peaks[[x.name]][which.peaks])
			}
			peaks
		},
		
		concentration = cbind(electrophoresis$peaks, estimated.variable = sapply(1:nrow(electrophoresis$peaks), function(peak.index) sum(electrophoresis$data$concentration[which(electrophoresis$data$peak == peak.index)]))),
		
		molarity = cbind(electrophoresis$peaks, estimated.variable = sapply(1:nrow(electrophoresis$peaks), function(peak.index) sum(electrophoresis$data$molarity[which(electrophoresis$data$peak == peak.index)])))
	)
	
	peaks <- subset(peaks, ! is.na(estimated.variable)) # remove NA's so they don't affect the x-limits and throw a warning
	
	result <- ggplot(peaks, aes_(as.name(variable), as.name("estimated.variable"), color = as.name("peak.observations"))) +
		geom_abline() +
		geom_point() +
		geom_smooth(method = "lm") +
		xlab(paste("software-reported", variable.label(electrophoresis, variable))) +
		ylab(paste("estimated", variable.label(electrophoresis, variable))) +
		facet_wrap(~ batch * well.number, labeller = labeller.electrophoresis(electrophoresis))
	
	if (log) result <- result + scale_x_log10() + scale_y_log10()
	
	result
}

