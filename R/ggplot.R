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
	variable # if not found just return the original variable name	
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

#' Plot electrophoresis data
#'
#' This function is a shortcut to plot the data from an \code{electrophoresis} object, wrapping \code{\link{ggplot}} similarly to \code{\link{qplot}}. The result is analogous to electropherograms produced by the Agilent software.
#'
#' @section Variable transformation:
#' Before plotting, unless the y-variable is fluorescence, it is scaled by the differentials in the x-value. Thus the units of the y-axis are divided by the units of the x-axis, e.g. molarity per length. This ensures that the area under the curve (width times height) represents the desired variable in the correct units. For example, if the x-variable is length in bp, the graph will be equivalent to a histogram with one bar for each possible molecule length in base pairs.
#'
#' @section Geoms and aesthetics:
#' The x- and y-variable names must be provided as \code{\link{character}} objects, e.g. \code{"length"} rather than \code{length}. However, additional aesthetics in \code{...} are passed directly to \code{\link[ggplot2]{aes}} so they must be provided as the variable names themselves, e.g. \code{color = sample.observations} rather than \code{color = "sample.observations"}.
#'
#' Both of the supported geoms for the main data, \code{\link[ggplot2]{geom_line}} and \code{\link[ggplot2]{geom_area}}, are given the specified x-variable as \code{x} and the appropriately transformed y-variable as \code{y}, as well as \code{group = sample.index} to ensure overlaid samples are plotted distinctly. If \code{facets} is null, then the geom is also given either \code{color = sample.name} (if \code{geom == "line"}) or \code{fill = sample.name} (if \code{geom == "area"}) to color-code overlaid samples. Thus if replicates of the same sample have the same name, they will also have the same color.
#'
#' If \code{peak.fill} is not NA and \code{facets} is not null and \code{electrophoresis$peaks} is not null, then the plot gets an extra \code{\link[ggplot2]{geom_area}} from only the subset of \code{electrophoresis$data} that is within the reported peaks, using the same settings as above except \code{fill = peak.fill}. If you use a faceting formula other than \code{NULL} that allows some samples to be overlaid in the same facet, it may be a good idea to set \code{peak.fill = NA} so their peaks don't overlap.
#'
#' If \code{region.alpha} is not NA and \code{facets} is not null and \code{electrophoresis$regions} is not null, then the plot gets a \code{\link[ggplot2]{geom_rect}} with \code{ymin = -Inf, ymax = Inf} and \code{xmin} and \code{xmax} set to the lower and upper boundaries of the regions, while \code{alpha = region.alpha}.
#'
#' @param electrophoresis An \code{electrophoresis} object.
#' @param x The name of the variable to use as the x-value of each point in the graph as a character vector. Usually one of \code{"time"}, \code{"aligned.time"}, \code{"distance"}, \code{"relative.distance"}, or \code{"length"}.
#' @param y The name of the variable to use as the y-value of each point in the graph, as a character vector. Usually one of \code{"fluorescence"}, \code{"concentration"}, or \code{"molarity"}.
#' @param ... Additional aesthetics passed to the geom for the main data (not the peaks or regions).
#' @param log Which variables to log-transform (\code{"x"}, \code{"y"}, or \code{"xy"}).
#' @param normalize Normalize the y-value in each sample with \code{\link{normalize.proportion}}.
#' @param facets Faceting formula to use. Picks \code{\link[ggplot2]{facet_wrap}} or \code{\link[ggplot2]{facet_grid}} depending on whether the formula is one- or two-sided. If \code{NULL}, overlay all samples in one color-coded graph.
#' @param margins Display marginal facets (via \code{\link[ggplot2]{facet_grid}}) if using a two-side faceting formula.
#' @param scales Scaling rules for the facets, passed to \code{\link[ggplot2]{facet_wrap}}.
#' @param geom Name of the geom to draw. Only \code{"line"} (\code{\link[ggplot2]{geom_line}}, to get continuous curves) and \code{"area"} (\code{\link[ggplot2]{geom_area}}, to fill the area under the curves) are supported.
#' @param include.ladder If \code{FALSE}, graph only the actual samples and not the ladder well(s).
#' @param include.markers If \code{FALSE}, graph only data between the marker peaks.
#' @param lower.marker.spread If normalizing the totals or excluding marker peaks, extend the lower marker peak by this amount (via \code{\link{between.markers}}).
#' @param xlim,ylim Limits of x- and y-axes 
#' @param peak.fill Color to fill the area under reported peaks. Set to \code{NA} to skip plotting the peaks.
#' @param region.alpha Alpha-transparency of the highlight in the reported regions of interested. Set to \code{NA} to skip plotting the regions.
#' @param area.alpha Alpha-transparency of the filled areas under the curves, if they are overlaid in one graph (\code{facet = FALSE} and \code{geom = "area"}), to make them visible through one another.
#' @param title,xlab,ylab Plot title, x-axis label, and y-axis label.
#'
#' @return A ggplot object containing several layers. You can draw it directly or customize it like any other ggplot object by adding more layers.
#'
#' @seealso \code{\link{sparkline.electrophoresis}}
#'
#' @export
#' @import ggplot2
qplot.electrophoresis <- function(electrophoresis,
	x = "length",
	y = "molarity",
	...,
	log = "",
	normalize = FALSE,
	facets = ~ sample.index,
	margins = FALSE,
	scales = "fixed",
	geom = "line",
	include.ladder = FALSE,
	include.markers = FALSE,
	lower.marker.spread = 5,
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
	
	# remove data outside the space between markers
	if (! include.markers) electrophoresis$data <- electrophoresis$data[which(between.markers(electrophoresis, lower.marker.spread)),]
	
	# remove data in unusable ranges
	electrophoresis$data <- electrophoresis$data[! is.na(electrophoresis$data[[x]]) & ! is.na(electrophoresis$data[[y]]),]
	if (log %in% c("x", "xy")) electrophoresis$data <- electrophoresis$data[electrophoresis$data[[x]] > 0,]
	if (log %in% c("y", "xy")) electrophoresis$data <- electrophoresis$data[electrophoresis$data[[y]] > 0,]
	
	# annotate data with metadata
	electrophoresis$data <- cbind(electrophoresis$data, electrophoresis$samples[electrophoresis$data$sample.index,])
	if (! is.null(electrophoresis$peaks)) electrophoresis$peaks <- cbind(electrophoresis$peaks, electrophoresis$samples[electrophoresis$peaks$sample.index,])
	if (! is.null(electrophoresis$regions)) electrophoresis$regions <- cbind(electrophoresis$regions, electrophoresis$samples[electrophoresis$regions$sample.index,])
	
	# normalize and scale y-values
	electrophoresis$data$y.normalized <- if (normalize) normalize.proportion(electrophoresis, y, lower.marker.spread) else electrophoresis$data[[y]]
	electrophoresis$data$y.scaled <- if (y == "fluorescence") electrophoresis$data$y.normalized else scale.by.differential(electrophoresis, x, "y.normalized") # don't scale fluorescence by differentials
	
	# remove data outside x-limits (after normalization so that's not distorted)
	electrophoresis$data <- electrophoresis$data[which(
		(is.na(xlim[1]) | electrophoresis$data[[x]] >= xlim[1]) &
		(is.na(xlim[2]) | electrophoresis$data[[x]] <= xlim[2])
	),]
	
	# also rename x-variable so aesthetics are easy
	electrophoresis$data$x.value <- electrophoresis$data[[x]]
	
	# create plot but don't add the geom yet
	this.plot <- ggplot(electrophoresis$data)
	
	# add regions
	if (! is.null(facets) & ! is.na(region.alpha) & ! is.null(electrophoresis$regions)) this.plot <- this.plot + geom_rect(aes_(xmin = as.name(paste0("lower.", x)), xmax = as.name(paste0("upper.", x)), ymin = -Inf, ymax = Inf), data = electrophoresis$regions, alpha = region.alpha)
	
	# finally add the geom (after the regions so it's in front)
	this.plot <- this.plot + switch(geom,
		line = geom_line(if (! is.null(facets))
			aes(x = x.value, y = y.scaled, group = sample.index, ...)
		else 
			aes(x = x.value, y = y.scaled, group = sample.index, color = sample.name, ...)
		),
		area = if (! is.null(facets))
			geom_area(aes(x = x.value, y = y.scaled, group = sample.index, ...))
		else
			geom_area(aes(x = x.value, y = y.scaled, group = sample.index, fill = sample.name, ...), alpha = area.alpha)
	)
	
	# add peaks
	if (! is.null(facets) & ! is.na(peak.fill) & ! is.null(electrophoresis$peaks)) {
		peak.data <- subset(cbind(electrophoresis$data, peak = in.peaks(electrophoresis)), ! is.na(peak))
		this.plot <- this.plot + geom_area(aes(x = x.value, y = y.scaled, group = peak), data = peak.data, fill = peak.fill)
	}
	
	# add faceting
	if (! is.null(facets)) this.plot <- this.plot +
		if (length(facets) == 2) # one-sided formula
			facet_wrap(facets, scales = scales, labeller = if (facets == ~ sample.index) labeller.electrophoresis(electrophoresis) else "label_value")
		else
			facet_grid(facets, margins = margins, scales = scales, labeller = if (facets == (. ~ sample.index) || facets == (sample.index ~ .)) labeller.electrophoresis(electrophoresis) else "label_value")
	
	# apply limits
	if (! all(is.na(xlim))) this.plot <- this.plot + lims(x = xlim)
	if (! all(is.na(ylim))) this.plot <- this.plot + lims(y = ylim)
	
	# apply log transformations
	if (log %in% c("x", "xy")) this.plot <- this.plot + scale_x_log10()
	if (log %in% c("y", "xy")) this.plot <- this.plot + scale_y_log10()
	
	# set labels and other settings for specific x & y variables
	this.plot <- this.plot + labs(
		x = if (! is.null(xlab)) xlab else variable.label(electrophoresis, x),
		y = if (! is.null(ylab)) ylab else variable.label(electrophoresis, (if (normalize) paste("proportion of", y) else y), if (y == "fluorescence") NULL else x),
		title = title
	)
	if (x %in% c("distance", "relative.distance")) this.plot <- this.plot + scale_x_reverse()
	
	this.plot
}

#' Sparklines of electrophoresis data
#'
#' This function is a shortcut that calls \code{\link{qplot.electrophoresis}} with customized settings to generate sparklines.
#'
#' In addition to the hardcoded default arguments for \code{\link{qplot.electrophoresis}}, the following settings are added with \code{\link[ggplot2]{theme}}: \preformatted{
#' axis.text.y = element_blank(),
#' axis.ticks.y = element_blank(),
#' panel.grid = element_blank(),
#' panel.background = element_blank(),
#' strip.background = element_blank(),
#' strip.text.y = element_text(angle = 0)
#'}
#'
#' @param ... Arguments passed to \code{\link{qplot.electrophoresis}}.
#' @param facets Passed to \code{\link{qplot.electrophoresis}} but the default is hardcoded for one column of sparklines.
#' @param scales,geom,peak.fill Passed to \code{\link{qplot.electrophoresis}} but the defaults are hardcoded for sparklines and probably do not make sense to change.
#'
#' @references
#' Tufte, Edward R. (1983) The Visual Display of Quantitative Information. Cheshire, Conn.: Graphics Press.
#'
#' @seealso \code{\link{qplot.electrophoresis}}
#'
#' @import ggplot2
#' @export
sparkline.electrophoresis <- function(
	...,
	facets = sample.index ~ .,
	scales = "free_y",
	geom = "line",
	peak.fill = NA
) qplot.electrophoresis(..., facets = facets, scales = scales, geom = geom, peak.fill = peak.fill) + theme(
	axis.text.y = element_blank(),
	axis.ticks.y = element_blank(),
	panel.grid = element_blank(),
	panel.background = element_blank(),
	strip.background = element_blank(),
	strip.text.y = element_text(angle = 0)
)

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
	x.name <- get.x.name(electrophoresis)
	which.ladders <- which(electrophoresis$samples$well.number == electrophoresis$samples$ladder.well)
	ladder.data <- subset(cbind(electrophoresis$data, peak = in.peaks(electrophoresis)), sample.index %in% which.ladders & ! is.na(peak))
	ladder.data$true.length <- electrophoresis$peaks$length[ladder.data$peak]
	ladder.peaks <- subset(electrophoresis$peaks, sample.index %in% which.ladders & ! is.na(length))
	
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
			x.name <- get.x.name(electrophoresis)			
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
		geom_smooth(method = "lm", formula = y ~ x) +
		xlab(paste("software-reported", variable.label(electrophoresis, variable))) +
		ylab(paste("estimated", variable.label(electrophoresis, variable))) +
		facet_wrap(~ sample.index, labeller = labeller.electrophoresis(electrophoresis))
	
	if (log) result <- result + scale_x_log10() + scale_y_log10()
	
	result
}

