library(ggplot2)

qplot.electrophoresis <- function(data, # returns a ggplot object, which can be extended by adding more features
	x = "length",
	y = "fluorescence",
	log = NULL, # which variables to log-transform, as in qplot(); defaults to "x" if x-value is length or none otherwise
	scales = "fixed", # scaling rules for the facets, passed to facet_wrap()
	geom = geom_line,
	include.ladder = FALSE,
	region.alpha = 0.2 # set to NA to stop showing regions
) {
	
	# remove ladders	
	if (! include.ladder) data$data <- subset(data$data, sample.observations != "Ladder")
	
	# create plot but don't add the geom yet (it should be added last to be the front layer)
	this.plot <- ggplot(data$data) + facet_wrap(~ well.number, scales = scales)
	
	# apply log transformations
	if (
		(is.null(log) && x == "length") ||
		(! is.null(log) && log %in% c("x", "xy"))
	) this.plot <- this.plot + scale_x_log10()
	if (! is.null(log) && log %in% c("y", "xy")) this.plot <- this.plot + scale_y_log10()
	
	# set labels and other settings for specific x & y variables
	if (x == "length") {
		this.plot <- this.plot + xlab("length (bases)")
	} else if (x == "time") {
		this.plot <- this.plot + xlab("time (s)")
	} else if (x == "distance") {
		this.plot <- this.plot + xlab("distance migrated") + scale_x_reverse()
	}
	if (y == "delta.molarity") this.plot <- this.plot + ylab("concentration (pM)")
	
	# add regions
	if (! is.na(region.alpha)) {
		this.plot <- this.plot + geom_rect(aes(xmin = lower.bound, xmax = upper.bound, ymin = -Inf, ymax = Inf), data = data$regions, alpha = region.alpha)
	}
	
	# finally add the geom
	this.plot <- this.plot + geom(aes_(x = as.name(x), y = as.name(y)))
	
	this.plot
}

