library(ggplot2)

qplot.electrophoresis <- function(data, # returns a ggplot object, which can be extended by adding more features
	x = "length",
	y = "fluorescence",
	log = NULL, # which variables to log-transform, as in qplot(); defaults to "x" if x-value is length or none otherwise
	scales = "fixed", # scaling rules for the facets, passed to facet_wrap()
	geom = geom_line, # another good option is geom_area
	include.ladder = FALSE, # show the ladder wells?
	between.markers = FALSE, # show only data between the marker peaks (so their signal doesn't ruin the scale)
	peak.fill = "darkred", # set to NA to stop showing peaks
	region.alpha = 0.2 # set to NA to stop showing regions
) {
	
	# remove ladders	
	if (! include.ladder) data$data <- subset(data$data, sample.observations != "Ladder")
	
	# remove data outside the space between markers
	if (between.markers) for (i in 1:nrow(data$peaks)) {
		if (data$peaks$peak.observations[i] == "Lower Marker") {
			data$data <- subset(data$data, ! (well.number == data$peaks$well.number[i] & distance >= 0.95 * data$peaks$lower.distance[i]))
		} else if (data$peaks$peak.observations[i] == "Upper Marker") {
			data$data <- subset(data$data, ! (well.number == data$peaks$well.number[i] & distance <= 1.05 * data$peaks$upper.distance[i]))
		}
	}
	
	# create facet labeler
	well.names <- as.character(data$samples$name)
	names(well.names) <- data$samples$well.number
	
	# create plot but don't add the geom yet
	this.plot <- ggplot(data$data) + facet_wrap(~ well.number, scales = scales, labeller = as_labeller(well.names))
	
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
	} else if (x == "relative.distance") {
		this.plot <- this.plot + xlab("distance migrated relative to markers") + scale_x_reverse()
	}
	if (y == "delta.molarity") this.plot <- this.plot + ylab("concentration (pM)")
	
	# add regions
	if (! is.na(region.alpha)) this.plot <- this.plot + geom_rect(aes_(xmin = as.name(paste0("lower.", x)), xmax = as.name(paste0("upper.", x)), ymin = -Inf, ymax = Inf), data = data$regions, alpha = region.alpha)
	
	# finally add the geom (after the regions so it's in front)
	this.plot <- this.plot + geom(aes_(x = as.name(x), y = as.name(y)))
	
	# add peaks
	if (! is.na(peak.fill)) this.plot <- this.plot + geom_area(aes_(x = as.name(x), y = as.name(y), group = as.name("peak")), data = subset(data$data, ! is.na(peak)), fill = peak.fill)
	
	this.plot
}

qc.mobility <- function(data, n.simulate = 100, line.color = "red") { # returns a ggplot object, which can be extended by adding more features
	ladder.data <- subset(data$data, sample.observations == "Ladder" & ! is.na(peak))
	ladder.data$true.length <- data$peaks$length[ladder.data$peak]
	simulated.data <- do.call(rbind, lapply(names(data$wells.by.ladder), function(ladder.well) {
		relative.distance.range <- range(subset(data$peaks, well.number == ladder.well)$relative.distance)
		relative.distance.diff <- diff(relative.distance.range)
		result <- data.frame(well.number = ladder.well, relative.distance = relative.distance.range[1] + relative.distance.diff * (0:(n.simulate - 1) / (n.simulate - 1)))
		result$estimated.length <- data$mobility.functions[[ladder.well]](result$relative.distance)
		result
	}))
	ggplot(ladder.data, aes(x = true.length, y = relative.distance, color = fluorescence)) +
		geom_point() + 
		geom_point(aes(x = length, y = relative.distance), data = subset(data$peaks, sample.observations == "Ladder"), color = line.color) + # overlay the reported peak positions
		geom_line(aes(x = estimated.length, y = relative.distance), data = simulated.data, col = line.color) + # overlay the simulated data from the mobility function +
		scale_y_reverse() +
		xlab("true length (bases)") +
		ylab("distance migrated relative to markers") +
		facet_wrap(~ well.number)
}

qc.molarity <- function(data) {
	peaks <- data$peaks
	peaks$estimated.molarity <- sapply(1:nrow(peaks), function(peak.index) sum(data$data$delta.molarity[which(data$data$peak == peak.index)])) # without the which() you get the NA's too
	
	# create facet labeler
	well.names <- as.character(data$samples$name)
	names(well.names) <- data$samples$well.number
	
	ggplot(peaks, aes(molarity, estimated.molarity)) +
		geom_point() +
		geom_abline() +
		xlab("true molarity (pM)") +
		ylab("estimated molarity (pM)") +
		facet_wrap(~ well.number, labeller = as_labeller(well.names))
}
