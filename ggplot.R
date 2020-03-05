library(ggplot2)

# ggplot2 labeller function to take batch and well.number factors (expecting ~ batch * well.number) and return facet labels that are the sample names
# batch names are not included in the labels
labeller.electrophoresis <- function(electrophoresis) function(factor.frame) list(
	apply(factor.frame, 1, function(labels) {
		which.sample <- which(electrophoresis$samples$batch == labels[1] & electrophoresis$samples$well.number == labels[2])
		stopifnot(length(which.sample) == 1)	
		as.character(electrophoresis$samples$sample.name[which.sample])
	})
)

qplot.electrophoresis <- function(electrophoresis, # returns a ggplot object, which can be extended by adding more features
	x = "length",
	y = "fluorescence",
	log = NULL, # which variables to log-transform, as in qplot(); defaults to "x" if x-value is length or none otherwise
	facet = TRUE, # if false, overlay everything on one graph and color-code the curves
	scales = "fixed", # scaling rules for the facets, passed to facet_wrap()
	geom = "line", # default gives geom_line; "area" for geom_area is another good option
	include.ladder = FALSE, # show the ladder wells?
	between.markers = FALSE, # show only data between the marker peaks (so their signal doesn't ruin the scale)
	peak.fill = "darkred", # set to NA to stop showing peaks
	region.alpha = 0.2, # set to NA to stop showing regions
	area.alpha = 0.2 # if facet = FALSE and geom = "area", alpha transparency of the filled areas
) {

	if (is.null(log)) log <- if (x == "length") "x" else NA
	
	# remove ladders	
	if (! include.ladder) electrophoresis$data <- subset(electrophoresis$data, sample.observations != "Ladder")
	
	# remove data in unusable ranges
	electrophoresis$data <- electrophoresis$data[! is.na(electrophoresis$data[[x]]) & ! is.na(electrophoresis$data[[y]]),]
	if (log %in% c("x", "xy")) electrophoresis$data <- electrophoresis$data[electrophoresis$data[[x]] > 0,]
	if (log %in% c("y", "xy")) electrophoresis$data <- electrophoresis$data[electrophoresis$data[[y]] > 0,]
	
	# remove data outside the space between markers
	if (between.markers) for (i in 1:nrow(electrophoresis$peaks)) {
		if (electrophoresis$peaks$peak.observations[i] == "Lower Marker") {
			electrophoresis$data <- subset(electrophoresis$data, ! (
				batch == electrophoresis$peaks$batch[i] &
				well.number == electrophoresis$peaks$well.number[i] &
				distance >= 0.95 * electrophoresis$peaks$lower.distance[i]
			))
		} else if (electrophoresis$peaks$peak.observations[i] == "Upper Marker") {
			electrophoresis$data <- subset(electrophoresis$data, ! (
				batch == electrophoresis$peaks$batch[i] &
				well.number == electrophoresis$peaks$well.number[i] &
				distance <= 1.05 * electrophoresis$peaks$upper.distance[i]
			))
		}
	}
	
	# create plot but don't add the geom yet
	this.plot <- ggplot(electrophoresis$data)
	
	# add regions
	if (facet & ! is.na(region.alpha) & ! is.null(electrophoresis$regions)) this.plot <- this.plot + geom_rect(aes_(xmin = as.name(paste0("lower.", x)), xmax = as.name(paste0("upper.", x)), ymin = -Inf, ymax = Inf), data = electrophoresis$regions, alpha = region.alpha)
	
	# finally add the geom (after the regions so it's in front)
	this.plot <- this.plot + switch(geom,
		line = geom_line(if (facet)
			aes_(x = as.name(x), y = as.name(y)) 
		else 
			aes_(x = as.name(x), y = as.name(y), color = as.name("name"))
		),
		area = if (facet)
			geom_area(aes_(x = as.name(x), y = as.name(y))) 
		else
			geom_area(aes_(x = as.name(x), y = as.name(y), fill = as.name("name")), alpha = area.alpha)
	)
	
	# add peaks
	if (facet & ! is.na(peak.fill) & ! is.null(electrophoresis$peaks)) this.plot <- this.plot + geom_area(aes_(x = as.name(x), y = as.name(y), group = as.name("peak")), data = subset(electrophoresis$data, ! is.na(peak)), fill = peak.fill)
	
	# add faceting
	if (facet) this.plot <- this.plot + facet_wrap(~ batch * well.number, scales = scales, labeller = labeller.electrophoresis(electrophoresis))
	
	# apply log transformations
	if (log %in% c("x", "xy")) this.plot <- this.plot + scale_x_log10()
	if (log %in% c("y", "xy")) this.plot <- this.plot + scale_y_log10()
	
	# set labels and other settings for specific x & y variables
	this.plot <- this.plot + switch(x,
		length = xlab("length (bases)"),
		time = xlab("time (s)"),
		distance = xlab("distance migrated") + scale_x_reverse(),
		relative.distance = xlab("distance migrated relative to markers") + scale_x_reverse()
	)
	if (y == "delta.molarity") this.plot <- this.plot + ylab("molarity (pM)")
	
	this.plot
}

qc.mobility <- function(electrophoresis, n.simulate = 100, line.color = "red") { # returns a ggplot object, which can be extended by adding more features
	ladder.data <- subset(electrophoresis$data, sample.observations == "Ladder" & ! is.na(peak))
	ladder.data$true.length <- electrophoresis$peaks$length[ladder.data$peak]
	good.peaks <- subset(electrophoresis$peaks, ! is.na(length))
	simulated.data <- do.call(rbind, lapply(names(electrophoresis$wells.by.ladder), function(batch) do.call(rbind, lapply(names(electrophoresis$wells.by.ladder[[batch]]), function(ladder.well) {
		relative.distance.range <- range(subset(good.peaks, batch == batch & well.number == ladder.well)$relative.distance)
		relative.distance.diff <- diff(relative.distance.range)
		result <- data.frame(batch, well.number = ladder.well, relative.distance = relative.distance.range[1] + relative.distance.diff * (0:(n.simulate - 1) / (n.simulate - 1)))
		result$estimated.length <- electrophoresis$mobility.functions[[batch]][[ladder.well]](result$relative.distance)
		result
	}))))
	ggplot(ladder.data, aes(x = true.length, y = relative.distance, color = fluorescence)) +
		geom_point() + 
		geom_point(aes(x = length, y = relative.distance), data = subset(good.peaks, sample.observations == "Ladder"), color = line.color) + # overlay the reported peak positions
		geom_line(aes(x = estimated.length, y = relative.distance), data = simulated.data, col = line.color) + # overlay the simulated data from the mobility function
		scale_y_reverse() +
		xlab("true length (bases)") +
		ylab("distance migrated relative to markers") +
		facet_wrap(~ batch * well.number)
}

qc.molarity <- function(electrophoresis, log = TRUE) {
	peaks <- electrophoresis$peaks
	peaks$estimated.molarity <- sapply(1:nrow(peaks), function(peak.index) sum(electrophoresis$data$delta.molarity[which(electrophoresis$data$peak == peak.index)])) # without the which() you get the NA's too
	peaks <- subset(peaks, ! is.na(molarity) & ! is.na(estimated.molarity)) # remove NA's so they don't affect the x-limits and throw a warning
	
	result <- ggplot(peaks, aes(molarity, estimated.molarity)) +
		geom_point() +
		geom_abline() +
		geom_smooth(method = "lm") +
		xlab("software-reported molarity (pM)") +
		ylab("estimated molarity (pM)") +
		facet_wrap(~ batch * well.number, labeller = labeller.electrophoresis(electrophoresis))
	
	if (log) result <- result + scale_x_log10() + scale_y_log10()
	
	result
}
