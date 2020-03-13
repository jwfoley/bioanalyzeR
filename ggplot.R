library(ggplot2)

# given an electrophoresis object and a variable name (like "length"), generate a nice axis label with appropriate units
axis.label <- function(electrophoresis, variable) switch(variable,
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
	between.markers = TRUE, # show only data between the marker peaks (so their signal doesn't ruin the scale)
	peak.fill = "darkred", # set to NA to stop showing peaks
	region.alpha = 0.2, # set to NA to stop showing regions
	area.alpha = 0.2 # if facet = FALSE and geom = "area", alpha transparency of the filled areas
) {

	if (is.null(log)) log <- if (x == "length") "x" else NA
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
	this.plot <- this.plot + xlab(axis.label(electrophoresis, x)) + ylab(axis.label(electrophoresis, y))
	if (x %in% c("distance", "relative.distance")) this.plot <- this.plot + scale_x_reverse()
	
	this.plot
}

qc.stdcrv <- function(electrophoresis, n.simulate = 100, line.color = "red") { # returns a ggplot object, which can be extended by adding more features
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
		xlab(paste("true", axis.label(electrophoresis, "length"))) +
		ylab(axis.label(electrophoresis, x.name)) +
		facet_wrap(~ batch * well.number)
	if (x.name == "relative.distance") this.plot <- this.plot + scale_y_reverse()
	
	this.plot
}

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
		xlab(paste("software-reported", axis.label(electrophoresis, variable))) +
		ylab(paste("estimated", axis.label(electrophoresis, variable))) +
		facet_wrap(~ batch * well.number, labeller = labeller.electrophoresis(electrophoresis))
	
	if (log) result <- result + scale_x_log10() + scale_y_log10()
	
	result
}

