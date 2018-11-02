library(ggplot2)

plot.electropherogram <- function(results, # returns a ggplot object, which can be extended by adding more features
	x = "length",
	scales = "free_y", # scaling rules for the facets, passed to facet_wrap()
	geom = geom_line,
	include.ladder = FALSE
) {
	
	if (! include.ladder) results <- subset(results, category != "Ladder")

	this.plot <- ggplot(results) +
		aes_(as.name(x), ~fluorescence) +
		geom() +
		facet_wrap(~ well.number, scales = scales) +
		ylab("fluorescence")
	
	if (x == "length") this.plot <- this.plot + scale_x_log10()
	
	# set labels
	if (x == "length") this.plot <- this.plot + xlab("length (bases)")
	if (x == "time") this.plot <- this.plot + xlab("time (s)")
	if (x == "distance") this.plot <- this.plot + xlab("distance migrated") + scale_x_reverse()
		
	this.plot
}


plot.molarity <- function(results, # returns a ggplot object, which can be extended by adding more features
	x = "length",
	scales = "free_y",
	include.ladder = FALSE
) {
	
	if (! include.ladder) results <- subset(results, category != "Ladder")
	
	results$xmin <- unlist(lapply(unique(results$well.number), function(well) {
		results.this.well <- subset(results, well.number == well)
		c(NA, results.this.well[1:(nrow(results.this.well) - 1), x])
	}))
		
	this.plot <- ggplot(results) +
		aes_(xmin = ~xmin, xmax = as.name(x), ymin = 0, ymax = ~delta.molarity) +
		geom_rect() +
		facet_wrap(~ well.number, scales = scales) +
		ylab("concentration (pM)")
	
	if (x == "length") this.plot <- this.plot + scale_x_log10()
	
	# set labels
	if (x == "length") this.plot <- this.plot + xlab("length (bases)")
	if (x == "time") this.plot <- this.plot + xlab("time (s)")
	if (x == "distance") this.plot <- this.plot + xlab("distance migrated") + scale_x_reverse()
	
	this.plot
}

