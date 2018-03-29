library(ggplot2)

plot.electropherogram <- function(results, # returns a ggplot object, which can be extended by adding more features
	x = "length",
	y = "fluorescence",
	scales = "free_y", # scaling rules for the facets, passed to facet_wrap()
	geom = geom_line,
	include.ladder = FALSE
) {
	if (! include.ladder) results <- subset(results, name != "Ladder")

	this.plot <- ggplot(results) +
	aes_string(x, y) +
	geom() +
	facet_wrap(~ well.number, scales = scales)
	
	if (x == "length") this.plot <- this.plot + scale_x_log10()
	
	# set labels
	if (x == "length") this.plot <- this.plot + xlab("length (bases)")
	if (x == "time") this.plot <- this.plot + xlab("time (s)")
	if (x == "distance") this.plot <- this.plot + xlab("distance migrated") + scale_x_reverse()
	if (y == "delta.molarity") this.plot <- this.plot + ylab("concentration (pM)")
		
	this.plot
}

