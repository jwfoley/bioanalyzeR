library(ggplot2)

plot.bioanalyzer <- function(results, # returns a ggplot object, which can be extended by adding more features
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
	facet_wrap(~ name, scales = scales)
	
	if (x == "length") this.plot <- this.plot + scale_x_log10()
	
	# set labels
	if (x == "length") this.plot <- this.plot + xlab("length (bases)")
	if (x == "time") this.plot <- this.plot + xlab("time (s)")
	if (y == "molarity") this.plot <- this.plot + ylab("concentration (pM)")
		
	this.plot
}

