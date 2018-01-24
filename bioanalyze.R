library(XML)
library(openssl)
library(ggplot2)

read.bioanalyzer <- function(xml.file, fit = "linear") {
	stopifnot(fit %in% c("linear", "spline"))

	xml.root <- xmlRoot(xmlParse(xml.file))
	
	# extract ladder specifications
	ladder <- data.frame(t(xmlSApply(xml.root[["Chips"]][["Chip"]][["AssayBody"]][["DAAssaySetpoints"]][["DAMAssayInfoMolecular"]][["LadderPeaks"]], function(standard) c(as.numeric(xmlValue(standard[["Size"]])), as.numeric(xmlValue(standard[["Concentration"]]))))), row.names = NULL)
	names(ladder) <- c("length", "concentration")
	ladder$molarity <- ladder$concentration / ladder$length / 0.00066 # magic number computed from Agilent standard curves

	# extract samples
	samples <- xml.root[["Chips"]][["Chip"]][["Files"]][["File"]][["Samples"]]
	results <- do.call(rbind, xmlApply(samples, function(this.sample) {
		if (xmlValue(this.sample[["HasData"]]) == "false") NULL else {
			signal.data <- this.sample[["DASignals"]][["DetectorChannels"]][[1]][["SignalData"]]
			n.values <- as.integer(xmlValue(signal.data[["NumberOfSamples"]]))
			result <- data.frame(
				index = as.integer(xmlValue(this.sample[["Index"]])),
				name = xmlValue(this.sample[["Name"]]),
				time = as.numeric(xmlValue(signal.data[["XStart"]])) + as.numeric(xmlValue(signal.data[["XStep"]])) * (1:n.values - 1),
				fluorescence = readBin(base64_decode(xmlValue(signal.data[["ProcessedSignal"]])), "numeric", size = 4, n = n.values, endian = "little"),
				stringsAsFactors = F
			)
			
			# fit standard curve for molecule length (done once per sample because that's how the results are provided, even though it should be the same standard curve for every sample)
			ladder.peaks <- this.sample[["DAResultStructures"]][["DARStandardCurve"]][[1]][["LadderPeaks"]]
			stopifnot(all.equal(xmlSApply(ladder.peaks, function(peak) as.numeric(xmlValue(peak[["FragmentSize"]]))), ladder$length, check.names = F)) # verify it matches the known ladder
			this.ladder <- cbind(ladder, data.frame(t(xmlSApply(ladder.peaks, function(peak) c(as.numeric(xmlValue(peak[["MigrationTime"]])), as.numeric(xmlValue(peak[["Area"]]))))), row.names = NULL))
			names(this.ladder)[ncol(ladder) + 1:2] <- c("time", "area")
			
			if (fit == "linear") {
				standard.curve.function <- approxfun(this.ladder$time, this.ladder$length)		
			} else if (fit == "spline") {
				standard.curve.function <- splinefun(this.ladder$time, this.ladder$length, method = "natural")
			}
			result$length <- standard.curve.function(result$time)
			
			# convert to molarity
			this.ladder$normalized.fluorescence <- this.ladder$area / this.ladder$length / this.ladder$time
			fluorescence.coefficient <- lm(molarity ~ normalized.fluorescence - 1, this.ladder)$coefficients
			length.derivative <- diff(result$length) / diff(result$time)
			result$molarity <- result$fluorescence * fluorescence.coefficient / result$time / c(NA, length.derivative) / result$length
		
			result
		}
	}))
	
	# format output nicely
	rownames(results) <- NULL
	results$name <- factor(results$name, levels = unique(results$name)) # make names into a factor that keeps them in the observed order
	
	results
}

plot.bioanalyzer <- function(results, # returns a ggplot object, which can be extended by adding more features
	x = "length",
	y = "fluorescence",
	scales = "free_y", # scaling rules for the facets, passed to facet_wrap()
	geom = geom_line,
	include_ladder = FALSE
) {
	if (! include_ladder) results <- subset(results, name != "Ladder")

	this.plot <- ggplot(results) +
	aes_string(x, y) +
	geom() +
	facet_wrap(~ name, scales = scales)
	
	if (x == "length") this.plot <- this.plot + scale_x_log10()
		
	this.plot
}

