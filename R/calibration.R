calculate.concentration <- function(electrophoresis, ladder.concentrations = NULL) {
	raw.x.name <- get.x.name(electrophoresis, T)
	data.calibration <- cbind(electrophoresis$data, do.call(rbind, by(electrophoresis$data, electrophoresis$data$sample.index, function(data.subset) data.frame(
		delta.fluorescence = c(NA, diff(data.subset$fluorescence)),
		delta.x = c(NA, diff(data.subset[[raw.x.name]]))
	), simplify = F)))
	if (raw.x.name == "distance") data.calibration$delta.x <- -data.calibration$delta.x # distances are stored in decreasing so the deltas are negative
	# estimate area under each measurement with the trapezoidal rule; to simplify math, each point's sum is for the trapezoid to the left of it
	data.calibration$area <- (2 * data.calibration$fluorescence - data.calibration$delta.fluorescence) * data.calibration$delta.x
	if (raw.x.name == "time") data.calibration$area <- data.calibration$area / data.calibration$time # compensate for time spent in the detector (faster-moving molecules get less signal)
	
	# calculate coefficients relating concentration to area under the curve
	mass.coefficients <- if (! is.null(ladder.concentrations)) {
		# if ladder concentrations are given, fit a mass coefficient on the non-marker ladder peaks
		which.markers <- lapply(1:nrow(electrophoresis$samples), function(sample.index) which(
			electrophoresis$peaks$sample.index == sample.index & (
				(
					electrophoresis$peaks$peak.observations == "Lower Marker" & 
					electrophoresis$peaks$concentration == ladder.concentrations[1] # verify this is the right peak (sometimes Bioanalyzer annotates more than one as the marker but it only assigns the hardcoded concentration to one)
				) | (
					electrophoresis$peaks$peak.observations == "Upper Marker" &
					electrophoresis$peaks$concentration == ladder.concentrations[length(ladder.concentrations)]
				)
			)
		))
		marker.areas <- lapply(which.markers, function(peak.indexes) sapply(peak.indexes, function(peak.index) sum(data.calibration$area[which(in.peak(electrophoresis, peak.index))])))
		result <- rep(NA, nrow(electrophoresis$samples))
		for (ladder.well in unique(electrophoresis$samples$ladder.well)) {
			ladder.index <- which(electrophoresis$samples$well.number == ladder.well)
			stopifnot(length(ladder.index) == 1)
			which.samples <- which(electrophoresis$samples$ladder.well == ladder.well)
			which.ladder.peaks <- which(
				electrophoresis$peaks$sample.index == ladder.index &
				electrophoresis$peaks$concentration %in% ladder.concentrations
			)
			which.ladder.peaks <- which.ladder.peaks[! electrophoresis$peaks$peak.observations[which.ladder.peaks] %in% c(LOWER.MARKER.NAMES, UPPER.MARKER.NAMES)]
			ladder.mass.coefficient <- mean(electrophoresis$peaks$concentration[which.ladder.peaks] / sapply(which.ladder.peaks, function(i) sum(data.calibration$area[which(in.peak(electrophoresis, i))])))
			
			# then modify the mass coefficient for each sample according to the ratio of its marker area(s) to the ladder's
			result[which.samples] <- ladder.mass.coefficient * sapply(marker.areas[which.samples], function(these.areas) mean(marker.areas[[ladder.index]] / these.areas)) # if there are two markers per sample, this gives the mean area ratio relative to their counterparts in the ladder well; if only one marker, the mean ratio is just the ratio 
		}
		
		result
		
	} else {
		# without known ladder concentrations, just calibrate by the upper marker if present, lower marker otherwise
		# this is the TapeStation's approach and that's the only known concentration it reports so it's the only way
		has.upper.marker <- any(electrophoresis$peaks$peak.observations %in% UPPER.MARKER.NAMES)
		which.markers <- sapply(1:nrow(electrophoresis$samples), function(sample.index) which(
			electrophoresis$peaks$sample.index == sample.index & (
				(has.upper.marker & electrophoresis$peaks$peak.observations %in% UPPER.MARKER.NAMES) |
				(! has.upper.marker & electrophoresis$peaks$peak.observations %in% LOWER.MARKER.NAMES)
			)
		))
		marker.concentrations <- electrophoresis$peaks$concentration[which.markers]
		stopifnot(all(marker.concentrations == marker.concentrations[1]))
		marker.concentrations / sapply(which.markers, function(i) sum(data.calibration$area[which(in.peak(electrophoresis, i))]))
	}

	mass.coefficients[electrophoresis$data$sample.index] * data.calibration$area	
}

