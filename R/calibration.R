calculate.length <- function(electrophoresis, fit = "spline") {
	x.name <- get.x.name(electrophoresis)
	lower.name <- paste0("lower.", x.name)
	upper.name <- paste0("upper.", x.name)
		
	# prepare new fields to be filled in piecemeal from each ladder
	electrophoresis$data$length <- NA
	electrophoresis$peaks$lower.length <- NA
	electrophoresis$peaks$upper.length <- NA
	if (! is.null(electrophoresis$regions)) {
		electrophoresis$regions[[lower.name]] <- NA
		electrophoresis$regions[[upper.name]] <- NA
	}
	electrophoresis$mobility.functions <- list()
	
	# fit a mobility model for each ladder and apply it to the appropriate samples
	for (batch in unique(electrophoresis$samples$batch)) {
		electrophoresis$mobility.functions[[batch]] <- list()
		in.this.batch <- electrophoresis$samples$batch == batch
		for (ladder.well in unique(electrophoresis$samples$ladder.well[which(in.this.batch)])) {
			which.ladder.index <- which(in.this.batch & electrophoresis$samples$well.number == ladder.well)
			peaks.ladder <- subset(electrophoresis$peaks, sample.index == which.ladder.index)
			which.samples <- which(in.this.batch & electrophoresis$samples$ladder.well == ladder.well)
			which.rows <- which(electrophoresis$data$sample.index %in% which.samples)
			which.peaks <- which(electrophoresis$peaks$sample.index %in% which.samples)
			which.regions <- which(electrophoresis$regions$sample.index %in% which.samples)
			
			peaks.ladder$x <- peaks.ladder[[x.name]]
			
			# fit standard curve for molecule length vs. x-value
			# do this in relative x space so it's effectively recalibrated for each sample's markers
			if (fit == "interpolation") {
				warning("linear interpolation gives ugly results for molarity estimation")
				standard.curve.function <- approxfun(peaks.ladder$x, peaks.ladder$length)
				standard.curve.inverse <- approxfun(peaks.ladder$length, peaks.ladder$x)
			} else if (fit == "spline") {
				standard.curve.function <- splinefun(peaks.ladder$x, peaks.ladder$length, method = "natural")
				standard.curve.inverse <- splinefun(peaks.ladder$length, peaks.ladder$x, method = "natural")
			} else if (fit == "regression") {
				if (x.name == "relative.distance") {
					mobility.model <- lm(x ~ log(length), peaks.ladder)
					standard.curve.function <- function(x) exp((x - mobility.model$coefficients[1]) / mobility.model$coefficients[2])
					standard.curve.inverse <- function(length) mobility.model$coefficients[1] + mobility.model$coefficients[2] * log(length)
				} else if (x.name == "aligned.time") {
					# if x-variable is time, we must correct for the fact that faster-moving molecules spend less time in front of the detector
					mobility.model <- lm(1/aligned.time ~ log(length), data = peaks.ladder)
					standard.curve.function <- function(aligned.time) exp((1 / aligned.time - mobility.model$coefficients[1]) / mobility.model$coefficients[2])
					standard.curve.inverse <- function(length) 1/(mobility.model$coefficients[1] + log(length) * mobility.model$coefficients[2])
				}
			}
			electrophoresis$mobility.functions[[1]][[ladder.well]] <- standard.curve.function
			
			# apply model to raw data
			electrophoresis$data$length[which.rows] <- standard.curve.function(electrophoresis$data[[x.name]][which.rows])
			electrophoresis$data$length[! (
				in.custom.region(electrophoresis$data, min(peaks.ladder$length), max(peaks.ladder$length)) &
				in.custom.region(electrophoresis$data, min(peaks.ladder$x), max(peaks.ladder$x), bound.variable = x.name)
			)] <- NA # avoid extrapolation
			
			# apply model to peaks
			electrophoresis$peaks$lower.length[which.peaks] <- standard.curve.function(electrophoresis$peaks[[upper.name]][which.peaks])
			electrophoresis$peaks$upper.length[which.peaks] <- standard.curve.function(electrophoresis$peaks[[lower.name]][which.peaks])
			
			# apply inverse model to regions
			if (! is.null(electrophoresis$regions)) {
				electrophoresis$regions[[lower.name]][which.regions] <- standard.curve.inverse(electrophoresis$regions$upper.length[which.regions])
				electrophoresis$regions[[upper.name]][which.regions] <- standard.curve.inverse(electrophoresis$regions$lower.length[which.regions])
			}
		}
	}
	
	electrophoresis
}


calculate.concentration <- function(electrophoresis, ladder.concentrations = NULL) {
	x.name <- get.x.name(electrophoresis)
	delta <- do.call(rbind, by(electrophoresis$data, electrophoresis$data$sample.index, function(data.subset) data.frame(
		fluorescence = c(NA, diff(data.subset$fluorescence)),
		x = c(NA, diff(data.subset[[x.name]])) # these will probably be constant but we'd better not assume
	), simplify = F))
	if (x.name == "relative.distance") delta$x <- -delta$x # distances are stored in decreasing order so the deltas are negative
	# estimate area under each measurement with the trapezoidal rule; to simplify math, each point's sum is for the trapezoid to the left of it
	area <- (2 * electrophoresis$data$fluorescence - delta$fluorescence) * delta$x
	if (x.name == "aligned.time") area <- area / electrophoresis$data[[x.name]] # compensate for time spent in the detector (faster-moving molecules get less signal)
	
	# calculate coefficients relating concentration to area under the curve
	has.upper.marker <- any(electrophoresis$peaks$peak.observations %in% UPPER.MARKER.NAMES)
	mass.coefficients <- if (! is.null(ladder.concentrations)) {
		# if ladder concentrations are given, fit a mass coefficient on the non-marker ladder peaks
		# exclude markers because they are at extreme values and Bioanalyzer RNA kits don't report the concentration of their marker anyway
		which.markers <- lapply(1:nrow(electrophoresis$samples), function(sample.index) which(
			electrophoresis$peaks$sample.index == sample.index & (
				(
					electrophoresis$peaks$peak.observations %in% LOWER.MARKER.NAMES & 
					electrophoresis$peaks$concentration == ladder.concentrations[1] # verify this is the right peak (sometimes Bioanalyzer annotates more than one as the marker but it only assigns the hardcoded concentration to one)
				) | (
					electrophoresis$peaks$peak.observations %in% UPPER.MARKER.NAMES &
					electrophoresis$peaks$concentration == ladder.concentrations[length(ladder.concentrations)]
				)
			)
		))
		marker.areas <- lapply(which.markers, function(peak.indexes) sapply(peak.indexes, function(peak.index) sum(area[which(in.peak(electrophoresis, peak.index))])))
		result <- rep(NA, nrow(electrophoresis$samples))
		for (ladder.well in unique(electrophoresis$samples$ladder.well)) {
			ladder.index <- which(electrophoresis$samples$well.number == ladder.well)
			stopifnot(length(ladder.index) == 1)
			which.samples <- which(electrophoresis$samples$ladder.well == ladder.well)
			non.marker.concentrations <- ladder.concentrations[-c(1, if (has.upper.marker) length(ladder.concentrations) else NULL)]
			which.ladder.peaks <- which(
				electrophoresis$peaks$sample.index == ladder.index &
				(! electrophoresis$peaks$peak.observations %in% c(LOWER.MARKER.NAMES, UPPER.MARKER.NAMES)) &
				electrophoresis$peaks$concentration %in% non.marker.concentrations
			)
			ladder.mass.coefficient <- mean(electrophoresis$peaks$concentration[which.ladder.peaks] / sapply(which.ladder.peaks, function(i) sum(area[which(in.peak(electrophoresis, i))])))
			
			# then modify the mass coefficient for each sample according to the ratio of its marker area(s) to the ladder's
			result[which.samples] <- ladder.mass.coefficient * sapply(marker.areas[which.samples], function(these.areas) mean(marker.areas[[ladder.index]] / these.areas)) # if there are two markers per sample, this gives the mean area ratio relative to their counterparts in the ladder well; if only one marker, the mean ratio is just the ratio 
		}
		
		result
		
	} else {
		# without known ladder concentrations, just calibrate by the upper marker if present, lower marker otherwise
		# this is the TapeStation's approach and that's the only known concentration it reports so it's the only way
		which.markers <- sapply(1:nrow(electrophoresis$samples), function(sample.index) which(
			electrophoresis$peaks$sample.index == sample.index & (
				(has.upper.marker & electrophoresis$peaks$peak.observations %in% UPPER.MARKER.NAMES) |
				(! has.upper.marker & electrophoresis$peaks$peak.observations %in% LOWER.MARKER.NAMES)
			)
		))
		marker.concentrations <- electrophoresis$peaks$concentration[which.markers]
		stopifnot(all(marker.concentrations == marker.concentrations[1]))
		marker.concentrations / sapply(which.markers, function(i) sum(area[which(in.peak(electrophoresis, i))]))
	}
	
	electrophoresis$data$concentration <- mass.coefficients[electrophoresis$data$sample.index] * area	
	electrophoresis
}

