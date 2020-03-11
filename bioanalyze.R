library(XML)
library(openssl)

read.bioanalyzer <- function(xml.file, fit = "spline") {
	stopifnot(fit %in% c("linear", "spline", "regression"))
	
	batch <- sub("\\.xml$", "", basename(xml.file))
	xml.root <- xmlRoot(xmlParse(xml.file))
	has.upper.marker <- as.logical(xmlValue(xml.root[["Chips"]][["Chip"]][["AssayBody"]][["DASampleSetpoints"]][["DAMAlignment"]][["Channel"]][["AlignUpperMarker"]]))
	defined.ladder.peaks <- xmlToDataFrame(xml.root[["Chips"]][["Chip"]][["AssayBody"]][["DAAssaySetpoints"]][["DAMAssayInfoMolecular"]][["LadderPeaks"]], colClasses = rep("numeric", 4))
	
	# read raw sample data
	result.list <- xmlApply(xml.root[["Chips"]][["Chip"]][["Files"]][["File"]][["Samples"]], function(this.sample) {
		if (xmlValue(this.sample[["HasData"]]) != "true") NULL else {
			# read metadata
			well.number <- as.integer(xmlValue(this.sample[["WellNumber"]]))
			sample.name <- trimws(xmlValue(this.sample[["Name"]]))
			is.ladder <- (trimws(xmlValue(this.sample[["Category"]])) == "Ladder")
			sample.observations <- trimws(xmlValue(this.sample[["Comment"]]))

			# read peaks
			peaks.raw <- xmlToDataFrame(this.sample[["DAResultStructures"]][["DARIntegrator"]][["Channel"]][["PeaksMolecular"]], stringsAsFactors = F)
			peaks <- data.frame(
				peak.observations =   trimws(peaks.raw$Observations),
				length =              as.numeric(peaks.raw$FragmentSize),
				time =                as.numeric(peaks.raw$MigrationTime),
				aligned.time =        as.numeric(peaks.raw$AlignedMigrationTime),
				lower.time =          as.numeric(peaks.raw$StartTime),
				upper.time =          as.numeric(peaks.raw$EndTime),
				lower.aligned.time =  as.numeric(peaks.raw$AlignedStartTime),
				upper.aligned.time =  as.numeric(peaks.raw$AlignedEndTime),
				area =                as.numeric(peaks.raw$Area),
				molarity =            as.numeric(peaks.raw$Molarity),
				stringsAsFactors =    F
			)
		
			# read signal
			signal.data <- this.sample[["DASignals"]][["DetectorChannels"]][[1]][["SignalData"]]
			n.values <- as.integer(xmlValue(signal.data[["NumberOfSamples"]]))
			raw.data <- data.frame(
				time = as.numeric(xmlValue(signal.data[["XStart"]])) + as.numeric(xmlValue(signal.data[["XStep"]])) * (1:n.values - 1),
				fluorescence = readBin(base64_decode(xmlValue(signal.data[["ProcessedSignal"]])), "numeric", size = 4, n = n.values, endian = "little"),
				stringsAsFactors = F
			)
			
			# align the observation times according to the markers in this sample		
			which.lower.marker <- which(peaks$peak.observations == "Lower Marker" & peaks$length == defined.ladder.peaks$Size[1]) # check the size too because sometimes the software annotates more than one as the same marker with no consequences
			stopifnot(length(which.lower.marker) == 1)
			if (has.upper.marker) {
				which.upper.marker <- which(peaks$peak.observations == "Upper Marker" & peaks$length == defined.ladder.peaks$Size[nrow(defined.ladder.peaks)])
				stopifnot(length(which.upper.marker) == 1)
				alignment.coefficient <- diff(peaks$aligned.time[c(which.lower.marker, which.upper.marker)]) / diff(peaks$time[c(which.lower.marker, which.upper.marker)])
				alignment.offset <- peaks$aligned.time[which.lower.marker] - alignment.coefficient * peaks$time[which.lower.marker]
				raw.data$aligned.time <- raw.data$time * alignment.coefficient + alignment.offset
			} else {
				alignment.coefficient <- peaks$aligned.time[which.lower.marker] / peaks$time[which.lower.marker]
				raw.data$aligned.time <- raw.data$time * alignment.coefficient
			}
			
			list(
				data = data.frame(batch, well.number, sample.name, is.ladder, sample.observations, raw.data, stringsAsFactors = F),
				samples = data.frame(batch, well.number, sample.name, is.ladder, sample.observations, stringsAsFactors = F),
				peaks = data.frame(batch, well.number, sample.name, is.ladder, sample.observations, peaks, stringsAsFactors = F)
			)
		}
	})
	result <- lapply(c("data", "samples", "peaks"), function(item) do.call(rbind, c(lapply(result.list, function(x) x[[item]]), make.row.names = F)))
	names(result) <- c("data", "samples", "peaks")
	
	# convert sample metadata into factors, ensuring all frames have the same levels and the levels are in the observed order
	for (field in c("batch", "well.number", "sample.name", "sample.observations")) {
		result$samples[,field] <- factor(result$samples[,field], levels = unique(result$samples[,field]))
		if (field %in% names(result$data)) result$data[,field] <- factor(result$data[,field], levels = levels(result$samples[,field]))
		if (field %in% names(result$peaks)) result$peaks[,field] <- factor(result$peaks[,field], levels = levels(result$samples[,field]))
	}
	# convert other text into factors without those restrictions
	result$peaks[,"peak.observations"] <- factor(result$peaks[,"peak.observations"])
	
	# analyze ladder
	which.ladder <- which(result$samples$is.ladder)
	stopifnot(length(which.ladder) == 1)
	peaks.ladder <- subset(result$peaks, well.number == result$samples$well.number[which.ladder])
	
	# fit standard curve for molecule length
	# do this with aligned times so it's effectively recalibrated for each sample's markers
	if (fit == "interpolate") {
		warning("linear interpolation gives ugly results for molarity estimation")
		standard.curve.function <- approxfun(peaks.ladder$aligned.time, peaks.ladder$length)
	} else if (fit == "spline") {
		standard.curve.function <- splinefun(peaks.ladder$aligned.time, peaks.ladder$length, method = "natural")
	} else if (fit == "regression") {
		mobility.model <- lm(1/aligned.time ~ log(length), data = peaks.ladder)
		standard.curve.function <- function(aligned.time) exp((1 / aligned.time - mobility.model$coefficients[1]) / mobility.model$coefficients[2])
	}
	result$data$length <- standard.curve.function(result$data$aligned.time)
	ladder.limits <- range(defined.ladder.peaks$Size)
	result$data$length[result$data$length < ladder.limits[1] | result$data$length > ladder.limits[2]] <- NA # avoid extrapolating
	result$peaks$lower.length <- standard.curve.function(result$peaks$lower.aligned.time)
	result$peaks$upper.length <- standard.curve.function(result$peaks$upper.aligned.time)
	
	# annotate which peak each data point is in, if any
	# WARNING: if peaks overlap, this will overwrite and each point will only be mapped to the last-occuring one!
	result$data$peak <- NA
	for (i in 1:nrow(result$peaks)) result$data$peak[result$data$well.number == result$peaks$well.number[i] & result$data$aligned.time >= result$peaks$lower.aligned.time[i] & result$data$aligned.time <= result$peaks$upper.aligned.time[i]] <- i
	
	# convert to molarity
	# the idea is that we must correct fluorescence area by migration time to account for the fact that faster-moving molecules spend less time in front of the detector (Agilent's TimeCorrectedArea apparently does this with the raw time, not the aligned time)
	# and then fluorescence is proportional to mass, which is molarity * length
	data.calibration <- cbind(result$data, do.call(rbind, lapply(result$samples$well.number, function(this.well) {
		result.this.well <- subset(result$data, well.number == this.well)
		data.frame(
			delta.fluorescence = c(NA, diff(result.this.well$fluorescence)),
			delta.time = c(NA, diff(result.this.well$time))
		)
	})))
	# estimate area under each measurement with the trapezoidal rule; to simplify math, each point's sum is for the trapezoid to the left of it
	data.calibration$area <- (2 * data.calibration$fluorescence - data.calibration$delta.fluorescence) * data.calibration$delta.time
	# correct area by migration time
	data.calibration$corrected.area <- data.calibration$area / data.calibration$time
	# calculate corrected area under peaks
	peaks.calibration <- cbind(result$peaks, mass = result$peaks$molarity * result$peaks$length)	
	# fit the coefficient of mass vs. corrected area, independently for each sample, according to the markers
	mass.coefficients <- sapply(result$samples$well.number, function(well) {
		which.markers <- which(result$peaks$well.number == well & ((
			result$peaks$peak.observations == "Lower Marker" & result$peaks$length == defined.ladder.peaks$Size[1]
		) | (
			result$peaks$peak.observations == "Upper Marker" & result$peaks$length == defined.ladder.peaks$Size[nrow(defined.ladder.peaks)]
		)))
		stopifnot(
			(has.upper.marker && length(which.markers) == 2) ||
			(! has.upper.marker && length(which.markers) == 1)
		)
		marker.areas <- sapply(which.markers, function(peak) sum(data.calibration$corrected.area[which(data.calibration$peak == peak)]))
		marker.masses <- peaks.calibration$mass[which.markers]
		lm(marker.masses ~ marker.areas - 1)$coefficients[1]
	})
	names(mass.coefficients) <- result$samples$well.number
	# finally apply this coefficient to get the molarity of each trapezoid
	result$data$molarity <- data.calibration$corrected.area * mass.coefficients[result$data$well.number] / result$data$length
	
	# construct well.by.ladder (analogous to TapeStation but not as useful here)
	wells.by.ladder <- list(list(result$samples$well.number))
	names(wells.by.ladder) <- batch
	names(wells.by.ladder[[1]]) <- result$samples$well.number[which.ladder]
	
	# construct mobility.functions
	mobility.functions <- list(list(standard.curve.function))
	names(mobility.functions) <- batch
	names(mobility.functions[[1]]) <- result$samples$well.number[which.ladder]
	
	# construct mass coefficients
	mass.coefficients <- list(mass.coefficients)
	names(mass.coefficients) <- batch
	
	structure(list(
		data = result$data,
		samples = result$samples,
		wells.by.ladder = wells.by.ladder,
		peaks = result$peaks,
		regions = NULL,
		mobility.functions = mobility.functions,
		mass.coefficients = mass.coefficients
	), class = "electrophoresis")
}

