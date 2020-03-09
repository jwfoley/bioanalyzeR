library(XML)
library(openssl)

read.bioanalyzer <- function(xml.file, fit = "spline") {
	stopifnot(fit %in% c("linear", "spline", "regression"))
	
	batch <- sub("\\.xml$", "", basename(xml.file))
	xml.root <- xmlRoot(xmlParse(xml.file))
	
	# read raw sample data
	samples <- xml.root[["Chips"]][["Chip"]][["Files"]][["File"]][["Samples"]]
	result.list <- xmlApply(samples, function(this.sample) {
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
			which.lower.marker <- which(peaks$peak.observations == "Lower Marker")
			stopifnot(length(which.lower.marker) == 1)
			which.upper.marker <- which(peaks$peak.observations == "Upper Marker")
			stopifnot(length(which.upper.marker) == 1)
			alignment.coefficient <- diff(peaks$aligned.time[c(which.lower.marker, which.upper.marker)]) / diff(peaks$time[c(which.lower.marker, which.upper.marker)])
			alignment.offset <- peaks$aligned.time[which.lower.marker] - alignment.coefficient * peaks$time[which.lower.marker]
			raw.data$aligned.time <- raw.data$time * alignment.coefficient + alignment.offset
			
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
	marker.aligned.times <- sapply(c("VirtualLowerMarkerTime", "VirtualUpperMarkerTime"), function (field) as.numeric(xmlValue(xml.root[["Chips"]][["Chip"]][["AssayBody"]][["DAAssaySetpoints"]][["DAMAssayInfoMolecular"]][[field]])))
	result$data$length <- sapply(result$data$aligned.time, function(x) if (x < marker.aligned.times[1] || x > marker.aligned.times[2]) NA else standard.curve.function(x)) # avoid extrapolating
	
	# convert to molarity
	peaks.ladder$normalized.area <- peaks.ladder$area / peaks.ladder$length / peaks.ladder$time
	fluorescence.coefficient <- lm(molarity ~ normalized.area - 1, peaks.ladder)$coefficients
	data.calibration <- cbind(result$data, do.call(rbind, lapply(result$samples$well.number, function(this.well) {
		result.this.well <- subset(result$data, well.number == this.well)
		data.frame(
			delta.fluorescence = c(NA, diff(result.this.well$fluorescence)),
			delta.aligned.time = c(NA, diff(result.this.well$aligned.time)),
			delta.length = c(NA, diff(result.this.well$length))
		)
		})))
	data.calibration$delta.area <- (2 * data.calibration$fluorescence - data.calibration$delta.fluorescence) / 2 * data.calibration$delta.aligned.time
	result$data$molarity <- data.calibration$delta.area * fluorescence.coefficient / data.calibration$aligned.time / data.calibration$delta.length * data.calibration$delta.aligned.time / data.calibration$length
	
	# annotate which peak each data point is in, if any
	# WARNING: if peaks overlap, this will overwrite and each point will only be mapped to the last-occuring one!
	result$data$peak <- NA
	for (i in 1:nrow(result$peaks)) result$data$peak[result$data$well.number == result$peaks$well.number[i] & result$data$aligned.time >= result$peaks$lower.aligned.time[i] & result$data$aligned.time <= result$peaks$upper.aligned.time[i]] <- i
	
	# construct well.by.ladder (analogous to TapeStation but not as useful here)
	wells.by.ladder <- list(list(result$samples$well.number))
	names(wells.by.ladder) <- batch
	names(wells.by.ladder[[1]]) <- result$samples$well.number[which.ladder]
	
	# construct mobility.functions
	mobility.functions <- list(list(standard.curve.function))
	names(mobility.functions) <- batch
	names(mobility.functions[[1]]) <- result$samples$well.number[which.ladder]
	
	structure(list(
		data = result$data,
		samples = result$samples,
		wells.by.ladder = wells.by.ladder,
		peaks = result$peaks,
		regions = NULL,
		mobility.functions = mobility.functions,
		mass.coefficients = list(fluorescence.coefficient)
	), class = "electrophoresis")
}

