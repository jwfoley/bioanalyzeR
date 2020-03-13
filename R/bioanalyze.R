#' @export
#' @importFrom XML xmlRoot xmlParse xmlValue xmlToDataFrame xmlApply
#' @importFrom openssl base64_decode
read.bioanalyzer <- function(xml.file, fit = "spline") {
	stopifnot(fit %in% c("linear", "spline", "regression"))
	
	batch <- sub("\\.xml$", "", basename(xml.file))
	xml.root <- xmlRoot(xmlParse(xml.file))
	chip.root <- xml.root[["Chips"]][["Chip"]]
	assay.info <- list(
		file.name =           xmlValue(chip.root[["Files"]][["File"]][["FileInformation"]][["FileName"]]),
		creation.date =       xmlValue(chip.root[["ChipInformation"]][["CreationDate"]]),
		assay.name =          xmlValue(chip.root[["AssayHeader"]][["Title"]]),
		assay.type =          xmlValue(chip.root[["AssayHeader"]][["Class"]]),
		length.unit =         xmlValue(chip.root[["AssayBody"]][["DAAssaySetpoints"]][["DAMAssayInfoMolecular"]][["SizeUnit"]]),
		concentration.unit =  xmlValue(chip.root[["AssayBody"]][["DAAssaySetpoints"]][["DAMAssayInfoMolecular"]][["ConcentrationUnit"]])
	)
	stopifnot(assay.info$assay.type %in% names(molecular.weight))
	# hardcode the molarity unit depending on concentration unit, otherwise the MW scales will be wrong
	assay.info$molarity.unit <- switch(assay.info$concentration.unit,
		"ng/µl" = "nM",
		"pg/µl" = "pM"
	)
	has.upper.marker <- as.logical(xmlValue(chip.root[["AssayBody"]][["DASampleSetpoints"]][["DAMAlignment"]][["Channel"]][["AlignUpperMarker"]]))
	defined.ladder.peaks <- xmlToDataFrame(chip.root[["AssayBody"]][["DAAssaySetpoints"]][["DAMAssayInfoMolecular"]][["LadderPeaks"]], colClasses = rep("numeric", 4))
	
	# read raw sample data
	result.list <- xmlApply(chip.root[["Files"]][["File"]][["Samples"]], function(this.sample) {
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
				concentration =       as.numeric(peaks.raw$Concentration),
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
			which.lower.marker <- which(peaks$peak.observations == "Lower Marker" & peaks$concentration == defined.ladder.peaks$Concentration[1]) # check the concentration too because sometimes the software annotates more than one as the same marker with no consequences, and sometimes the size is off by a tiny bit, but the concentration is hardcoded
			stopifnot(length(which.lower.marker) == 1)
			if (has.upper.marker) {
				which.upper.marker <- which(peaks$peak.observations == "Upper Marker" & peaks$concentration == defined.ladder.peaks$Concentration[nrow(defined.ladder.peaks)])
				stopifnot(length(which.upper.marker) == 1)
				alignment.coefficient <- diff(peaks$aligned.time[c(which.lower.marker, which.upper.marker)]) / diff(peaks$time[c(which.lower.marker, which.upper.marker)])
				alignment.offset <- peaks$aligned.time[which.lower.marker] - alignment.coefficient * peaks$time[which.lower.marker]
			} else {
				alignment.coefficient <- peaks$aligned.time[which.lower.marker] / peaks$time[which.lower.marker]
				alignment.offset <- 0
			}
			raw.data$aligned.time <- raw.data$time * alignment.coefficient + alignment.offset
			
			list(
				data = data.frame(batch, well.number, sample.name, is.ladder, sample.observations, raw.data, stringsAsFactors = F),
				samples = data.frame(batch, well.number, sample.name, is.ladder, sample.observations, stringsAsFactors = F),
				peaks = data.frame(batch, well.number, sample.name, is.ladder, sample.observations, peaks, stringsAsFactors = F),
				alignment.values = c(alignment.coefficient, alignment.offset)
			)
		}
	})
	result <- do.call(rbind, c(lapply(result.list, function(x) x$data), make.row.names = F))
	samples <- do.call(rbind, c(lapply(result.list, function(x) x$samples), make.row.names = F))
	peaks <- do.call(rbind, c(lapply(result.list, function(x) x$peaks), make.row.names = F))
	alignment.values <- lapply(result.list, function(x) x$alignment.values)
	
	# convert sample metadata into factors, ensuring all frames have the same levels and the levels are in the observed order
	for (field in c("batch", "well.number", "sample.name", "sample.observations")) {
		samples[,field] <- factor(samples[,field], levels = unique(samples[,field]))
		if (field %in% names(result)) result[,field] <- factor(result[,field], levels = levels(samples[,field]))
		if (field %in% names(peaks)) peaks[,field] <- factor(peaks[,field], levels = levels(samples[,field]))
	}
	# convert other text into factors without those restrictions
	peaks[,"peak.observations"] <- factor(peaks[,"peak.observations"])
	
	# read smear regions
	# they are only defined once for the whole assay, but for compatibility with TapeStation data they must be defined repeatedly for each sample (TapeStation can have different regions for different samples)
	regions.raw <- xmlToDataFrame(chip.root[["AssayBody"]][["DASampleSetpoints"]][["DAMSmearAnalysis"]][["Channel"]][["RegionsMolecularSetpoints"]], stringsAsFactors = F)
	regions <- if (nrow(regions.raw) == 0) NULL else data.frame(samples[rep(1:nrow(samples), each = nrow(regions.raw)),], lower.length = as.numeric(regions.raw$StartBasePair), upper.length = as.numeric(regions.raw$EndBasePair), row.names = NULL)
	
	# analyze ladder
	which.ladder <- which(samples$is.ladder)
	stopifnot(length(which.ladder) == 1)
	peaks.ladder <- subset(peaks, well.number == samples$well.number[which.ladder])
	
	# fit standard curve for molecule length
	# do this with aligned times so it's effectively recalibrated for each sample's markers
	if (fit == "interpolate") {
		warning("linear interpolation gives ugly results for molarity estimation")
		standard.curve.function <- approxfun(peaks.ladder$aligned.time, peaks.ladder$length)
		standard.curve.inverse <- approxfun(peaks.ladder$length, peaks.ladder$aligned.time)
	} else if (fit == "spline") {
		standard.curve.function <- splinefun(peaks.ladder$aligned.time, peaks.ladder$length, method = "natural")
		standard.curve.inverse <- splinefun(peaks.ladder$length, peaks.ladder$aligned.time, method = "natural")
	} else if (fit == "regression") {
		mobility.model <- lm(1/aligned.time ~ log(length), data = peaks.ladder)
		standard.curve.function <- function(aligned.time) exp((1 / aligned.time - mobility.model$coefficients[1]) / mobility.model$coefficients[2])
		standard.curve.inverse <- function(length) 1/(mobility.model$coefficients[1] + log(length) * mobility.model$coefficients[2])
	}
	result$length <- standard.curve.function(result$aligned.time)
	ladder.limits <- range(defined.ladder.peaks$Size)
	result$length[result$length < ladder.limits[1] | result$length > ladder.limits[2]] <- NA # avoid extrapolating
	peaks$lower.length <- standard.curve.function(peaks$lower.aligned.time)
	peaks$upper.length <- standard.curve.function(peaks$upper.aligned.time)
	if (! is.null(regions)) {
		regions$lower.aligned.time <- standard.curve.inverse(regions$lower.length)
		regions$upper.aligned.time <- standard.curve.inverse(regions$upper.length)
		
		# convert inferred aligned times of regions back to raw times
		regions$lower.time <- NA
		regions$upper.time <- NA
		for (i in 1:nrow(samples)) {
			which.regions <- which(regions$well.number == samples$well.number[i])
			regions$lower.time[which.regions] <- (regions$lower.aligned.time[which.regions] - alignment.values[[i]][2])/alignment.values[[i]][1]
			regions$upper.time[which.regions] <- (regions$upper.aligned.time[which.regions] - alignment.values[[i]][2])/alignment.values[[i]][1]
		}
	}
	
	# annotate which peak each data point is in, if any
	# WARNING: if peaks overlap, this will overwrite and each point will only be mapped to the last-occuring one!
	result$peak <- NA
	for (i in 1:nrow(peaks)) result$peak[result$well.number == peaks$well.number[i] & result$aligned.time >= peaks$lower.aligned.time[i] & result$aligned.time <= peaks$upper.aligned.time[i]] <- i
	
	# convert to concentration and molarity
	# the idea is that we must correct fluorescence area by migration time to account for the fact that faster-moving molecules spend less time in front of the detector (Agilent's TimeCorrectedArea apparently does this with the raw time, not the aligned time)
	# and then fluorescence is proportional to concentration, which is molarity * length
	data.calibration <- cbind(result, do.call(rbind, lapply(samples$well.number, function(this.well) {
		result.this.well <- subset(result, well.number == this.well)
		data.frame(
			delta.fluorescence = c(NA, diff(result.this.well$fluorescence)),
			delta.time = c(NA, diff(result.this.well$time))
		)
	})))
	# estimate area under each measurement with the trapezoidal rule; to simplify math, each point's sum is for the trapezoid to the left of it
	data.calibration$area <- (2 * data.calibration$fluorescence - data.calibration$delta.fluorescence) * data.calibration$delta.time
	# correct area by migration time
	data.calibration$corrected.area <- data.calibration$area / data.calibration$time
	# fit the coefficient of mass vs. corrected area, using only the (non-marker) ladder peaks
	# this is because in the RNA kits, there is only one marker and its concentration is reported as zero, so we can't directly use it to calibrate the other samples
	# instead, we calculate one coefficient for the ladder and then scale each sample's coefficient by the area of its marker peak relative to the area of the ladder's marker peak
	peaks.calibration <- cbind(peaks, corrected.area = sapply(1:nrow(peaks), function(peak) sum(data.calibration$corrected.area[which(data.calibration$peak == peak)])))
	ladder.peaks <- subset(peaks.calibration, peak.observations == "Ladder Peak")
	ladder.mass.coefficient <- mean(ladder.peaks$concentration / ladder.peaks$corrected.area)
	marker.areas <- lapply(samples$well.number, function(well) {
		peaks.calibration$corrected.area[which(peaks$well.number == well & ((
			peaks$peak.observations == "Lower Marker" & peaks$concentration == defined.ladder.peaks$Concentration[1]
		) | (
			peaks$peak.observations == "Upper Marker" & peaks$concentration == defined.ladder.peaks$Concentration[nrow(defined.ladder.peaks)]
		)))]
	})
	mass.coefficients <- ladder.mass.coefficient * sapply(marker.areas, function(these.areas) mean(marker.areas[[which(samples$is.ladder)]] / these.areas)) # if there are two markers per sample, this gives the mean area ratio relative to their counterparts in the ladder well; if only one marker, the mean ratio is just the ratio 
	names(mass.coefficients) <- samples$well.number
	# apply this coefficient to get the concentration of each trapezoid
	result$concentration <- data.calibration$corrected.area * mass.coefficients[result$well.number]
	# finally scale by molecular weight to get the molarity
	result$molarity <- result$concentration / molecular.weight[[assay.info$assay.type]](result$length)
	
	# construct well.by.ladder (analogous to TapeStation but not as useful here)
	wells.by.ladder <- list(list(samples$well.number))
	names(wells.by.ladder) <- batch
	names(wells.by.ladder[[1]]) <- samples$well.number[which.ladder]
	
	# construct mobility.functions
	mobility.functions <- list(list(standard.curve.function))
	names(mobility.functions) <- batch
	names(mobility.functions[[1]]) <- samples$well.number[which.ladder]
	
	# construct mass.coefficients
	mass.coefficients <- list(mass.coefficients)
	names(mass.coefficients) <- batch
	
	# construct assay.info
	assay.info <- list(assay.info)
	names(assay.info) <- batch
	
	structure(list(
		data = result,
		assay.info = assay.info,
		samples = samples,
		wells.by.ladder = wells.by.ladder,
		peaks = peaks,
		regions = regions,
		mobility.functions = mobility.functions,
		mass.coefficients = mass.coefficients
	), class = "electrophoresis")
}

