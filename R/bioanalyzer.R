#' Read a Bioanalyzer XML file
#'
#' This function reads an XML file exported from the Bioanalyzer software, then fills out the results with estimates of molecule length, concentration, and molarity.
#'
#' Spline fitting seems to perform reasonably well on all data. Agilent appears to use linear interpolation with DNA data and log-linear regression on RNA data, so you could choose those options if you want to reproduce the results of the software more precisely. However, linear interpolation creates sudden spikes in the derivative that make the concentration and molarity estimates unstable; spline fitting is basically a smoother version of that. Log-linear regression is the standard theoretical approach but does not actually fit the data very well; more sophisticated parametric models may be added in the future.
#'
#' @param xml.file The filename of an XML file exported from the Bioanalyzer software. The filename is expected to end in \code{.xml} and the name before that extension is used as the name of the batch.
#' @param fit The method used to fit the mobility model of molecule length vs. migration distance, one of \code{"interpolation"} (linear interpolation via \code{\link{approxfun}}), \code{"spline"} (splines via \code{\link{splinefun}}), or \code{"regression"} (log-linear regression via \code{\link{lm}} with the model \code{relative.distance ~ log(length)}).
#'
#' @return An \code{electrophoresis} object containing the data from this Bioanalyzer run.
#' 
#' @seealso \code{\link{read.electrophoresis}}, \code{\link{read.tapestation}}
#'
#' @export
#' @importFrom XML xmlRoot xmlParse xmlValue xmlToDataFrame xmlApply
#' @importFrom base64enc base64decode
read.bioanalyzer <- function(xml.file, fit = "spline") {
	stopifnot(fit %in% c("interpolation", "spline", "regression"))
	
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
			sample.comment <- trimws(xmlValue(this.sample[["Comment"]]))
			sample.observations <- trimws(xmlValue(this.sample[["ResultFlagCommonLabel"]]))

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
				fluorescence = readBin(base64decode(xmlValue(signal.data[["ProcessedSignal"]])), "numeric", size = 4, n = n.values, endian = "little"),
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
				data = raw.data,
				samples = data.frame(batch, well.number, sample.name, sample.observations, sample.comment, is.ladder, stringsAsFactors = F),
				peaks = peaks,
				alignment.values = c(alignment.coefficient, alignment.offset)
			)
		}
	})
	result <- structure(list(
		data = do.call(rbind, c(lapply(1:length(result.list), function(i) cbind(sample.index = i, result.list[[i]]$data)), make.row.names = F)),
		assay.info = NULL,
		samples = do.call(rbind, c(lapply(result.list, function(x) x$samples), make.row.names = F)),
		peaks = do.call(rbind, c(lapply(1:length(result.list), function(i) if (is.null(result.list[[i]]$peaks)) NULL else cbind(sample.index = i, result.list[[i]]$peaks)), make.row.names = F)),
		regions = NULL,
		mobility.functions = NULL,
		mass.coefficients = NULL
	), class = "electrophoresis")
	alignment.values <- lapply(result.list, function(x) x$alignment.values)
		
	# construct assay.info
	result$assay.info <- list(assay.info)
	names(result$assay.info) <- batch
	
	# convert sample metadata into factors, ensuring all frames have the same levels and the levels are in the observed order
	for (field in c("batch", "well.number", "sample.name", "sample.observations", "sample.comment")) result$samples[,field] <- factor(result$samples[,field], levels = unique(result$samples[,field]))
	
	# convert other text into factors without those restrictions
	result$peaks[,"peak.observations"] <- factor(result$peaks[,"peak.observations"])
	
	# read smear regions
	# they are only defined once for the whole assay, but for compatibility with TapeStation data they must be defined repeatedly for each sample (TapeStation can have different regions for different samples)
	regions.raw <- xmlToDataFrame(chip.root[["AssayBody"]][["DASampleSetpoints"]][["DAMSmearAnalysis"]][["Channel"]][["RegionsMolecularSetpoints"]], stringsAsFactors = F)
	if (nrow(regions.raw) > 0) result$regions <- data.frame(sample.index = rep(1:nrow(result$samples), each = nrow(regions.raw)), lower.length = as.numeric(regions.raw$StartBasePair), upper.length = as.numeric(regions.raw$EndBasePair), row.names = NULL)
	
	# analyze ladder
	which.ladder <- which(result$samples$is.ladder)
	stopifnot(length(which.ladder) == 1)
	result$samples <- result$samples[,! names(result$samples) == "is.ladder"]
	peaks.ladder <- subset(result$peaks, sample.index == which.ladder & peak.observations %in% c("Lower Marker", "Ladder Peak", "Upper Marker"))
	result$samples$ladder.well <- factor(which.ladder, levels = levels(result$samples$well.number))
	
	# fit standard curve for molecule length
	# do this with aligned times so it's effectively recalibrated for each sample's markers
	if (fit == "interpolation") {
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
	result$mobility.functions <- list(list(standard.curve.function))
	names(result$mobility.functions[[1]]) <- which.ladder
	names(result$mobility.functions) <- batch
	
	# apply mobility function
	result$data$length <- standard.curve.function(result$data$aligned.time)
	result$data$length[! in.custom.region(result$data, min(defined.ladder.peaks$Size), max(defined.ladder.peaks$Size))] <- NA # avoid extrapolating
	result$peaks$lower.length <- standard.curve.function(result$peaks$lower.aligned.time)
	result$peaks$upper.length <- standard.curve.function(result$peaks$upper.aligned.time)
	if (! is.null(result$regions)) {
		result$regions$lower.aligned.time <- standard.curve.inverse(result$regions$lower.length)
		result$regions$upper.aligned.time <- standard.curve.inverse(result$regions$upper.length)
		
		# convert inferred aligned times of regions back to raw times
		result$regions$lower.time <- NA
		result$regions$upper.time <- NA
		for (i in 1:nrow(result$samples)) {
			which.regions <- which(result$regions$sample.index == i)
			result$regions$lower.time[which.regions] <- (result$regions$lower.aligned.time[which.regions] - alignment.values[[i]][2])/alignment.values[[i]][1]
			result$regions$upper.time[which.regions] <- (result$regions$upper.aligned.time[which.regions] - alignment.values[[i]][2])/alignment.values[[i]][1]
		}
	}
	
	# convert to concentration and molarity
	# the idea is that we must correct fluorescence area by migration time to account for the fact that faster-moving molecules spend less time in front of the detector (Agilent's TimeCorrectedArea apparently does this with the raw time, not the aligned time)
	# and then fluorescence is proportional to concentration, which is molarity * length
	data.calibration <- cbind(result$data, peak = in.peaks(result), do.call(rbind, lapply(1:nrow(result$samples), function(i) {
		result.this.well <- subset(result$data, sample.index == i)
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
	peaks.calibration <- cbind(result$peaks, corrected.area = sapply(1:nrow(result$peaks), function(peak) sum(data.calibration$corrected.area[which(data.calibration$peak == peak)])))
	ladder.peaks <- subset(peaks.calibration, peak.observations == "Ladder Peak")
	ladder.mass.coefficient <- mean(ladder.peaks$concentration / ladder.peaks$corrected.area)
	marker.areas <- lapply(1:nrow(result$samples), function(i) {
		peaks.calibration$corrected.area[which(result$peaks$sample.index == i & ((
			result$peaks$peak.observations == "Lower Marker" & result$peaks$concentration == defined.ladder.peaks$Concentration[1]
		) | (
			result$peaks$peak.observations == "Upper Marker" & result$peaks$concentration == defined.ladder.peaks$Concentration[nrow(defined.ladder.peaks)]
		)))]
	})
	result$mass.coefficients <- ladder.mass.coefficient * sapply(marker.areas, function(these.areas) mean(marker.areas[[which.ladder]] / these.areas)) # if there are two markers per sample, this gives the mean area ratio relative to their counterparts in the ladder well; if only one marker, the mean ratio is just the ratio 
	# apply this coefficient to get the concentration of each trapezoid
	result$data$concentration <- data.calibration$corrected.area * result$mass.coefficients[result$data$sample.index]
	# finally scale by molecular weight to get the molarity
	result$data$molarity <- result$data$concentration / molecular.weight(result$data$length, assay.info$assay.type) * 1E6 # we're converting ng/uL to nmol/L or pg/uL to pmol/L so we need to scale by 1E6
	
	result
}

