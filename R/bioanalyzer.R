#' @describeIn read.electrophoresis Read a Bioanalyzer XML file
#'
#' @inheritParams calibrate.electrophoresis
#'
#' @export
#' @importFrom XML xmlRoot xmlParse xmlValue xmlToDataFrame xmlApply
#' @importFrom base64enc base64decode
read.bioanalyzer <- function(xml.file, method = "hyman") {
	batch <- sub("\\.xml(\\.gz)?$", "", basename(xml.file))
	xml.root <- xmlRoot(xmlParse(xml.file))
	chip.root <- xml.root[["Chips"]][["Chip"]]
	assay.info <- list(
		file.name =           xmlValue(chip.root[["Files"]][["File"]][["FileInformation"]][["FileName"]]),
		creation.date =       xmlValue(chip.root[["ChipInformation"]][["CreationDate"]]),
		assay.name =          xmlValue(chip.root[["AssayHeader"]][["Title"]]),
		assay.type =          xmlValue(chip.root[["AssayHeader"]][["Class"]]),
		length.unit =         xmlValue(chip.root[["AssayBody"]][["DAAssaySetpoints"]][["DAMAssayInfoMolecular"]][["SizeUnit"]]),
		concentration.unit =  xmlValue(chip.root[["AssayBody"]][["DAAssaySetpoints"]][["DAMAssayInfoMolecular"]][["ConcentrationUnit"]]),
		molarity.unit =       NULL
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
			suppressWarnings(RIN <- as.numeric(xmlValue(this.sample[["DAResultStructures"]][["DARRIN"]][["Channel"]][["RIN"]])))

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
				time = as.numeric(xmlValue(signal.data[["XStart"]])) + as.numeric(xmlValue(signal.data[["XStep"]])) * (seq(n.values) - 1),
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
				samples = data.frame(batch, well.number, sample.name, sample.observations, sample.comment, RIN, is.ladder, stringsAsFactors = F),
				peaks = peaks,
				alignment.values = c(alignment.coefficient, alignment.offset)
			)
		}
	})
	for (i in length(result.list):1) if (is.null(result.list[[i]])) result.list[[i]] <- NULL # delete empty elements; go in reverse because indexes will change as elements are deleted
	result <- structure(list(
		data = do.call(rbind, c(lapply(seq_along(result.list), function(i) cbind(sample.index = i, result.list[[i]]$data)), make.row.names = F)),
		assay.info = setNames(list(assay.info), batch),
		samples = do.call(rbind, c(lapply(result.list, function(x) x$samples), make.row.names = F)),
		peaks = do.call(rbind, c(lapply(seq_along(result.list), function(i) if (is.null(result.list[[i]]$peaks)) NULL else cbind(sample.index = i, result.list[[i]]$peaks)), make.row.names = F)),
		regions = NULL,
		mobility.functions = NULL
	), class = "electrophoresis")
	if (all(is.na(result$samples$RIN))) result$samples$RIN <- NULL
	alignment.values <- lapply(result.list, function(x) x$alignment.values)
	
	# convert sample metadata into factors, ensuring all frames have the same levels and the levels are in the observed order
	for (field in c("batch", "well.number", "sample.name", "sample.observations", "sample.comment")) result$samples[[field]] <- factor(result$samples[[field]], levels = unique(result$samples[[field]]))
	
	# convert other text into factors without those restrictions
	result$peaks[,"peak.observations"] <- factor(result$peaks[,"peak.observations"])
	
	# read smear regions
	# they are only defined once for the whole assay, but for compatibility with TapeStation data they must be defined repeatedly for each sample (TapeStation can have different regions for different samples)
	regions.raw <- xmlToDataFrame(chip.root[["AssayBody"]][["DASampleSetpoints"]][["DAMSmearAnalysis"]][["Channel"]][["RegionsMolecularSetpoints"]], stringsAsFactors = F)
	if (nrow(regions.raw) > 0) result$regions <- data.frame(
		sample.index = rep(seq(nrow(result$samples)), each = nrow(regions.raw)),
		lower.length = as.numeric(regions.raw$StartBasePair),
		upper.length = as.numeric(regions.raw$EndBasePair),
		row.names = NULL
	)
	
	# analyze ladder
	which.ladder <- result$samples$well.number[result$samples$is.ladder]
	stopifnot(length(which.ladder) == 1)
	result$samples$is.ladder <- NULL
	result$samples$ladder.well <- factor(which.ladder, levels = levels(result$samples$well.number))
	
	# perform calibrations
	result <- calculate.molarity(calculate.concentration(calculate.length(result, method), defined.ladder.peaks$Concentration))
	
	# convert inferred aligned times of regions back to raw times
	if (! is.null(result$regions)) {
		result$regions$lower.time <- NA
		result$regions$upper.time <- NA
		for (i in seq(nrow(result$samples))) {
			which.regions <- which(result$regions$sample.index == i)
			result$regions$lower.time[which.regions] <- (result$regions$lower.aligned.time[which.regions] - alignment.values[[i]][2])/alignment.values[[i]][1]
			result$regions$upper.time[which.regions] <- (result$regions$upper.aligned.time[which.regions] - alignment.values[[i]][2])/alignment.values[[i]][1]
		}	
	}
	
	result
}

