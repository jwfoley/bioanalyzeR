# filename suffixes
SUFFIX <- list(
	ELECTROPHEROGRAM = " Electropherogram.csv",
	PEAKS = " Peak Table.csv",
	REGIONS = " Smear Analysis Result.csv",
	CALIBRATION = " Size Calibration.csv"
)

# reformatting units
CONCENTRATION.UNITS <- list(
	`ng/ul` = "ng/µl",
	`ng/uL` = "ng/µl"
)
MOLARITY.UNITS <- list(
	`nmole/L` = "nM"
)

# guessing assay type from length column label
ASSAY.TYPE <- list(
	bp = "DNA",
	nt = "RNA"
)


#' Read a ProSize electropherogram
#'
#' This function reads an electropherogram (fluorescence vs. migration table) from the ProSize software saved in CSV format. The electropherogram must have been exported with migration times, not estimated sizes.
#'
#' @param csv.file The filename of an electropherogram CSV exported by ProSize. The filename can be a URL.
#'
#' @return A list containing a data frame of the raw fluorescence data and a data frame of the sample metadata (a partial \code{electrophoresis} object).
#'
#' @seealso \code{\link{read.prosize}}, \code{\link{read.prosize.peaks}}, \code{\link{read.prosize.regions}}
#'
#' @export
read.prosize.electropherogram <- function(csv.file) {
	data.raw <- read.csv(csv.file, check.names = F)
	stopifnot("electropherogram must be exported with times, not sizes" = colnames(data.raw)[1] == "Time (sec)")
	batch <- sub(paste0(SUFFIX$ELECTROPHEROGRAM, "$"), "", basename(csv.file))
	
	sample.long.names <- names(data.raw)[-1]
	n.samples <- length(sample.long.names)
	well.numbers <- sub(":.*", "", sample.long.names)
	sample.names <- sub("^[^:]+: *", "", sample.long.names)
	which.ladder <- grep("ladder", sample.names, ignore.case = T)
	if (length(which.ladder) > 1) {
		which.ladder <- which.ladder[length(which.ladder)]
		warning(paste("multiple ladders detected; only well", well.numbers[which.ladder], "used for calibration"))
	}
	
	list(
		data = data.frame(
			sample.index = rep(seq(n.samples), each = nrow(data.raw)),
			aligned.time = data.raw[,1],
			fluorescence = unlist(data.raw[,-1], use.names = F),
			row.names = NULL
		),
		samples = data.frame(
			batch,
			well.number = well.numbers,
			well.row = substr(well.numbers, 1, 1),
			well.col = substr(well.numbers, 2, 3),
			sample.name = sample.names,
			ladder.well = well.numbers[which.ladder]
		)
	)
}

#' Read a ProSize peak table
#'
#' This function reads a peak table from the ProSize software saved in CSV format. The peak table must have been exported in the "alternate" format and it must include the "From" and "To" columns.
#'
#' @param csv.file The filename of a peak table CSV exported by ProSize. The filename can be a URL.
#'
#' @return A list containing a list of assay metadata and a data frame of peaks (a partial \code{electrophoresis} object).
#'
#' @seealso \code{\link{read.prosize}}, \code{\link{read.prosize.electropherogram}}, \code{\link{read.prosize.regions}}
#'
#' @export
read.prosize.peaks <- function(csv.file) {
	peaks.raw <- read.csv(csv.file, check.names = F)
	batch <- sub(paste0(SUFFIX$PEAKS, "$"), "", basename(csv.file))
	
	# parse units and rename columns for easy reference later
	cols <- list(
		length = which(startsWith(colnames(peaks.raw), "Size (")),
		conc = which(colnames(peaks.raw) %in% names(CONCENTRATION.UNITS)),
		molarity = which(colnames(peaks.raw) %in% names(MOLARITY.UNITS)),
		percent = which(startsWith(colnames(peaks.raw), "% (Conc.) (")),
		lower.length = which(startsWith(colnames(peaks.raw), "From (")),
		upper.length = which(startsWith(colnames(peaks.raw), "To ("))
	)
	stopifnot(
		"missing or duplicated size column" = length(cols$length) == 1,
		"missing or duplicated concentration column" = length(cols$conc) == 1,
		"missing or duplicated molarity column" = length(cols$molarity) == 1,
		"missing or duplicated percent concentration column" = length(cols$percent) == 1,
		"missing or duplicated From column" = length(cols$lower.length) == 1,
		"missing or duplicated To column" = length(cols$upper.length) == 1
	)
	length.unit <- sub("\\)$", "", sub("^Size \\(", "", colnames(peaks.raw)[cols$length]))
	conc.unit <- CONCENTRATION.UNITS[[colnames(peaks.raw)[cols$conc]]]
	molarity.unit <- MOLARITY.UNITS[[colnames(peaks.raw)[cols$molarity]]]
	colnames(peaks.raw)[unlist(cols)] <- names(cols)
	
	# handle bad values
	peaks.raw <- subset(peaks.raw, ! is.na(conc))
	for (field in c("length", "lower.length", "upper.length")) if (class(peaks.raw[,field]) == "character") { # unparsed because of bad values
		peaks.raw[startsWith(peaks.raw[,field], ">"),field] <- NA
		peaks.raw[,field] <- as.integer(peaks.raw[,field])
	}
	
	# guess which peaks are markers: they don't have percent concentrations
	peak.observations <- rep("", nrow(peaks.raw))
	is.marker <- is.na(peaks.raw$percent)
	marker.lengths <- unique(peaks.raw$length[is.marker])
	peak.observations[is.marker & peaks.raw$length == marker.lengths[1]] <- "Lower Marker"
	if (length(marker.lengths) > 1) {
		stopifnot("too many markers" = length(marker.lengths) == 2)
		peak.observations[is.marker & peaks.raw$length == marker.lengths[2]] <- "Upper Marker"
	}
	
	list(
		assay.info = setNames(list(list(
			creation.date = batch,
			assay.type = ASSAY.TYPE[[length.unit]],
			length.unit = length.unit,
			concentration.unit = conc.unit,
			molarity.unit = molarity.unit
		)), batch),
		peaks = data.frame(
			batch,
			well.number = sub(":", "", peaks.raw$Well), # for empty sample name the colon may be in the well number (ProSize bug)
			sample.name = peaks.raw$`Sample ID`,
			peak.observations,
			length = peaks.raw$length,
			lower.length = peaks.raw$lower.length,
			upper.length = peaks.raw$upper.length,
			concentration = peaks.raw$conc,
			molarity = peaks.raw$molarity
		)
	)
}

#' Read a ProSize smear analysis
#'
#' This function reads smear analysis table from the ProSize software saved in CSV format. The smear analysis must have been exported in the "alternate" format.
#'
#' @param csv.file The filename of a smear analysis CSV exported by ProSize. The filename can be a URL.
#'
#' @return A list containing a list of assay metadata and a data frame of regions (a partial \code{electrophoresis} object).
#'
#' @seealso \code{\link{read.prosize}}, \code{\link{read.prosize.electropherogram}}, \code{\link{read.prosize.peaks}}
#'
#' @export
read.prosize.regions <- function(csv.file) {
	regions.raw <- subset(read.csv(csv.file, check.names = F), Range != "") # regions aren't reported for ladder well but they get a partially empty line
	batch <- sub(paste0(SUFFIX$REGIONS, "$"), "", basename(csv.file))
	
	# parse units and rename columns for easy reference later
	cols <- list(
		conc = which(colnames(regions.raw) %in% names(CONCENTRATION.UNITS)),
		molarity = which(colnames(regions.raw) %in% names(MOLARITY.UNITS))
	)
	stopifnot(
		"missing or duplicated concentration column" = length(cols$conc) == 1,
		"missing or duplicated molarity column" = length(cols$molarity) == 1
	)
	conc.unit <- CONCENTRATION.UNITS[[colnames(regions.raw)[cols$conc]]]
	molarity.unit <- MOLARITY.UNITS[[colnames(regions.raw)[cols$molarity]]]
	colnames(regions.raw)[unlist(cols)] <- names(cols)
	length.unit <- unique(sub(".* ", "", regions.raw$Range))
	stopifnot("conflicting units detected" = length(length.unit) == 1)
	
	list(
		assay.info = setNames(list(list(
			creation.date = batch,
			assay.type = ASSAY.TYPE[[length.unit]],
			length.unit = length.unit,
			concentration.unit = conc.unit,
			molarity.unit = molarity.unit
		)), batch),
		regions = data.frame(
			batch,
			well.number = sub(":", "", regions.raw$Well), # for empty sample name the colon may be in the well number (ProSize bug)
			sample.name = regions.raw$`Sample ID`,
			lower.length = as.integer(sub(" bp to .*$", "", regions.raw$Range)),
			upper.length = as.integer(sub(" bp$", "", sub("^.* bp to ", "", regions.raw$Range))),
			concentration = regions.raw$conc,
			molarity = regions.raw$molarity
		)
	)
}

#' @describeIn read.electrophoresis Read ProSize CSV files
#'
#' @export
read.prosize <- function(csv.file, fit = "spline") {
	root.path <- sub(paste0(SUFFIX$ELECTROPHEROGRAM, "$"), "", csv.file)
	batch <- basename(root.path)
	
	data.import <- read.prosize.electropherogram(paste0(root.path, SUFFIX$ELECTROPHEROGRAM))
	master.samples <- split(data.import$samples[,c("batch", "well.number", "sample.name")], seq(nrow(data.import$samples))) # table of sample metadata to verify other files have at least a subset of these
	stopifnot(
		"conflicting batch names in electropherogram" = length(unique(data.import$samples$batch)) == 1,
		"duplicate well numbers in electropherogram" = anyDuplicated(data.import$samples$well.number) == 0
	)
	sample.index.lookup <- setNames(as.list(seq(nrow(data.import$samples))), data.import$samples$well.number)
	
	# add peaks safely
	peaks.import <- read.prosize.peaks(paste0(root.path, SUFFIX$PEAKS))
	stopifnot("conflicting sample metadata in peak table" = all(split(peaks.import$peaks[,c("batch", "well.number", "sample.name")], seq(nrow(peaks.import$peaks))) %in% master.samples))
	peaks <- data.frame(sample.index = unlist(sample.index.lookup[peaks.import$peaks$well.number]), peaks.import$peaks[,! colnames(peaks.import$peaks) %in% c("batch", "well.number", "sample.name")])
	
	# add regions safely
	regions.csv <- paste0(root.path, SUFFIX$REGIONS)
	regions <- if (file.exists(regions.csv)) {
		regions.import <- read.prosize.regions(regions.csv)
		stopifnot(
			"conflicting sample metadata in smear analysis table" = all(split(regions.import$regions[,c("batch", "well.number", "sample.name")], seq(nrow(regions.import$regions))) %in% master.samples),
			"conflicting batch metadata between smear analysis and peaks" = all.equal(regions.import$assay.info, peaks.import$assay.info)
		)
		data.frame(sample.index = unlist(sample.index.lookup[regions.import$regions$well.number]), regions.import$regions[,! colnames(regions.import$regions) %in% c("batch", "well.number", "sample.name")])
	} else NULL
	
	# build the electrophoresis object	
	result <- structure(list(
		data = data.import$data,
		assay.info = peaks.import$assay.info,
		samples = data.import$samples,
		peaks = peaks,
		regions = regions
	), class = "electrophoresis")
	
	# convert sample metadata into factors, ensuring all frames have the same levels and the levels are in the observed order
	for (field in c("batch", "well.number", "sample.name")) result$samples[,field] <- factor(result$samples[,field], levels = unique(result$samples[,field]))
	result$samples$ladder.well <- factor(result$samples$ladder.well, levels = levels(result$samples$well.number))
	# convert well row and column into factors but use the range of all possible rows/columns as levels
	result$samples$well.row <- factor(result$samples$well.row, levels = LETTERS[1:8])
	result$samples$well.col <- factor(result$samples$well.col, levels = 1:12)
	# convert other text into factors without those restrictions
	result$peaks$peak.observations <- factor(result$peaks$peak.observations)
	
	# read ladder calibration and reverse ProSize's calibration to get peak times
	ladder.peaks <- read.csv(paste0(root.path, SUFFIX$CALIBRATION), check.names = F)
	which.ladder <- which(result$samples$well.number == result$samples$ladder.well)
	stopifnot("multiple ladders" = length(which.ladder) == 1) # assume only one ladder
	which.ladder.peaks <- result$peaks$sample.index == which.ladder
	stopifnot("conflicting ladder peaks in calibration" = all(ladder.peaks[,1] == result$peaks$length[which.ladder.peaks]))
	reverse.calibration <- approxfun(ladder.peaks[,1], ladder.peaks[,2], rule = 2) # extrapolates values outside the marker peaks to be the same! so basically each marker peak's coordinates are missing its outer part (lower marker's lower.aligned.time = aligned.time, upper marker's upper.aligned.time = aligned.time)
	result$peaks$aligned.time <- reverse.calibration(result$peaks$length)
	result$peaks$lower.aligned.time <- reverse.calibration(result$peaks$lower.length)
	result$peaks$upper.aligned.time <- reverse.calibration(result$peaks$upper.length)
	
	calculate.molarity(calculate.concentration(calculate.length(result, fit)))
}

