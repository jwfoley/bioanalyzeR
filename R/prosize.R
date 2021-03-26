# filename suffixes
SUFFIX <- list(
	ELECTROPHEROGRAM = " Electropherogram.csv",
	PEAKS = " Peak Table.csv",
	REGIONS = " Smear Analysis Result.csv",
	CALIBRATION = " Size Calibration.csv",
	QUALITY = " Quality Table.csv"
)

# reformatting units
CONCENTRATION.UNITS <- list(
	`ng/ul` = "ng/µl",
	`ng/uL` = "ng/µl",
	`pg/ul` = "pg/µl",
	`pg/uL` = "pg/µl"
)
MOLARITY.UNITS <- list(
	`nmole/L` = "nM",
	`pmole/L` = "pM"
)

# guessing assay type from length column label
ASSAY.TYPE <- list(
	bp = "DNA",
	nt = "RNA"
)

QUALITY.METRICS <- c(
	"GQN",
	"RQN",
	"28S/18S",
	"% RNA Contamination",
	"% microRNA"
)


#' Read a ProSize electropherogram
#'
#' This function reads an electropherogram (fluorescence vs. migration table) from the ProSize software saved in CSV format. The electropherogram must have been exported with migration times, not estimated sizes.
#'
#' @param electropherogram.csv An electropherogram CSV exported by ProSize.
#'
#' @return A list containing a data frame of the raw fluorescence data and a data frame of the sample metadata (a partial \code{\link{electrophoresis}} object).
#'
#' @seealso \code{\link{read.prosize}}, \code{\link{read.prosize.peaks}}, \code{\link{read.prosize.regions}}
#'
#' @export
read.prosize.electropherogram <- function(electropherogram.csv) {
	data.raw <- read.csv(electropherogram.csv, check.names = F)
	stopifnot("electropherogram must be exported with times, not sizes" = colnames(data.raw)[1] == "Time (sec)")
	
	sample.long.names <- names(data.raw)[-1]
	n.samples <- length(sample.long.names)
	well.numbers <- sub(":.*", "", sample.long.names)
	sample.names <- sub("^[^:]+: *", "", sample.long.names)
	
	electrophoresis(
		data = data.frame(
			sample.index = rep(seq(n.samples), each = nrow(data.raw)),
			aligned.time = data.raw[,1],
			fluorescence = unlist(data.raw[,-1], use.names = F),
			row.names = NULL
		),
		samples = data.frame(
			well.number = well.numbers,
			well.row = substr(well.numbers, 1, 1),
			well.col = substr(well.numbers, 2, 3),
			sample.name = sample.names
		)
	)
}

#' Read a ProSize peak table
#'
#' This function reads a peak table from the ProSize software saved in CSV format. The peak table must have been exported in the "alternate" format and it must include the "From" and "To" columns.
#'
#' @param peaks.csv A peak table CSV exported by ProSize.
#'
#' @return A list containing a list of assay metadata and a data frame of peaks (a partial \code{\link{electrophoresis}} object).
#'
#' @seealso \code{\link{read.prosize}}, \code{\link{read.prosize.electropherogram}}, \code{\link{read.prosize.regions}}
#'
#' @export
read.prosize.peaks <- function(peaks.csv) {
	peaks.raw <- read.csv(peaks.csv, check.names = F)
	
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
	
	electrophoresis(
		assay.info = list(
			assay.type = ASSAY.TYPE[[length.unit]],
			length.unit = length.unit,
			concentration.unit = conc.unit,
			molarity.unit = molarity.unit
		),
		peaks = data.frame(
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
#' @param smear.csv A smear analysis CSV exported by ProSize.
#'
#' @return A list containing a list of assay metadata and a data frame of regions (a partial \code{\link{electrophoresis}} object).
#'
#' @seealso \code{\link{read.prosize}}, \code{\link{read.prosize.electropherogram}}, \code{\link{read.prosize.peaks}}
#'
#' @export
read.prosize.regions <- function(smear.csv) {
	regions.raw <- subset(read.csv(smear.csv, check.names = F), Range != "") # regions aren't reported for ladder well but they get a partially empty line
	
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
	
	electrophoresis(
		assay.info = list(
			assay.type = ASSAY.TYPE[[length.unit]],
			length.unit = length.unit,
			concentration.unit = conc.unit,
			molarity.unit = molarity.unit
		),
		regions = data.frame(
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
#' @inheritParams read.prosize.electropherogram
#' @inheritParams read.prosize.peaks
#' @inheritParams read.prosize.regions
#' @param calibration.csv A size calibration CSV exported by ProSize.
#' @param quality.csv A quality table CSV exported by ProSize.
#' @param batch Name of the batch.
#'
#' @export
read.prosize <- function(
	electropherogram.csv,
	calibration.csv = NULL,
	peaks.csv = NULL,
	smear.csv = NULL,
	quality.csv = NULL,
	method = "hyman",
	batch = NULL
) {
	if ("character" %in% class(electropherogram.csv)) { # character or subclass
		root.path <- sub(paste0(SUFFIX$ELECTROPHEROGRAM, "$"), "", electropherogram.csv)
		if (is.null(calibration.csv)) calibration.csv <- paste0(root.path, SUFFIX$CALIBRATION)
		if (is.null(peaks.csv)) peaks.csv <- paste0(root.path, SUFFIX$PEAKS)
		if (is.null(smear.csv)) {
			smear.csv <- paste0(root.path, SUFFIX$REGIONS)
			if (! file.exists(smear.csv)) smear.csv <- NULL
		}
		if (is.null(quality.csv)) {
			quality.csv <- paste0(root.path, SUFFIX$QUALITY)
			if (! file.exists(quality.csv)) quality.csv <- NULL
		}
		if (is.null(batch)) batch <- basename(root.path)
	}
	
	data.import <- read.prosize.electropherogram(electropherogram.csv)
	master.samples <- split(data.import$samples[,c("well.number", "sample.name")], seq(nrow(data.import$samples))) # table of sample metadata to verify other files have at least a subset of these
	stopifnot("duplicate well numbers in electropherogram" = anyDuplicated(data.import$samples$well.number) == 0)
	sample.index.lookup <- setNames(as.list(seq(nrow(data.import$samples))), data.import$samples$well.number)
	
	# add peaks safely
	peaks.import <- read.prosize.peaks(peaks.csv)
	stopifnot("conflicting sample metadata in peak table" = all(split(peaks.import$peaks[,c("well.number", "sample.name")], seq(nrow(peaks.import$peaks))) %in% master.samples))
	peaks <- data.frame(sample.index = unlist(sample.index.lookup[peaks.import$peaks$well.number]), peaks.import$peaks[,! colnames(peaks.import$peaks) %in% c("well.number", "sample.name")])
	
	# add regions safely
	regions <- if (! is.null(smear.csv)) {
		regions.import <- read.prosize.regions(smear.csv)
		stopifnot(
			"conflicting sample metadata in smear analysis table" = all(split(regions.import$regions[,c("well.number", "sample.name")], seq(nrow(regions.import$regions))) %in% master.samples),
			"conflicting batch metadata between smear analysis and peaks" = all.equal(regions.import$assay.info, peaks.import$assay.info)
		)
		data.frame(sample.index = unlist(sample.index.lookup[regions.import$regions$well.number]), regions.import$regions[,! colnames(regions.import$regions) %in% c("well.number", "sample.name")])
	} else NULL
	
	# build the electrophoresis object	
	result <- electrophoresis(
		data = data.import$data,
		assay.info = setNames(list(peaks.import$assay.info), batch),
		samples = cbind(batch, data.import$samples),
		peaks = peaks,
		regions = regions
	)
	
	# add quality safely
	if (! is.null(quality.csv)) {
		quality.raw <- read.csv(quality.csv, check.names = F)
		quality.raw$Well <- sub(":", "", quality.raw$Well) # for empty sample name the colon may be in the well number (ProSize bug)
		stopifnot("conflicting sample metadata in quality table" =
			all(data.import$samples$well.number == quality.raw$Well) &&
			all(data.import$samples$sample.name == quality.raw$`Sample ID`)
		)
		for (field in intersect(names(quality.raw), QUALITY.METRICS)) result$samples[,field] <- quality.raw[,field]
	}
	
	# read ladder calibration and reverse ProSize's calibration to get peak times
	ladder.peaks <- read.csv(calibration.csv, check.names = F)
	which.ladder <- which(by(result$peaks, result$peaks$sample.index, function(sample.peaks) nrow(sample.peaks) == nrow(ladder.peaks) && all(sample.peaks$length == ladder.peaks[,1])))
	if (length(which.ladder) > 1) {
		warning(paste0("multiple ladders found in wells ", cat(result$samples$well.number[which.ladder]), "; using only ", result$samples$well.number[which.ladder[length(which.ladder)]]))
		which.ladder <- which.ladder[length(which.ladder)]
	}
	stopifnot("no ladder found" = length(which.ladder) == 1)
	which.ladder.peaks <- result$peaks$sample.index == which.ladder
	reverse.calibration <- approxfun(ladder.peaks[,1], ladder.peaks[,2], rule = 2) # extrapolates values outside the marker peaks to be the same! so basically each marker peak's coordinates are missing its outer part (lower marker's lower.aligned.time = aligned.time, upper marker's upper.aligned.time = aligned.time)
	result$peaks$aligned.time <- reverse.calibration(result$peaks$length)
	result$peaks$lower.aligned.time <- reverse.calibration(result$peaks$lower.length)
	result$peaks$upper.aligned.time <- reverse.calibration(result$peaks$upper.length)
	
	# convert sample metadata into factors, ensuring all frames have the same levels and the levels are in the observed order
	for (field in c("batch", "well.number", "sample.name")) result$samples[,field] <- factor(result$samples[,field], levels = unique(result$samples[,field]))
	result$samples$ladder.well <- result$samples$well.number[which.ladder]
	# convert well row and column into factors but use the range of all possible rows/columns as levels
	result$samples$well.row <- factor(result$samples$well.row, levels = LETTERS[1:8])
	result$samples$well.col <- factor(result$samples$well.col, levels = 1:12)
	# convert other text into factors without those restrictions
	result$peaks$peak.observations <- factor(result$peaks$peak.observations)
	
	calculate.molarity(calculate.concentration(calculate.length(result, method)))
}

#' @describeIn read.electrophoresis Read a ProSize ZIP file
#'
#' @param prosize.zip Path of a ZIP file containing the required CSVs exported by ProSize.
#'
#' @export
read.prosize.zip <- function(prosize.zip, method = "hyman") {
	all.files <- unzip(prosize.zip, list = T)$Name
	electropherogram.csv.name <- grep(SUFFIX$ELECTROPHEROGRAM, all.files, value = T)
	stopifnot("missing or duplicated electropherogram" = length(electropherogram.csv.name) == 1)
	root.name <- sub(paste0(SUFFIX$ELECTROPHEROGRAM, "$"), "", electropherogram.csv.name)
	filenames <- list(
		electropherogram = electropherogram.csv.name,
		calibration = paste0(root.name, SUFFIX$CALIBRATION),
		peaks = paste0(root.name, SUFFIX$PEAKS),
		regions = paste0(root.name, SUFFIX$REGIONS),
		quality = paste0(root.name, SUFFIX$QUALITY)
	)
	read.prosize(
		unz(prosize.zip, filenames$electropherogram),
		unz(prosize.zip, filenames$calibration),
		unz(prosize.zip, filenames$peaks),
		if (filenames$regions %in% all.files) unz(prosize.zip, filenames$regions) else NULL,
		if (filenames$quality %in% all.files) unz(prosize.zip, filenames$quality) else NULL,
		method,
		sub("\\.zip$", "", basename(prosize.zip))
	)
}

