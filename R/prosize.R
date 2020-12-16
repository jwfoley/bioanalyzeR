# filename suffixes
SUFFIX <- list(
	ELECTROPHEROGRAM = " Electropherogram.csv",
	PEAKS = " Peak Table.csv",
	REGIONS = " Smear Analysis Result.csv",
	CALIBRATION = " Size Calibration.csv"
)

# reformatting units
UNIT <- list(
	`nmole/L` = "nM",
	`ng/ul` = "ng/µl",
	`ng/uL` = "ng/µl"
)

# guessing assay type from length column label
ASSAY.TYPE <- list(
	bp = "DNA",
	nt = "RNA"
)

# columns in peak table
PEAKS.COL <- list(
	WELL.NUMBER = 1,
	SAMPLE.NAME = 2,
	LENGTH = 4,
	CONC.PERCENT = 5,
	MOLARITY = 6,
	CONC = 7,
	AREA = 8
)

# columns in region table
REGIONS.COL <- list(
	WELL.NUMBER = 1,
	SAMPLE.NAME = 2,
	RANGE = 3,
	CONC = 4,
	PERCENT.TOTAL = 5,
	MOLARITY = 6,
	AVG.LENGTH = 7
)

#' @export
read.prosize.electropherogram <- function(csv.file) {
	data.raw <- read.csv(csv.file, check.names = F)
	stopifnot("electropherogram must be exported with times, not sizes" = colnames(data.raw)[1] == "Time (sec)")
	batch <- sub(paste0(SUFFIX$ELECTROPHEROGRAM, "$"), "", basename(csv.file))
	
	sample.long.names <- names(data.raw)[-1]
	n.samples <- length(sample.long.names)
	well.numbers <- sub(":.*", "", sample.long.names)
	sample.names <- sub("^[^:]+: ", "", sample.long.names)
	which.ladder <- which(sample.names == "Ladder")
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
			sample.name = sample.names,
			ladder.well = well.numbers[which.ladder]
		)
	)
}

#' @export
read.prosize.peaks <- function(csv.file) {
	peaks.raw <- read.csv(csv.file, check.names = F)
	batch <- sub(paste0(SUFFIX$PEAKS, "$"), "", basename(csv.file))
	length.unit <- sub("\\)$", "", sub("^Size \\(", "", colnames(peaks.raw)[PEAKS.COL$LENGTH]))
	
	# guess which peaks are markers: they don't have percent concentrations
	peak.observations <- rep("", nrow(peaks.raw))
	is.marker <- is.na(peaks.raw[,PEAKS.COL$CONC.PERCENT])
	marker.lengths <- unique(peaks.raw[is.marker,PEAKS.COL$LENGTH])
	# this assumes there are always two markers: is that true?
	stopifnot("wrong number of markers detected" = length(marker.lengths) == 2)
	peak.observations[is.marker & peaks.raw[,PEAKS.COL$LENGTH] == marker.lengths[1]] <- "Lower Marker"
	peak.observations[is.marker & peaks.raw[,PEAKS.COL$LENGTH] == marker.lengths[2]] <- "Upper Marker"
	
	list(
		assay.info = setNames(list(list(
			creation.date = batch,
			assay.type = ASSAY.TYPE[[length.unit]],
			length.unit = length.unit,
			concentration.unit = UNIT[[colnames(peaks.raw)[PEAKS.COL$CONC]]],
			molarity.unit = UNIT[[colnames(peaks.raw)[PEAKS.COL$MOLARITY]]]
		)), batch),
		peaks = data.frame(
			batch,
			well.number = peaks.raw[,PEAKS.COL$WELL.NUMBER],
			sample.name = peaks.raw[,PEAKS.COL$SAMPLE.NAME],
			peak.observations,
			length = peaks.raw[,PEAKS.COL$LENGTH],
			aligned.time = NA,
			lower.aligned.time = NA,
			upper.aligned.time = NA,
			area = peaks.raw[,PEAKS.COL$AREA],
			concentration = peaks.raw[,PEAKS.COL$CONC],
			molarity = peaks.raw[,PEAKS.COL$MOLARITY]
		)
	)
}

#' @export
read.prosize.regions <- function(csv.file) {
	regions.raw <- subset(read.csv(csv.file, check.names = F), Range != "") # regions aren't reported for ladder well but they get a partially empty line
	batch <- sub(paste0(SUFFIX$REGIONS, "$"), "", basename(csv.file))
	length.unit <- unique(sub(".* ", "", regions.raw[,REGIONS.COL$RANGE]))
	stopifnot("conflicting units detected" = length(length.unit) == 1)
	
	list(
		assay.info = setNames(list(list(
			creation.date = batch,
			assay.type = ASSAY.TYPE[[length.unit]],
			length.unit = length.unit,
			concentration.unit = UNIT[[colnames(regions.raw)[REGIONS.COL$CONC]]],
			molarity.unit = UNIT[[colnames(regions.raw)[REGIONS.COL$MOLARITY]]]
		)), batch),
		regions = data.frame(
			batch,
			well.number = regions.raw[,REGIONS.COL$WELL.NUMBER],
			sample.name = regions.raw[,REGIONS.COL$SAMPLE.NAME],
			lower.length = as.integer(sub(" bp to .*$", "", regions.raw[,REGIONS.COL$RANGE])),
			upper.length = as.integer(sub(" bp$", "", sub("^.* bp to ", "", regions.raw[,REGIONS.COL$RANGE]))),
			average.length = regions.raw[,REGIONS.COL$AVG.LENGTH],
			concentration = regions.raw[,REGIONS.COL$CONC],
			molarity = regions.raw[,REGIONS.COL$MOLARITY],
			proportion.of.total = regions.raw[,REGIONS.COL$PERCENT.TOTAL] / 100
		)
	)
}

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
	
	# add ladder calibration
	ladder.peaks <- read.csv(paste0(root.path, SUFFIX$CALIBRATION), check.names = F)
	which.ladder <- which(result$samples$well.number == result$samples$ladder.well)
	stopifnot("multiple ladders" = length(which.ladder) == 1) # assume only one ladder
	which.ladder.peaks <- result$peaks$sample.index == which.ladder
	stopifnot("conflicting ladder peaks in calibration" = all(ladder.peaks$`Ladder Size (bp)` == result$peaks$length[which.ladder.peaks]))
	result$peaks$aligned.time[which.ladder.peaks] <- ladder.peaks$`Time (sec)`
	
	calculate.molarity(calculate.concentration(calculate.length(result, fit)))
}

