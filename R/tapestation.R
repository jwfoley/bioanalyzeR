# hardcoded row number of aligned data CSVs
NROW.ALIGNED <- 760


#' Read a TapeStation CSV file
#'
#' This function reads a CSV file of electrophoresis data exported from the TapeStation software.
#'
#' Because the CSV file contains only the raw fluorescence data and not the metadata, this function is less useful by itself than when it is called inside \code{\link{read.tapestation}}.
#'
#' @param csv.file The filename of an electropherogram CSV file exported from the TapeStation software. The file may be compressed with \code{gzip} and the filename is expected to end in \code{_Electropherogram.csv} or \code{_Electropherogram.csv.gz}. The filename can be a remote URL. The electropherogram is expected to be unaligned.
#'
#' @return A data frame with one row for each reading of each gel lane.
#'
#' @seealso \code{\link{read.tapestation}}, \code{\link{read.tapestation.xml}}
#'
#' @export
read.tapestation.csv <- function(csv.file) {
	raw.data <- read.csv(csv.file, encoding = "latin1")
	if (nrow(raw.data) == NROW.ALIGNED) warning(paste("this looks like \"aligned\" data; did you export the unaligned data correctly?", csv.file))
	data.frame(
		sample.index = rep(seq(ncol(raw.data)), each = nrow(raw.data)),
		distance = rep(rev(seq(nrow(raw.data))), ncol(raw.data)) / nrow(raw.data),
		fluorescence = unlist(raw.data),
		row.names = seq(prod(dim(raw.data)))
	)
}


#' Read a TapeStation XML file
#'
#' This function reads an XML file exported from the TapeStation software.
#'
#' Because the XML file contains only metadata and not the raw fluorescence data, this function is less useful by itself than when it is called inside \code{\link{read.tapestation}}.
#'
#' @param xml.file The filename of an XML file exported from the TapeStation software. The file may be compressed with \code{gzip} and the filename is expected to end in \code{.xml} or \code{.xml.gz}; the name before that extension is used as the name of the batch. The filename can be a remote URL.
#'
#' @return A list of some of the components of an \code{\link{electrophoresis}} object
#' 
#' @seealso \code{\link{read.tapestation}}
#'
#' @export
#' @importFrom XML xmlRoot xmlParse xmlValue xmlApply xmlChildren xmlToDataFrame
read.tapestation.xml <- function(xml.file) {
 	batch <- sub("\\.xml(\\.gz)?$", "", basename(xml.file))
 	xml.root <- xmlRoot(xmlParse(xml.file))
 	
 	assay.info <- list(
 		file.name =           xmlValue(xml.root[["FileInformation"]][["FileName"]]),
 		creation.date =       xmlValue(xml.root[["FileInformation"]][["RunEndDate"]]),
 		assay.name =          xmlValue(xml.root[["FileInformation"]][["Assay"]]),
 		assay.type =          NULL,
 		length.unit =         NULL,
 		concentration.unit =  NULL,
 		molarity.unit =       NULL
 	)
 	# try to guess the assay type from the name
 	if (grepl("RNA", assay.info$assay.name)) {
 		assay.info$assay.type <- "RNA"
	} else if (grepl("D", assay.info$assay.name)) {
		assay.info$assay.type <- "DNA"
	} else {
		stop("unrecognized assay name")
	}
	assay.info$length.unit <- xmlValue(xml.root[["Assay"]][["Units"]][["MolecularWeightUnit"]])
	assay.info$concentration.unit <- xmlValue(xml.root[["Assay"]][["Units"]][["ConcentrationUnit"]])
	# hardcode the molarity unit depending on concentration unit, otherwise the MW scales will be wrong
	assay.info$molarity.unit <- switch(assay.info$concentration.unit,
		"ng/µl" = "nM",
		"pg/µl" = "pM"
	)
 	
	result.list <- xmlApply(xml.root[["Samples"]], function(sample.xml) {
		well.number <- xmlValue(sample.xml[["WellNumber"]])
		sample.name <- trimws(xmlValue(sample.xml[["Comment"]]))
		sample.observations <- trimws(xmlValue(sample.xml[["Observations"]]))
		if (sample.observations == "Marker(s) not detected") warning(paste(sample.observations, "for well", well.number, sample.name))
		is.ladder <- grepl("Ladder(?! run as sample)(?! sizing changed)", sample.observations, perl = T) # use negative lookaheads avoid false positives from known warning messages ("Comment" can be edited by user to anything so only the automatic "Observations" can be trusted to identify ladder, and they may contain various extraneous messages)
		if (sample.name == "") sample.name <- well.number
		well.row <- ifelse(is.ladder && sample.name == "Electronic Ladder", NA, substr(well.number, 1, 1))
		well.col <- ifelse(is.ladder && sample.name == "Electronic Ladder", NA, substr(well.number, 2, 3))
		suppressWarnings(RINe <- as.numeric(xmlValue(sample.xml[["RNA"]][["RINe"]])))
		suppressWarnings(DIN <- as.numeric(xmlValue(sample.xml[["DIN"]])))
		
		reagent.id <- xmlValue(sample.xml[["ScreenTapeID"]])
		
		suppressWarnings( # will throw warnings if missing values are coerced to NA but we can live with that
			peaks <- if (length(xmlChildren(sample.xml[["Peaks"]])) == 0) NULL else {
				peaks.raw <- xmlToDataFrame(sample.xml[["Peaks"]], stringsAsFactors = F)
				data.frame(
					peak.observations =  trimws(peaks.raw$Observations),
					peak.comment =       trimws(peaks.raw$Comment),
					length =             as.integer(peaks.raw$Size),
					distance =           as.numeric(peaks.raw$RunDistance) / 100,
					lower.distance =     as.numeric(peaks.raw$FromPercent) / 100,
					upper.distance =     as.numeric(peaks.raw$ToPercent) / 100,
					concentration =      as.numeric(peaks.raw$CalibratedQuantity),
					molarity =           as.numeric(peaks.raw$Molarity),
					stringsAsFactors =   F
				)
			}
		)
		
		suppressWarnings( # will throw warnings if missing values are coerced to NA but we can live with that
			regions <- if (length(xmlChildren(sample.xml[["Regions"]])) == 0) NULL else {
				regions.raw <- xmlToDataFrame(sample.xml[["Regions"]], stringsAsFactors = F)
				data.frame(
					region.comment =       trimws(regions.raw$Comment),
					lower.length =         as.integer(regions.raw$From),
					upper.length =         as.integer(regions.raw$To),
					average.length =       as.integer(regions.raw$AverageSize),
					concentration =        as.integer(regions.raw$Concentration),
					molarity =             as.numeric(regions.raw$Molarity),
					proportion.of.total =  as.numeric(regions.raw$PercentOfTotal) / 100,
					stringsAsFactors =     F
				)
			}
		)
		
		electrophoresis(
			samples = data.frame(batch, well.number, well.row, well.col, sample.name, sample.observations, reagent.id, RINe, DIN, is.ladder, stringsAsFactors = F),
			peaks = peaks,
			regions = regions
		)
	})
	has.peaks <- ! all(sapply(result.list, function(x) is.null(x$peaks)))
	has.regions <- ! all(sapply(result.list, function(x) is.null(x$regions)))
	result <- electrophoresis(
		samples = do.call(rbind, c(lapply(result.list, function(x) x$samples), make.row.names = F)),
		peaks = if (! has.peaks) NULL else do.call(rbind, c(lapply(seq_along(result.list), function(i) if (is.null(result.list[[i]]$peaks)) NULL else cbind(sample.index = i, result.list[[i]]$peaks)), make.row.names = F)),
		regions = if (! has.regions) NULL else do.call(rbind, c(lapply(seq_along(result.list), function(i) if (is.null(result.list[[i]]$regions)) NULL else cbind(sample.index = i, result.list[[i]]$regions)), make.row.names = F)),
		assay.info = assay.info
	)
	if (all(is.na(result$samples$RINe))) result$samples$RINe <- NULL
	if (all(is.na(result$samples$DIN))) result$samples$DIN <- NULL
	
	# convert sample metadata into factors, ensuring all frames have the same levels and the levels are in the observed order
	for (field in c("batch", "well.number", "sample.name", "reagent.id", "sample.observations")) result$samples[,field] <- factor(result$samples[,field], levels = unique(result$samples[,field]))
	# convert well row and column into factors but use the range of all possible rows/columns as levels
	result$samples$well.row <- factor(result$samples$well.row, levels = LETTERS[1:8])
	result$samples$well.col <- factor(result$samples$well.col, levels = 1:12)
	# convert other text into factors without those restrictions
	if (has.peaks) for (field in c("peak.observations", "peak.comment")) result$peaks[,field] <- factor(result$peaks[,field])
	if (has.regions) for (field in c("region.comment")) result$regions[,field] <- factor(result$regions[,field])
	
	# determine ladder scheme
	ladder.wells <- result$samples$well.number[which(result$samples$is.ladder)]
	result$samples$ladder.well <- factor(NA, levels = levels(result$samples$well.number))
	# missing ladder! leave NA
	if (length(ladder.wells) == 0) {
		return(result)
	# scheme: only one ladder for the whole run
	} else if (length(ladder.wells) == 1) {
		result$samples$ladder.well <- ladder.wells
	# scheme: one electronic ladder for the whole run but it's displayed more than once
	} else if (
		length(unique(result$samples$reagent.id[result$samples$well.number %in% ladder.wells])) == 1 &&
		"Electronic Ladder" %in% result$samples$sample.name[result$samples$well.number %in% ladder.wells] # not all copies get this name for some reason but at least one should
	) {
		result$samples$ladder.well <- ladder.wells[1]
	} else if (all.equal(as.character(result$samples$reagent.id[result$samples$well.number %in% ladder.wells]), levels(result$samples$reagent.id))) { 
		for (ladder.well in ladder.wells) result$samples$ladder.well[result$samples$reagent.id == result$samples$reagent.id[result$samples$well.number == ladder.well]] <- ladder.well
	}
	
	result
}


#' @describeIn read.electrophoresis Read a TapeStation XML and CSV file pair
#'
#' @inheritParams read.tapestation.xml
#' @inheritParams read.tapestation.csv
#' @inheritParams calibrate.electrophoresis
#'
#' @export
read.tapestation <- function(xml.file, csv.file = NULL, method = "hyman", extrapolate = FALSE) {
	# find the raw data file
	if (is.null(csv.file)) { # none provided
		# note the use of file.exists breaks the promise of using remote URLs
		candidate.files <- sapply(c(
			csv1 = "_Electropherogram.csv",
			csv2 = "_Electropherogram.csv.gz"
		), function(suffix) sub("\\.xml(\\.gz)?$", suffix, xml.file))
		found.files <- sapply(candidate.files, file.exists)
		stopifnot("found conflicting CSV files" = sum(found.files) <= 1)
		stopifnot("couldn't find CSV file" = sum(found.files) > 0)
		csv.file <- candidate.files[found.files]
	} 
	
	parsed.data <- read.tapestation.xml(xml.file)
	stopifnot("multiple batches provided" = length(unique(parsed.data$samples$batch)) == 1)
	batch <- parsed.data$samples$batch[1]
	result <- electrophoresis(
		data = read.tapestation.csv(csv.file),
		assay.info = setNames(list(parsed.data$assay.info), batch),
		samples = parsed.data$samples,
		peaks = parsed.data$peaks,
		regions = parsed.data$regions
	)
	
	# remove duplicate ladders
	if (sum(result$samples$is.ladder) > 1 && length(unique(result$samples$ladder.well)) == 1) {
		result <- subset(result, ! (is.ladder & well.number != ladder.well))
	}
	
	# calculate relative distances
	is.lower.marker <- result$peaks$peak.observations %in% LOWER.MARKER.NAMES
	marker.distances <- data.frame(lower = sapply(seq(nrow(result$samples)), function(i) {
		distance <- result$peaks$distance[is.lower.marker & result$peaks$sample.index == i]
		if (length(distance) == 0) {
			return(NA)
		} else if (length(distance) == 1) {
			return(distance)
		} else {
			stop(paste("multiple lower marker peaks for sample", i))
		}
	}))
	is.upper.marker <- result$peaks$peak.observations %in% UPPER.MARKER.NAMES
	if (sum(is.upper.marker) == 0) { # kit lacks upper marker
		marker.distances$upper <- 0 # effectively normalizes only to lower marker
	} else {
		marker.distances$upper <- sapply(seq(nrow(result$samples)), function(i) {
			distance <- result$peaks$distance[is.upper.marker & result$peaks$sample.index == i]
			if (length(distance) == 0) {
				return(NA)
			} else if (length(distance) == 1) {
				return(distance)
			} else {
				stop(paste("multiple upper marker peaks for sample", i))
			}
		})
	}
	marker.distances$range <- marker.distances$lower - marker.distances$upper
	result$data$relative.distance <- (result$data$distance - marker.distances$upper[result$data$sample.index]) / marker.distances$range[result$data$sample.index]
	result$peaks$relative.distance <- (result$peaks$distance - marker.distances$upper[result$peaks$sample.index]) / marker.distances$range[result$peaks$sample.index]
	result$peaks$lower.relative.distance <- (result$peaks$lower.distance - marker.distances$upper[result$peaks$sample.index]) / marker.distances$range[result$peaks$sample.index]
	result$peaks$upper.relative.distance <- (result$peaks$upper.distance - marker.distances$upper[result$peaks$sample.index]) / marker.distances$range[result$peaks$sample.index]
	
	# abort early if there aren't ladders to use
	if (sum(! is.na(result$samples$ladder.well)) == 0) {
		warning("no ladder identified so calibrations are not performed")
		return(result)
	}
	
	# perform calibrations
	result <- calculate.molarity(calculate.concentration(calculate.length(result, method, extrapolate)))
	
	# convert inferred relative distance of regions back to raw distance
	if (! is.null(result$regions)) {
		result$regions$lower.distance <- result$regions$lower.relative.distance * marker.distances$range[result$regions$sample.index] + marker.distances$upper[result$regions$sample.index]
		result$regions$upper.distance <- result$regions$upper.relative.distance * marker.distances$range[result$regions$sample.index] + marker.distances$upper[result$regions$sample.index]
	}
	
	result$samples$is.ladder <- NULL # remove extraneous temporary field
	
	result
}

