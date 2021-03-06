# hardcoded colors
RGB.HIGHLIGHT <-  c(209, 228, 250)  # highlight around selected lane is light blue
RGB.GOOD <-       c(138, 208, 160)  # label for high RIN is green
RGB.MEDIUM <-     c(255, 237, 101)  # label for medium RIN is yellow
RGB.BAD <-        c(255, 106,  71)  # label for low RIN is red
RGB.EMPTY <-      c(255, 255, 255)  # image gets pasted into all-white background, probably

# hardcoded margin widths, in pixels
WARNING.PAD <- 10 # extra pixels to discard on both sides of a warning label in a gel lane


# find all pixels in an RGB array from readPNG, with values in [0,1], that match a given RGB trio, with values in [0, 255]
find.matching.pixels.mat <- function(rgb.mat, rgb.values) {
	rgb.fractions <- rgb.values / 255
	return(rgb.mat[,,1] == rgb.fractions[1] & rgb.mat[,,2] == rgb.fractions[2] & rgb.mat[,,3] == rgb.fractions[3])
}

# version for a "vector" (actually a matrix with one spatial dimension and one RGB channel dimension)
find.matching.pixels.vec <- function(rgb.vec, rgb.values) {
	rgb.fractions <- rgb.values / 255
	return(rgb.vec[,1] == rgb.fractions[1] & rgb.vec[,2] == rgb.fractions[2] & rgb.vec[,3] == rgb.fractions[3])
}


#' Read a TapeStation gel image
#'
#' (DEPRECATED) This function reads a gel image exported from the TapeStation software and saved in PNG format. The gel image must include a blue highlight around one lane in order for the function to identify the boundaries of the gel area.
#'
#' DEPRECATED: Reading raw data from gel images is deprecated and produces noisy results. Please reopen your data in version 4.1 or higher of the TapeStation Analysis software and export a CSV instead: \url{https://explore.agilent.com/Software-Download-TapeStation-Systems}
#'
#' Because the gel image alone contains little metadata, this function returns only a simple data frame containing the fluorescence intensity vs. migration distance at every point in every lane of the gel (numbered from left to right). It is less useful by itself than when it is called inside \code{\link{read.tapestation}}.
#'
#' Note: Fluorescence is reported only for one column of pixels down the center of each lane. If there are annotations on the gel image, such as the yellow warning symbol, those pixels are reported as NA, but usually they are at the top of the image outside the reportable range anyway.
#'
#' @param gel.image.file The filename of a TapeStation gel image with blue highlight, in PNG format. The filename can be a URL.
#' @param n.lanes The number of lanes in the gel image.
#'
#' @return A data frame with one row for each vertical pixel of each gel lane.
#' 
#' @seealso \code{\link{read.tapestation}}, \code{\link{read.tapestation.xml}}, \code{\link{read.tapestation.csv}}
#'
#' @export
#' @importFrom png readPNG
read.tapestation.gel.image <- function(gel.image.file, n.lanes) {
	warning("Reading gel images is deprecated; export CSV instead")
	
	gel.image.con <- file(gel.image.file, "rb", raw = T) # workaround to allow URLs 
	gel.image.rgb <- readPNG(readBin(gel.image.con, what = "raw", n = 25E6)) # safe overestimate of the maximum size (slightly above 4K resolution @ 24 bits uncompressed)
	close(gel.image.con) # if not explicitly closed, R gives a warning
	# note: this is in the form (y, x, channel); [1,1,] is the upper left corner
	
	# find gel boundaries
	highlight.cols <- which(find.matching.pixels.vec(gel.image.rgb[1,,], RGB.HIGHLIGHT)) # use only the first pixel row to find the highlight
	stopifnot("highlight is not continuous in first pixel row" = all(diff(highlight.cols) == 1))
	highlighted.lane.width <- length(highlight.cols)
	highlighted.subset <- gel.image.rgb[,highlight.cols,]
	highlight.rows <- which(rowSums(find.matching.pixels.mat(highlighted.subset, RGB.HIGHLIGHT)) == highlighted.lane.width) # find all pixel rows with full highlight (will miss ones with annotation text over them)
	top.highlight.rows <- highlight.rows[highlight.rows < nrow(gel.image.rgb) / 2]
	end.of.top.highlight <- top.highlight.rows[length(top.highlight.rows)] # assume it's the last row in the top half
	subposition.is.quality.label <- find.matching.pixels.mat(highlighted.subset, RGB.GOOD) | find.matching.pixels.mat(highlighted.subset, RGB.MEDIUM) | find.matching.pixels.mat(highlighted.subset, RGB.BAD)
	start.of.bottom.highlight <- if (any(subposition.is.quality.label)) which(rowSums(subposition.is.quality.label) > 0)[1] else highlight.rows[length(top.highlight.rows) + 1] # quality label supersedes any blue highlight
	highlight.border.offsets <- which(find.matching.pixels.vec(highlighted.subset[end.of.top.highlight + 1,,], RGB.HIGHLIGHT)) # sometimes will be empty if there's only one pixel of border and the color is off because of antialiasing, but we can live with that much error
	stopifnot(
		"border is not contiguous" = length(highlight.border.offsets < 2) || all(diff(highlight.border.offsets) == 1),
		"border is not on edge" = length(highlight.border.offsets) == 0 || (highlight.border.offsets[1] %in% 0:1 || highlight.border.offsets[length(highlight.border.offsets)] %in% (highlighted.lane.width - 1:0)) # assume the border is on one edge or the other, allowing one pixel column of error due to antialiasing
	)
	lane.center <- ((if (1 %in% highlight.border.offsets) highlight.border.offsets[length(highlight.border.offsets)] else 0) +(highlighted.lane.width - length(highlight.border.offsets)) / 2) / highlighted.lane.width # approximate x-position of the center of the lane, from the left, as a proportion of the total width
	
	# identify empty left and right margins
	margin.columns <- colSums((! find.matching.pixels.mat(gel.image.rgb, RGB.EMPTY))) == 0
	left.margin <- if (head(margin.columns, 1)) head(which(diff(margin.columns) == -1), 1) else 0
	right.margin <- if (tail(margin.columns, 1)) ncol(gel.image.rgb) - tail(which(diff(margin.columns) == 1), 1) else 0
	
	# extract fluorescence values by lane
	average.lane.width <- (ncol(gel.image.rgb) - left.margin - right.margin) / n.lanes
	lane.pixels <- gel.image.rgb[
		(start.of.bottom.highlight - 1):(end.of.top.highlight + 1), # reverse rows to put fastest migration first like Bioanalyzer
		left.margin + 1 + round(average.lane.width * (seq(n.lanes) - 1 + lane.center)),
	]
	n.readings <- nrow(lane.pixels)
	bad.pixels <- ! (lane.pixels[,,1] == lane.pixels[,,2] & lane.pixels[,,1] == lane.pixels[,,3]) # find non-grayscale pixels, indicating annotations that block the data
	for (col in which(colSums(bad.pixels) > 0)) {
		bad.rows <- which(bad.pixels[,col])
		lane.pixels[max((min(bad.rows) - WARNING.PAD), 1):min((max(bad.rows) + WARNING.PAD), n.readings), col,] <- NA # set bad pixels and pad around them to NA (don't let pad go outside the range)
	}
	
	data.frame(
		sample.index =  rep(seq(n.lanes), each = n.readings),
		distance =      n.readings:1 / n.readings,
		fluorescence =  as.vector(1 - lane.pixels[,,1]) # subtract from 1 because it's a negative; use only red channel
	)
}


#' Read a TapeStation CSV file
#'
#' This function reads a CSV file of electrophoresis data exported from the TapeStation software.
#'
#' Because the CSV file contains only the raw fluorescence data and not the metadata, this function is less useful by itself than when it is called inside \code{\link{read.tapestation}}.
#'
#' @param csv.file The filename of a CSV file exported from the TapeStation software. The file may be compressed with \code{gzip} and the filename is expected to end in \code{_Electropherogram.csv} or \code{_Electropherogram.csv.gz}. The filename can be a remote URL.
#'
#' @return A data frame with one row for each reading of each gel lane.
#'
#' @seealso \code{\link{read.tapestation}}, \code{\link{read.tapestation.xml}}, \code{\link{read.tapestation.gel.image}}
#'
#' @export
read.tapestation.csv <- function(csv.file) {
	raw.data <- read.csv(csv.file, encoding = "latin1")
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
#' @seealso \code{\link{read.tapestation}}, \code{\link{read.tapestation.gel.image}}
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
		is.ladder <- grepl("Ladder", sample.observations) & ! grepl("Ladder run as sample", sample.observations)
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


#' @describeIn read.electrophoresis Read a TapeStation XML and PNG file pair
#'
#' @inheritParams read.tapestation.xml
#' @inheritParams read.tapestation.csv
#' @inheritParams read.tapestation.gel.image
#' @inheritParams calibrate.electrophoresis
#'
#' @export
read.tapestation <- function(xml.file, csv.file = NULL, gel.image.file = NULL, method = "hyman") {
	# find the raw data file
	if (is.null(csv.file) && is.null(gel.image.file)) { # none provided
		# note the use of file.exists breaks the promise of using remote URLs; removing PNG support may resolve this
		candidate.files <- sapply(c(
			gel.image = ".png",
			csv1 = "_Electropherogram.csv",
			csv2 = "_Electropherogram.csv.gz"
		), function(suffix) sub("\\.xml(\\.gz)?$", suffix, xml.file))
		found.files <- sapply(candidate.files, file.exists)
		stopifnot("found conflicting CSV or PNG files" = sum(found.files) <= 1)
		stopifnot("couldn't find CSV or PNG file" = sum(found.files) > 0)
		if (found.files[1]) gel.image.file <- candidate.files[1] else csv.file <- candidate.files[found.files]
	} 
	stopifnot("conflicting CSV and PNG files provided" = is.null(csv.file) || is.null(gel.image.file))
	
	parsed.data <- read.tapestation.xml(xml.file)
	stopifnot("multiple batches provided" = length(unique(parsed.data$samples$batch)) == 1)
	batch <- parsed.data$samples$batch[1]
	result <- electrophoresis(
		data = if (is.null(gel.image.file)) read.tapestation.csv(csv.file) else read.tapestation.gel.image(gel.image.file, nrow(parsed.data$samples)),
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
	result <- calculate.molarity(calculate.concentration(calculate.length(result, method)))
	
	# convert inferred relative distance of regions back to raw distance
	if (! is.null(result$regions)) {
		result$regions$lower.distance <- result$regions$lower.relative.distance * marker.distances$range[result$regions$sample.index] + marker.distances$upper[result$regions$sample.index]
		result$regions$upper.distance <- result$regions$upper.relative.distance * marker.distances$range[result$regions$sample.index] + marker.distances$upper[result$regions$sample.index]
	}
	
	result$samples$is.ladder <- NULL # remove extraneous temporary field
	
	result
}

