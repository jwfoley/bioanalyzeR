# hardcoded colors
RGB.UPPER.MARKER <-  c(128,   0, 128)  # upper marker is purple
RGB.LOWER.MARKER <-  c(  0, 128,   0)  # lower marker is green
RGB.HIGHLIGHT <-     c(209, 228, 250)  # highlight around selected lane is light blue
RGB.GOOD <-          c(138, 208, 160)  # label for high RIN is green
RGB.MEDIUM <-        c(255, 237, 101)  # label for medium RIN is yellow
RGB.BAD <-           c(255, 106,  71)  # label for low RIN is red

# hardcoded margin widths, in pixels
LEFT.MARGIN <- 8
RIGHT.MARGIN <- 15


# find all pixels in an RGB array from readPNG, with values in [0,1], that match a given RGB trio, with values in [0, 255]
find.matching.pixels <- function(rgb.image, rgb.values) {
	rgb.fractions <- rgb.values / 255
	return(rgb.image[,,1] == rgb.fractions[1] & rgb.image[,,2] == rgb.fractions[2] & rgb.image[,,3] == rgb.fractions[3])
}


#' Read a TapeStation gel image
#'
#' This function reads a gel image exported from the TapeStation software and saved in PNG format. The gel image must include a blue highlight around one lane in order for the function to identify the boundaries of the gel area.
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
#' @seealso \code{\link{read.tapestation}}, \code{\link{read.tapestation.xml}}
#'
#' @export
#' @importFrom png readPNG
read.tapestation.gel.image <- function(gel.image.file, n.lanes) {
	gel.image.con <- file(gel.image.file, "rb", raw = T) # workaround to allow URLs 
	gel.image.rgb <- readPNG(readBin(gel.image.con, what = "raw", n = 25E6)) # safe overestimate of the maximum size (slightly above 4K resolution @ 24 bits uncompressed)
	close(gel.image.con) # if not explicitly closed, R gives a warning
	# note: this is in the form (y, x, channel); [1,1,] is the upper left corner
	
	# stop if markers are overlaid
	if (any(find.matching.pixels(gel.image.rgb, RGB.LOWER.MARKER) | find.matching.pixels(gel.image.rgb, RGB.UPPER.MARKER))) stop('Marker bands detected in gel image. Please re-export it after unchecking "Show Marker Annotations".')
	
	# find gel boundaries
	position.is.highlight <- find.matching.pixels(gel.image.rgb, RGB.HIGHLIGHT)
	highlight.cols <- which(position.is.highlight[1,]) # use only the first pixel row to find the highlight
	stopifnot(all(diff(highlight.cols) == 1)) # assume the highlight is continuous in the first pixel row
	highlighted.lane.width <- length(highlight.cols)
	highlight.subset <- position.is.highlight[,highlight.cols]
	highlight.rows <- which(rowSums(highlight.subset) == highlighted.lane.width) # find all pixel rows with full highlight (will miss ones with annotation text over them)
	top.highlight.rows <- highlight.rows[highlight.rows < nrow(gel.image.rgb) / 2]
	end.of.top.highlight <- top.highlight.rows[length(top.highlight.rows)] # assume it's the last row in the top half
	highlighted.positions <- gel.image.rgb[,highlight.cols,]
	subposition.is.quality.label <- find.matching.pixels(highlighted.positions, RGB.GOOD) | find.matching.pixels(highlighted.positions, RGB.MEDIUM) | find.matching.pixels(highlighted.positions, RGB.BAD)
	start.of.bottom.highlight <- if (any(subposition.is.quality.label)) which(rowSums(subposition.is.quality.label) > 0)[1] else highlight.rows[length(top.highlight.rows) + 1] # quality label supersedes any blue highlight
	highlight.border.offsets <- which(highlight.subset[end.of.top.highlight + 1,]) # sometimes will be empty if there's only one pixel of border and the color is off because of antialiasing, but we can live with that much error
	stopifnot(length(highlight.border.offsets < 2) || all(diff(highlight.border.offsets) == 1)) # assume the border is contiguous
	stopifnot(length(highlight.border.offsets) == 0 || (highlight.border.offsets[1] == 1 || highlight.border.offsets[length(highlight.border.offsets)] == highlighted.lane.width)) # assume the border is on one edge or the other
	lane.center <- ((if (1 %in% highlight.border.offsets) highlight.border.offsets[length(highlight.border.offsets)] else 0) +(highlighted.lane.width - length(highlight.border.offsets)) / 2) / highlighted.lane.width # approximate x-position of the center of the lane, from the left, as a proportion of the total width
	
	# extract fluorescence values by lane
	average.lane.width <- (ncol(gel.image.rgb) - LEFT.MARGIN - RIGHT.MARGIN) / n.lanes
	fluorescence.matrix <- 1 - gel.image.rgb[ # subtract from 1 because it's negative
		(start.of.bottom.highlight - 1):(end.of.top.highlight + 1), # reverse so it will be in order from fastest to slowest migration like Bioanalyzer
		LEFT.MARGIN + 1 + round(average.lane.width * (1:n.lanes - 1 + lane.center)),
		1 # use only the red channel
	]
	
	data.frame(
		sample.index =  rep(1:ncol(fluorescence.matrix), each = nrow(fluorescence.matrix)),
		distance =      nrow(fluorescence.matrix):1 / nrow(fluorescence.matrix),
		fluorescence =  as.vector(fluorescence.matrix)
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
#' @return A list of some of the components of an \code{electrophoresis} object
#' 
#' @seealso \code{\link{read.tapestation}}, \code{\link{read.tapestation.gel.image}}
#'
#' @export
#' @importFrom XML xmlRoot xmlParse xmlValue xmlApply xmlChildren xmlToDataFrame
read.tapestation.xml <- function(xml.file) {
 	batch <- sub("\\.xml(\\.gz)?$", "", basename(xml.file))
 	xml.root <- xmlRoot(xmlParse(xml.file))
 	
 	assay.info <- list(
 		file.name =      xmlValue(xml.root[["FileInformation"]][["FileName"]]),
 		creation.date =  xmlValue(xml.root[["FileInformation"]][["RunEndDate"]]),
 		assay.name =     xmlValue(xml.root[["FileInformation"]][["Assay"]])
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
		if (sample.name == "") sample.name <- well.number
		sample.observations <- trimws(xmlValue(sample.xml[["Observations"]]))
		if (sample.observations == "Marker(s) not detected") {
			warning(paste(sample.observations, "for well", well.number, sample.name))
			return(NULL)
		}
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
		
		list(
			samples = data.frame(batch, well.number, sample.name, sample.observations, reagent.id, RINe, DIN, stringsAsFactors = F),
			peaks = peaks,
			regions = regions
		)
	})
	has.peaks <- ! all(sapply(result.list, function(x) is.null(x$peaks)))
	has.regions <- ! all(sapply(result.list, function(x) is.null(x$regions)))
	result <- list(
		samples = do.call(rbind, c(lapply(result.list, function(x) x$samples), make.row.names = F)),
		peaks = if (! has.peaks) NULL else do.call(rbind, c(lapply(1:length(result.list), function(i) if (is.null(result.list[[i]]$peaks)) NULL else cbind(sample.index = i, result.list[[i]]$peaks)), make.row.names = F)),
		regions = if (! has.regions) NULL else do.call(rbind, c(lapply(1:length(result.list), function(i) if (is.null(result.list[[i]]$regions)) NULL else cbind(sample.index = i, result.list[[i]]$regions)), make.row.names = F)),
		assay.info = assay.info
	)
	if (all(is.na(result$samples$RINe))) result$samples$RINe <- NULL
	if (all(is.na(result$samples$DIN))) result$samples$DIN <- NULL
	
	# convert sample metadata into factors, ensuring all frames have the same levels and the levels are in the observed order
	for (field in c("batch", "well.number", "sample.name", "reagent.id", "sample.observations")) result$samples[,field] <- factor(result$samples[,field], levels = unique(result$samples[,field]))
	# convert other text into factors without those restrictions
	if (has.peaks) for (field in c("peak.observations", "peak.comment")) result$peaks[,field] <- factor(result$peaks[,field])
	if (has.regions) for (field in c("region.comment")) result$regions[,field] <- factor(result$regions[,field])
	
	# determine ladder scheme
	ladder.wells <- result$samples$well.number[result$samples$sample.observations == "Ladder"]
	result$samples$ladder.well <- factor(NA, levels = levels(result$samples$well.number))
	# scheme: only one ladder for the whole run
	if (length(ladder.wells) == 1) {
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
#' @export
read.tapestation <- function(xml.file, gel.image.file = NULL, fit = "spline") {
	stopifnot(fit %in% c("interpolation", "spline", "regression"))
	if (is.null(gel.image.file)) gel.image.file <- sub("\\.xml(\\.gz)?$", ".png", xml.file)
	
	parsed.data <- read.tapestation.xml(xml.file)
	stopifnot(length(unique(parsed.data$samples$batch)) == 1)
	batch <- parsed.data$samples$batch[1]
	result <- structure(list(
		data = read.tapestation.gel.image(gel.image.file, nrow(parsed.data$samples)),
		assay.info = list(parsed.data$assay.info),
		samples = parsed.data$samples,
		peaks = parsed.data$peaks,
		regions = parsed.data$regions,
		mobility.functions = NULL,
		mass.coefficients = NULL
	), class = "electrophoresis")
	names(result$assay.info) <- batch
	
	# remove duplicate ladders
	if (sum(result$samples$sample.observations == "Ladder") > 1 && length(unique(result$samples$ladder.well)) == 1) {
		result <- subset(result, ! (sample.observations == "Ladder" & well.number != ladder.well))
	}
	
	# calculate relative distances
	is.lower.marker <- result$peaks$peak.observations %in% c("Lower Marker", "edited Lower Marker")
	marker.distances <- data.frame(lower = sapply(1:nrow(result$samples), function(i) {
		distance <- result$peaks$distance[is.lower.marker & result$peaks$sample.index == i]
		if (length(distance) == 0) {
			return(NA)
		} else if (length(distance) == 1) {
			return(distance)
		} else {
			stop(paste("multiple lower marker peaks for sample", i))
		}
	}))
	is.upper.marker <- result$peaks$peak.observations %in% c("Upper Marker", "edited Upper Marker")
	if (sum(is.upper.marker) == 0) { # kit lacks upper marker
		marker.distances$upper <- 0 # effectively normalizes only to lower marker
	} else {
		marker.distances$upper <- sapply(1:nrow(result$samples), function(i) {
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
		warning("unknown ladder scheme so lengths and molarities are not calculated")
		return(result)
	}
	
	# prepare more fields to be filled in piecemeal from each ladder
	result$data$length <- NA
	result$data$concentration <- NA
	result$data$molarity <- NA
	result$peaks$lower.length <- NA
	result$peaks$upper.length <- NA
	if (! is.null(result$regions)) {
		result$regions$lower.relative.distance <- NA
		result$regions$upper.relative.distance <- NA
	}
	result$mobility.functions <- list(list())
	names(result$mobility.functions) <- batch
	result$mass.coefficients <- rep(NA, nrow(result$samples))
	data.calibration <- cbind(result$data, do.call(rbind, lapply(1:nrow(result$samples), function(i) {
		which.this.sample <- result$data$sample.index == i
		data.frame(
			delta.fluorescence = c(NA, diff(result$data$fluorescence[which.this.sample])),
			delta.distance = c(NA, -diff(result$data$distance[which.this.sample]))
		)
	})))
	# estimate area under each measurement with the trapezoidal rule; to simplify math, each point's sum is for the trapezoid to the left of it
	data.calibration$area <- (2 * data.calibration$fluorescence - data.calibration$delta.fluorescence) * data.calibration$delta.distance
	for (ladder.well in unique(result$samples$ladder.well)) {
		which.ladder.index <- which(result$samples$well.number == ladder.well)
		peaks.ladder <- subset(result$peaks, sample.index == which.ladder.index)
		which.samples <- which(result$samples$ladder.well == ladder.well)
		which.rows <- which(result$data$sample.index %in% which.samples)
		which.peaks <- which(result$peaks$sample.index %in% which.samples)
		which.regions <- which(result$regions$sample.index %in% which.samples)
		
		# fit standard curve for molecule length vs. migration distance
		# do this in relative-distance space so it's effectively recalibrated for each sample's markers
		if (fit == "interpolation") {
			warning("linear interpolation gives ugly results for molarity estimation")
			standard.curve.function <- approxfun(peaks.ladder$relative.distance, peaks.ladder$length)
			standard.curve.inverse <- approxfun(peaks.ladder$length, peaks.ladder$relative.distance)
		} else if (fit == "spline") {
			standard.curve.function <- splinefun(peaks.ladder$relative.distance, peaks.ladder$length, method = "natural")
			standard.curve.inverse <- splinefun(peaks.ladder$length, peaks.ladder$relative.distance, method = "natural")
		} else if (fit == "regression") {
			mobility.model <- lm(relative.distance ~ log(length), peaks.ladder)
			standard.curve.function <- function(relative.distance) exp((relative.distance - mobility.model$coefficients[1]) / mobility.model$coefficients[2])
			standard.curve.inverse <- function(length) mobility.model$coefficients[1] + mobility.model$coefficients[2] * log(length)
		}
		result$mobility.functions[[1]][[ladder.well]] <- standard.curve.function
		
		# apply model to raw data
		result$data$length[which.rows] <- standard.curve.function(result$data$relative.distance[which.rows])
		result$data$length[! in.custom.region(result$data, min(peaks.ladder$length), max(peaks.ladder$length))] <- NA # avoid extrapolation
		
		# apply model to peaks
		result$peaks$lower.length[which.peaks] <- standard.curve.function(result$peaks$upper.relative.distance[which.peaks])
		result$peaks$upper.length[which.peaks] <- standard.curve.function(result$peaks$lower.relative.distance[which.peaks])
		
		# apply inverse model to regions
		if (! is.null(result$regions)) {
			result$regions$lower.relative.distance[which.regions] <- standard.curve.inverse(result$regions$upper.length[which.regions])
			result$regions$upper.relative.distance[which.regions] <- standard.curve.inverse(result$regions$lower.length[which.regions])
		}
		
		# convert to molarity
		# first solve for a coefficient that relates known masses (known molarities scaled by length) to area under the electropherogram peaks, then apply that coefficient to each individual measurement
		# this is done by fitting a one-parameter model on the non-marker peaks of the ladder, because the marker peaks in all samples are unreadable (blocked by green and purple bands)
		# preferably it would be scaled to each sample's own markers, but that's impossible because the markers are unreadable! and we can't even use their Agilent-reported molarity estimates or areas to scale relative to the ladder because those are always normalized to the upper marker! (which is often the one that's contaminated by sample anyway)
		# so all we can do is normalize to the ladder's non-marker peaks, therefore each sample will be randomly off by some constant scaling factor, but at least molarity comparisons within a sample ought to be accurate
		peaks.ladder$area <- sapply(1:nrow(peaks.ladder), function(i) sum(subset(data.calibration, sample.index == which.ladder.index & distance >= peaks.ladder$lower.distance[i] & distance <= peaks.ladder$upper.distance[i])$area))
		mass.coefficient <- lm(concentration ~ area - 1, data = peaks.ladder)$coefficients[1]
		result$mass.coefficients[which.samples] <- mass.coefficient
		result$data$concentration[which.rows] <- mass.coefficient * data.calibration$area[which.rows]
		result$data$molarity[which.rows] <- result$data$concentration[which.rows] / molecular.weight(result$data$length[which.rows], parsed.data$assay.info$assay.type) * 1E6 # we're converting ng/uL to nmol/L or pg/uL to pmol/L so we need to scale by 1E6
	}
	
	# convert inferred relative distances of regions back to raw distances
	if (! is.null(result$regions)) {
		result$regions$lower.distance <- result$regions$lower.relative.distance * marker.distances$range[result$regions$sample.index] + marker.distances$upper[result$regions$sample.index]
		result$regions$upper.distance <- result$regions$upper.relative.distance * marker.distances$range[result$regions$sample.index] + marker.distances$upper[result$regions$sample.index]
	}
	
	result
}

