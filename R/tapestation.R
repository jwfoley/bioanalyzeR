# hardcoded colors
RGB.UPPER.MARKER <-  c(128,   0, 128)  # upper marker is purple
RGB.LOWER.MARKER <-  c(  0, 128,   0)  # lower marker is green
RGB.HIGHLIGHT <-     c(209, 228, 250)  # highlight around selected lane is light blue
RGB.GOOD <-          c(138, 208, 160)  # label for high RIN is green
RGB.MEDIUM <-        c(255, 237, 101)  # label for medium RIN is yellow
RGB.BAD <-           c(255, 106,  71)  # label for low RIN is red


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
#' Note: This function attempts to find the marker bands by their special color. But because this color is overlaid on the gel image, it is impossible to read the fluorescence intensity inside the markers, so it is reported as NA.
#'
#' @param gel.image.file The filename of a TapeStation gel image with blue highlight, in PNG format.
#'
#' @return A data frame with one row for each vertical pixel of each gel lane.
#' 
#' @seealso \code{\link{read.tapestation}}, \code{\link{read.tapestation.xml}}
#'
#' @export
#' @importFrom png readPNG
read.tapestation.gel.image <- function(gel.image.file) {
	gel.image.rgb <- readPNG(gel.image.file) # this is in the form (y, x, channel); [1,1,] is the upper left corner
	
	# find lower marker
	pixel.is.lower.marker <- find.matching.pixels(gel.image.rgb, RGB.LOWER.MARKER)
	x.is.lower.marker <- apply(pixel.is.lower.marker, 2, any)
	
	# use lower marker bands to identify a single representative x-value for each gel lane
	x.gel <- which(diff(c(FALSE, x.is.lower.marker)) == 1) # guess the start x-positions of the lanes based on where there are gaps between lower markers; add this FALSE so that column 1 will test positive if necessary, and this also offsets all the indices correctly
	lane.spacings <- diff(x.gel)
	if (diff(range(lane.spacings)) > 1) { # the guessed lanes are not evenly spaced, even tolerating 1 pixel of antialiasing error
		lane.spacing.freq <- table(lane.spacings)
		estimated.lane.width <- as.integer(names(lane.spacing.freq)[which.max(lane.spacing.freq)]) # guess that the most common width is the correct one (assuming there aren't lots of problems in this batch!)
		n.lanes <- floor((dim(gel.image.rgb)[2] - x.gel[1]) / estimated.lane.width) # assume the first edge is definitely called correctly (should only be whitespace to its left) and any whitespace to the right is narrower than a lane
		x.gel <- x.gel[1] + estimated.lane.width * 1:n.lanes - round(estimated.lane.width / 2) # aim for the centers of the lanes because we might have slight error
	}
	gel.image.rgb.reduced <- gel.image.rgb[,x.gel,]
			
	# find gel boundaries
	position.is.highlight <- find.matching.pixels(gel.image.rgb.reduced, RGB.HIGHLIGHT)
	lane.with.borders <- which(position.is.highlight[1,]) # assuming there will be a highlight in the top pixel row
	stopifnot(length(lane.with.borders) == 1) # need one highlighted lane to find gel borders
	position.is.label <- find.matching.pixels(gel.image.rgb.reduced, RGB.GOOD) |
		find.matching.pixels(gel.image.rgb.reduced, RGB.MEDIUM) |
		find.matching.pixels(gel.image.rgb.reduced, RGB.BAD)
	border.transition <- diff((position.is.highlight | position.is.label)[,lane.with.borders])
	y.gel.start <- which(border.transition == -1)[1] - 1
	y.gel.end <- which(border.transition == 1)[1] - 1

	# now finally extract the intensities!
	result.rgb <- gel.image.rgb.reduced[y.gel.end:y.gel.start,,] # only the actual data values; reverse order so it goes bottom to top like peak calls and distance
	fluorescence.matrix <- 1 - result.rgb[,,1] # only get red fluorescence because all channels are equal in the places we care about; subtract from 1 because it's a negative (red decreases in the protein gels too even though they're blue instead of black)
	fluorescence.matrix[result.rgb[,,1] != result.rgb[,,2]] <- NA # set non-data pixels (obscured by marker band color) to NA; assume red channel always equals green channel, but not necessary blue because protein gels use blue
	data.frame(
		sample.index =  rep(1:length(x.gel), each = nrow(fluorescence.matrix)),
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
#' @param xml.file The filename of an XML file exported from the TapeStation software. The filename is expected to end in \code{.xml} and the name before that extension is used as the name of the batch.
#'
#' @return A list of some of the components of an \code{electrophoresis} object
#' 
#' @seealso \code{\link{read.tapestation}}, \code{\link{read.tapestation.gel.image}}
#'
#' @export
#' @importFrom XML xmlRoot xmlParse xmlValue xmlApply xmlChildren xmlToDataFrame
read.tapestation.xml <- function(xml.file) {
 	batch <- sub("\\.xml$", "", basename(xml.file))
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
	if (is.null(gel.image.file)) gel.image.file <- sub("\\.xml$", ".png", xml.file)
	
	parsed.data <- read.tapestation.xml(xml.file)
	stopifnot(length(unique(parsed.data$samples$batch)) == 1)
	batch <- parsed.data$samples$batch[1]
	gel.data <- read.tapestation.gel.image(gel.image.file)
	stopifnot(length(unique(gel.data$sample.index)) == nrow(parsed.data$samples))
	result <- structure(list(
		data = gel.data,
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

