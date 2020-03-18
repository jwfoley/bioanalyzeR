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
		gel.lane =      rep(1:length(x.gel), each = nrow(fluorescence.matrix)),
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
		is.ladder <- sample.observations == "Ladder"
		if (sample.observations == "Marker(s) not detected") {
			warning(paste(sample.observations, "for well", well.number, sample.name))
			return(NULL)
		}
		
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
			sample.info = data.frame(batch, well.number, sample.name, reagent.id, is.ladder, sample.observations, stringsAsFactors = F),
			peaks = if (is.null(peaks)) NULL else data.frame(batch, well.number, sample.name, reagent.id, is.ladder, sample.observations, peaks, stringsAsFactors = F),
			regions = if (is.null(regions)) NULL else data.frame(batch, well.number, sample.name, reagent.id, is.ladder, sample.observations, regions, stringsAsFactors = F)
		)
	})
	
	peaks.list <- lapply(result.list, function(x) x$peaks)
	regions.list <- lapply(result.list, function(x) x$regions)
	result <- list(
		samples = do.call(rbind, c(lapply(result.list, function(x) x$sample.info), make.row.names = F)),
		peaks = if (all(unlist(lapply(peaks.list, is.null)))) NULL else do.call(rbind, c(peaks.list, make.row.names = F)),
		regions = if (all(unlist(lapply(regions.list, is.null)))) NULL else do.call(rbind, c(regions.list, make.row.names = F)),
		assay.info = assay.info
	)
	
	# convert sample metadata into factors, ensuring all frames have the same levels and the levels are in the observed order
	for (field in c("batch", "well.number", "sample.name", "reagent.id", "sample.observations")) {
		result$samples[,field] <- factor(result$samples[,field], levels = unique(result$samples[,field]))
		result$peaks[,field] <- factor(result$peaks[,field], levels = levels(result$samples[,field]))
		result$regions[,field] <- factor(result$regions[,field], levels = levels(result$samples[,field]))
	}
	# convert other text into factors without those restrictions
	for (field in c("peak.observations", "peak.comment")) result$peaks[,field] <- factor(result$peaks[,field])
	for (field in c("region.comment")) result$regions[,field] <- factor(result$regions[,field])
	
	result
}


#' Read TapeStation data
#'
#' This function reads TapeStation run data from an exported XML file and a PNG gel image, then fills out the results with estimates of molecule length, concentration, and molarity.
#'
#' Spline fitting seems to perform reasonably well on all data. Agilent appears to use linear interpolation with DNA data and log-linear regression on RNA data, so you could choose those options if you want to reproduce the results of the software more precisely. However, linear interpolation creates sudden spikes in the derivative that make the concentration and molarity estimates unstable; spline fitting is basically a smoother version of that. Log-linear regression is the standard theoretical approach but does not actually fit the data very well; more sophisticated parametric models may be added in the future.
#'
#' @param xml.file The filename of an XML file exported from the TapeStation software. The filename is expected to end in \code{.xml} and the name before that extension is used as the name of the batch.
#' @param gel.image.file The filename of a TapeStation gel image with blue highlight, in PNG format. If \code{NULL}, the gel image file is expected to have the same name as the XML file with a different extension, e.g. \code{experiment1.xml} and \code{experiment1.png}, so if you name your files in that pattern you don't need to fill out this argument.
#' @param fit The method used to fit the mobility model of molecule length vs. migration distance, one of \code{"interpolation"} (linear interpolation via \code{\link{approxfun}}), \code{"spline"} (splines via \code{\link{splinefun}}), or \code{"regression"} (log-linear regression via \code{\link{lm}} with the model \code{relative.distance ~ log(length)}).
#'
#' @return An \code{electrophoresis} object containing the data from this TapeStation run.
#' 
#' @seealso \code{\link{read.bioanalyzer}}, \code{\link{read.tapestation.gel.image}}, \code{\link{read.tapestation.xml.file}}
#'
#' @export
read.tapestation <- function(xml.file, gel.image.file = NULL, fit = "spline") {
	stopifnot(fit %in% c("interpolation", "spline", "regression"))
	if (is.null(gel.image.file)) gel.image.file <- sub("\\.xml$", ".png", xml.file)
	
	parsed.data <- read.tapestation.xml(xml.file)
	stopifnot(length(unique(parsed.data$samples$batch)) == 1)
	batch <- parsed.data$samples$batch[1]
	gel.data <- read.tapestation.gel.image(gel.image.file)
	stopifnot(length(unique(gel.data$gel.lane)) == nrow(parsed.data$samples))
	result <- structure(list(
		data = cbind(parsed.data$samples[gel.data$gel.lane,], gel.data[,colnames(gel.data) != "gel.lane"]),
		assay.info = list(parsed.data$assay.info),
		samples = parsed.data$samples,
		wells.by.ladder = NULL,
		peaks = parsed.data$peaks,
		regions = parsed.data$regions,
		mobility.functions = NULL,
		mass.coefficients = NULL
	), class = "electrophoresis")
	names(result$assay.info) <- batch
	
	# calculate relative distances
	lower.marker.peaks <- subset(result$peaks, peak.observations %in% c("Lower Marker", "edited Lower Marker"))
	marker.distances <- data.frame(lower = sapply(unique(result$peaks$well.number), function(this.well.number) {
		distance <- subset(lower.marker.peaks, well.number == this.well.number)$distance
		if (length(distance) == 0) {
			return(NA)
		} else if (length(distance) == 1) {
			return(distance)
		} else {
			stop(paste("multiple lower marker peaks for well", this.well.number))
		}
	}), row.names = unique(result$peaks$well.number))
	upper.marker.peaks <- subset(result$peaks, peak.observations %in% c("Upper Marker", "edited Upper Marker"))
	if (nrow(upper.marker.peaks) == 0) { # kit lacks upper marker
		marker.distances$upper <- 0 # effectively normalizes only to lower marker
	} else {
		marker.distances$upper <- sapply(rownames(marker.distances), function(this.well.number) {
			distance <- subset(upper.marker.peaks, well.number == this.well.number)$distance
			if (length(distance) == 0) {
				return(NA)
			} else if (length(distance) == 1) {
				return(distance)
			} else {
				stop(paste("multiple lower marker peaks for well", this.well.number))
			}
		})
	}
	marker.distances$range <- marker.distances$lower - marker.distances$upper
	result$data$relative.distance <- (result$data$distance - marker.distances$upper[result$data$well.number]) / marker.distances$range[result$data$well.number]
	result$peaks$relative.distance <- (result$peaks$distance - marker.distances$upper[result$peaks$well.number]) / marker.distances$range[result$peaks$well.number]
	result$peaks$lower.relative.distance <- (result$peaks$lower.distance - marker.distances$upper[result$peaks$well.number]) / marker.distances$range[result$peaks$well.number]
	result$peaks$upper.relative.distance <- (result$peaks$upper.distance - marker.distances$upper[result$peaks$well.number]) / marker.distances$range[result$peaks$well.number]
	
	# determine ladder scheme
	ladder.wells <- as.character(subset(result$samples, is.ladder)$well.number)
	result$wells.by.ladder <- list(list())
	names(result$wells.by.ladder) <- batch
	rows.by.ladder <- list()
	# scheme: no ladder	
	if (length(ladder.wells) == 0) {
		warning("warning: no ladder specified so lengths and molarities are not calculated")
	# scheme: only one ladder for the whole run
	} else if (length(ladder.wells) == 1) {
		result$wells.by.ladder[[1]][[ladder.wells]] <- unique(result$data$well.number)
		rows.by.ladder[[ladder.wells]] <- 1:nrow(result$data) # all rows of the frame are associated with this ladder
	# scheme: one electronic ladder for the whole run but it's displayed more than once
	} else if (
		length(unique(subset(result$peaks, well.number %in% ladder.wells)$reagent.id)) == 1 &&
		unique(subset(result$peaks, well.number %in% ladder.wells)$sample.name) == "Electronic Ladder"
	) {
		result$wells.by.ladder[[1]][[ladder.wells[1]]] <- unique(result$data$well.number)
		rows.by.ladder[[ladder.wells[1]]] <- 1:nrow(result$data) # all rows of the frame are associated with this ladder
	# scheme: one ladder per reagent.id (ScreenTape)
	} else if (all.equal(unique(result$peaks$reagent.id), unique(subset(result$peaks, well.number %in% ladder.wells)$reagent.id))) { # 
		rows.by.ladder <- lapply(ladder.wells, function(well) which(result$data$reagent.id == unique(subset(result$peaks, well.number == well)$reagent.id)))
		names(rows.by.ladder) <- ladder.wells
		result$wells.by.ladder[[1]] <- lapply(rows.by.ladder, function(x) unique(result$data[x,]$well.number))
	# scheme: something unexpected
	}	else {
		warning("warning: unknown ladder scheme so lengths and molarities are not calculated")
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
	result$mass.coefficients <- list(rep(NA, nrow(result$samples)))
	names(result$mass.coefficients) <- batch
	names(result$mass.coefficients[[1]]) <- result$samples$well.number
	data.calibration <- cbind(result$data, do.call(rbind, lapply(result$samples$well.number, function(this.well) {
		result.this.well <- subset(result$data, well.number == this.well)
		data.frame(
			delta.fluorescence = c(NA, diff(result.this.well$fluorescence)),
			delta.distance = c(NA, -diff(result.this.well$distance))
		)
	})))
	# estimate area under each measurement with the trapezoidal rule; to simplify math, each point's sum is for the trapezoid to the left of it
	data.calibration$area <- (2 * data.calibration$fluorescence - data.calibration$delta.fluorescence) * data.calibration$delta.distance
	for (ladder.well in names(rows.by.ladder)) { # use names(rows.by.ladder) because it only exists if we're in a recognized scheme
		peaks.ladder <- subset(result$peaks, well.number == ladder.well)
		which.rows <- rows.by.ladder[[ladder.well]]
		which.peaks <- which(result$peaks$well.number %in% result$wells.by.ladder[[1]][[ladder.well]])
		which.regions <- which(result$regions$well.number %in% result$wells.by.ladder[[1]][[ladder.well]])
		
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
		peaks.ladder$area <- sapply(1:nrow(peaks.ladder), function(i) sum(subset(data.calibration, well.number == ladder.well & distance >= peaks.ladder$lower.distance[i] & distance <= peaks.ladder$upper.distance[i])$area))
		mass.coefficient <- lm(concentration ~ area - 1, data = peaks.ladder)$coefficients[1]
		result$mass.coefficients[[1]][result$wells.by.ladder[[1]][[ladder.well]]] <- mass.coefficient
		result$data$concentration[which.rows] <- mass.coefficient * data.calibration$area[which.rows]
		result$data$molarity[which.rows] <- result$data$concentration[which.rows] / molecular.weight(result$data$length[which.rows], parsed.data$assay.info$assay.type) * 1E6 # we're converting ng/uL to nmol/L or pg/uL to pmol/L so we need to scale by 1E6
	}
	
	# convert inferred relative distances of regions back to raw distances
	if (! is.null(result$regions)) {
		result$regions$lower.distance <- result$regions$lower.relative.distance * marker.distances$range[result$regions$well.number] + marker.distances$upper[result$regions$well.number]
		result$regions$upper.distance <- result$regions$upper.relative.distance * marker.distances$range[result$regions$well.number] + marker.distances$upper[result$regions$well.number]
	}
	
	rownames(result$data) <- NULL # clean up row names again
		
	result
}
