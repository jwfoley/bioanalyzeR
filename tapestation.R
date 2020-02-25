library(XML)
library(png)

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


read.tapestation.gel.image <- function(gel.image.files) {

	combined.results <- do.call(rbind, lapply(gel.image.files, function(gel.image.file) {

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
		result.rgb <- gel.image.rgb.reduced[y.gel.start:y.gel.end,,] # only the actual data values
		fluorescence.matrix <- 1 - result.rgb[,,1] # only get red fluorescence because all channels are equal in the places we care about; subtract from 1 because it's a negative (red decreases in the protein gels too even though they're blue instead of black)
		fluorescence.matrix[result.rgb[,,1] != result.rgb[,,2]] <- NA # set non-data pixels (obscured by marker band color) to NA; assume red channel always equals green channel, but not necessary blue because protein gels use blue
		results <- data.frame(
			batch =         sub("\\.png$", "", basename(gel.image.file)),
			gel.lane =      rep(1:length(x.gel), each = nrow(fluorescence.matrix)),
			distance =      1:nrow(fluorescence.matrix) / nrow(fluorescence.matrix),
			fluorescence =  as.vector(fluorescence.matrix)
		)
	}))
	
	rownames(combined.results) <- NULL
	combined.results$batch <- factor(combined.results$batch, levels = unique(combined.results$batch)) # make batches into a factor that keeps them in the observed order
	
	combined.results
}


read.tapestation.xml <- function(...) {
	result.list <- do.call(c, lapply(list(...), function(xml.file) {
		xml.root <- xmlRoot(xmlParse(xml.file))
		xmlApply(xml.root[["Samples"]], function(sample.xml) {
			well.number <- xmlValue(sample.xml[["WellNumber"]])
			name <- trimws(xmlValue(sample.xml[["Comment"]]))
			sample.observations <- trimws(xmlValue(sample.xml[["Observations"]]))
			if (sample.observations == "Marker(s) not detected") {
				warning(paste(sample.observations, "for well", well.number, name))
				return(NULL)
			}
			well.row <- gsub("[^A-HL]", "", well.number) # assumes all row names are in A through H, or L for ladder
			well.col <- as.integer(gsub("[^[:digit:]]", "", well.number)) # assumes all column names are numbers!
			
			reagent.id <- xmlValue(sample.xml[["ScreenTapeID"]])
			
			suppressWarnings( # will throw warnings if missing values are coerced to NA but we can live with that
				peaks <- if (length(xmlChildren(sample.xml[["Peaks"]])) == 0) NULL else do.call(rbind, c(xmlApply(sample.xml[["Peaks"]], function(peak.xml) data.frame(
					peak.observations =  trimws(xmlValue(peak.xml[["Observations"]])),
					peak.comment =       trimws(xmlValue(peak.xml[["Comment"]])),
					length =             as.integer(xmlValue(peak.xml[["Size"]])),
					distance =           as.numeric(xmlValue(peak.xml[["RunDistance"]])) / 100,
					lower.bound =        as.numeric(xmlValue(peak.xml[["FromPercent"]])) / 100,
					upper.bound =        as.numeric(xmlValue(peak.xml[["ToPercent"]])) / 100, 
					area =               as.numeric(xmlValue(peak.xml[["Area"]])),
					molarity =           as.numeric(xmlValue(peak.xml[["Molarity"]]))
				)), make.row.names = F))
			)
			
			suppressWarnings( # will throw warnings if missing values are coerced to NA but we can live with that
				regions <- if (length(xmlChildren(sample.xml[["Regions"]])) == 0) NULL else do.call(rbind, c(xmlApply(sample.xml[["Regions"]], function(region.xml) data.frame(
					peak.comment =         trimws(xmlValue(region.xml[["Comment"]])),
					lower.bound =          as.integer(xmlValue(region.xml[["From"]])),
					upper.bound =          as.integer(xmlValue(region.xml[["To"]])),
					average.length =       as.integer(xmlValue(region.xml[["AverageSize"]])),
					area =                 as.numeric(xmlValue(region.xml[["Area"]])),
					concentration =        as.numeric(xmlValue(region.xml[["Concentration"]])),
					molarity =             as.numeric(xmlValue(region.xml[["Molarity"]])),
					proportion.of.total =  as.numeric(xmlValue(region.xml[["PercentOfTotal"]])) / 100				
				)), make.row.names = F))
			)
			
			list(
				peaks = if (is.null(peaks)) NULL else cbind(batch = sub("\\.xml$", "", basename(xml.file)), well.number, well.row, well.col, name, reagent.id, sample.observations, peaks),
				regions = if (is.null(regions)) NULL else cbind(batch = sub("\\.xml$", "", basename(xml.file)), well.number, well.row, well.col, name, reagent.id, sample.observations, regions)
			)
		})
	}))
	
	peaks.list <- lapply(result.list, function(x) x$peaks)
	regions.list <- lapply(result.list, function(x) x$regions)
	list(
		peaks = if (all(unlist(lapply(peaks.list, is.null)))) NULL else do.call(rbind, c(peaks.list, make.row.names = F)),
		regions = if (all(unlist(lapply(regions.list, is.null)))) NULL else do.call(rbind, c(regions.list, make.row.names = F))
	)
}


read.tapestation <- function(xml.file, gel.image.file = NULL, fit = "spline") {
	stopifnot(fit %in% c("interpolate", "spline", "regression"))
	if (is.null(gel.image.file)) gel.image.file <- sub("\\.xml$", ".png", xml.file)
	
	parsed.data <- read.tapestation.xml(xml.file)
	peaks <- parsed.data$peaks
	result <- read.tapestation.gel.image(gel.image.file)
	stopifnot(length(unique(result$gel.lane)) == length(unique(peaks$well.number)))
	
	# get sample names and observations (might not be unique but we assume well numbers are)
	sample.table <- unique(peaks[,c("batch", "reagent.id", "name", "sample.observations", "well.number", "well.row", "well.col")])
	stopifnot(nrow(sample.table) == length(unique(result$gel.lane)))
	
	result <- cbind(sample.table[result$gel.lane,], subset(result, select = -batch))
	
	# calculate relative distances
	lower.marker.peaks <- subset(peaks, peak.observations %in% c("Lower Marker", "edited Lower Marker"))
	marker.distances <- data.frame(lower = sapply(unique(peaks$well.number), function(this.well.number) {
		distance <- subset(lower.marker.peaks, well.number == this.well.number)$distance
		if (length(distance) == 0) {
			return(NA)
		} else if (length(distance) == 1) {
			return(distance)
		} else {
			stop(paste("multiple lower marker peaks for well", this.well.number))
		}
	}), row.names = unique(peaks$well.number))
	upper.marker.peaks <- subset(peaks, peak.observations %in% c("Upper Marker", "edited Upper Marker"))
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
	result$relative.distance <- (result$distance - marker.distances$upper[result$well.number]) / marker.distances$range[result$well.number]
	peaks$relative.distance <- (peaks$distance - marker.distances$upper[peaks$well.number]) / marker.distances$range[peaks$well.number]
	
	# allocate new variables, which might or might not get populated
	result$length <- NA
	result <- cbind(result, do.call(rbind, lapply(unique(result$well.number), function(this.well) {
		result.this.well <- subset(result, well.number == this.well)
		data.frame(
			delta.distance = c(NA, diff(result.this.well$distance)),
			delta.fluorescence = c(NA, diff(result.this.well$fluorescence))
		)
	})))
	result$delta.area <- (2 * result$fluorescence - result$delta.fluorescence) / 2 * result$delta.distance
	result$delta.mass <- NA
	result$delta.molarity <- NA
	
	# determine ladder scheme
	ladder.wells <- as.character(unique(subset(peaks, sample.observations == "Ladder")$well.number))
	ladder.rows <- list()
	
	# scheme: no ladder	
	if (length(ladder.wells) == 0) {
		warning("warning: no ladder specified so lengths and molarities are not calculated")
	
	# scheme: only one ladder for the whole run
	} else if (length(ladder.wells) == 1) {
		ladder.rows[[ladder.wells]] <- 1:nrow(result) # all rows of the frame are associated with this ladder
	
	# scheme: one electronic ladder for the whole run but it's displayed more than once
	} else if (
		length(unique(subset(peaks, well.number %in% ladder.wells)$reagent.id)) == 1 &&
		unique(subset(peaks, well.number %in% ladder.wells)$name) == "Electronic Ladder"
	) {
		ladder.rows[[ladder.wells[1]]] <- 1:nrow(result) # all rows of the frame are associated with this ladder
	
	# scheme: one ladder per reagent.id (ScreenTape)
	} else if (all.equal(unique(peaks$reagent.id), unique(subset(peaks, well.number %in% ladder.wells)$reagent.id))) { # 
			ladder.rows <- lapply(ladder.wells, function(well) which(result$reagent.id == as.character(unique(subset(peaks, well.number == well)$reagent.id))))
			names(ladder.rows) <- ladder.wells
	
	# scheme: something unexpected
	}	else {
			warning("warning: unknown ladder scheme so lengths and molarities are not calculated")
	}
	
	for (ladder.well in names(ladder.rows)) { # use names(ladder.rows) because it only exists if we're in a recognized scheme
		
		# fit standard curve for molecule length
		# do this in relative-distance space so it's effectively recalibrated for each sample's markers

		peaks.ladder <- subset(peaks, well.number == ladder.well)
		which.rows <- ladder.rows[[ladder.well]]
		if (fit == "interpolate") {
			warning("linear interpolation gives ugly results for molarity estimation")
			standard.curve.function <- approxfun(peaks.ladder$relative.distance, peaks.ladder$length)
		} else if (fit == "spline") {
			standard.curve.function <- splinefun(peaks.ladder$relative.distance, peaks.ladder$length, method = "natural")
		} else if (fit == "regression") {
			mobility.model <- lm(relative.distance ~ log(length), peaks.ladder)
			standard.curve.function <- function(distance) exp((distance - mobility.model$coefficients[1]) / mobility.model$coefficients[2])
		}
		result$length[which.rows] <- standard.curve.function(result$relative.distance[which.rows])
		
		# convert to molarity
		peaks.ladder$mass <- peaks.ladder$length * peaks.ladder$molarity
		peaks.ladder$estimated.area <- c(NA, sapply(2:(nrow(peaks.ladder) - 1), function(i) {
			sum(subset(result, well.number == ladder.well & distance >= peaks.ladder$lower.bound[i] & distance <= peaks.ladder$upper.bound[i])$delta.area)
		}), NA)
		fluorescence.model <- lm(mass ~ estimated.area - 1, data = peaks.ladder)
		result$delta.mass[which.rows] <- fluorescence.model$coefficients[1] * result$delta.area[which.rows]
		result$delta.molarity[which.rows] <- result$delta.mass[which.rows] / result$length[which.rows]
	}
	
	rownames(result) <- NULL # clean up row names again
	
	structure(list(data = result, peaks = peaks, regions = parsed.data$regions), class = "tapestation")
}

