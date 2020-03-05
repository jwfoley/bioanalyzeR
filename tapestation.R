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


read.tapestation.xml <- function(xml.file) {
 	batch <- sub("\\.xml$", "", basename(xml.file))
	result.list <- xmlApply(xmlRoot(xmlParse(xml.file))[["Samples"]], function(sample.xml) {
		well.number <- xmlValue(sample.xml[["WellNumber"]])
		sample.name <- trimws(xmlValue(sample.xml[["Comment"]]))
		if (sample.name == "") sample.name <- well.number
		sample.observations <- trimws(xmlValue(sample.xml[["Observations"]]))
		if (sample.observations == "Marker(s) not detected") {
			warning(paste(sample.observations, "for well", well.number, sample.name))
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
				lower.distance =     as.numeric(xmlValue(peak.xml[["FromPercent"]])) / 100,
				upper.distance =     as.numeric(xmlValue(peak.xml[["ToPercent"]])) / 100, 
				area =               as.numeric(xmlValue(peak.xml[["Area"]])),
				molarity =           as.numeric(xmlValue(peak.xml[["Molarity"]])),
				stringsAsFactors =   F
			)), make.row.names = F))
		)
		
		suppressWarnings( # will throw warnings if missing values are coerced to NA but we can live with that
			regions <- if (length(xmlChildren(sample.xml[["Regions"]])) == 0) NULL else do.call(rbind, c(xmlApply(sample.xml[["Regions"]], function(region.xml) data.frame(
				region.comment =       trimws(xmlValue(region.xml[["Comment"]])),
				lower.length =         as.integer(xmlValue(region.xml[["From"]])),
				upper.length =         as.integer(xmlValue(region.xml[["To"]])),
				average.length =       as.integer(xmlValue(region.xml[["AverageSize"]])),
				area =                 as.numeric(xmlValue(region.xml[["Area"]])),
				concentration =        as.numeric(xmlValue(region.xml[["Concentration"]])),
				molarity =             as.numeric(xmlValue(region.xml[["Molarity"]])),
				proportion.of.total =  as.numeric(xmlValue(region.xml[["PercentOfTotal"]])) / 100,
				stringsAsFactors =     F
			)), make.row.names = F))
		)
		
		list(
			sample.info = data.frame(batch, well.number, well.row, well.col, sample.name, reagent.id, sample.observations, stringsAsFactors = F),
			peaks = if (is.null(peaks)) NULL else data.frame(batch, well.number, well.row, well.col, sample.name, reagent.id, sample.observations, peaks, stringsAsFactors = F),
			regions = if (is.null(regions)) NULL else data.frame(batch, well.number, well.row, well.col, sample.name, reagent.id, sample.observations, regions, stringsAsFactors = F)
		)
	})
	
	peaks.list <- lapply(result.list, function(x) x$peaks)
	regions.list <- lapply(result.list, function(x) x$regions)
	result <- list(
		samples = do.call(rbind, c(lapply(result.list, function(x) x$sample.info), make.row.names = F)),
		peaks = if (all(unlist(lapply(peaks.list, is.null)))) NULL else do.call(rbind, c(peaks.list, make.row.names = F)),
		regions = if (all(unlist(lapply(regions.list, is.null)))) NULL else do.call(rbind, c(regions.list, make.row.names = F))
	)
	
	# convert sample metadata into factors, ensuring all frames have the same levels and the levels are in the observed order
	for (field in colnames(result$samples)) {
		result$samples[,field] <- factor(result$samples[,field], levels = unique(result$samples[,field]))
		result$peaks[,field] <- factor(result$peaks[,field], levels = levels(result$samples[,field]))
		result$regions[,field] <- factor(result$regions[,field], levels = levels(result$samples[,field]))
	}
	# convert other text into factors without those restrictions
	for (field in c("peak.observations", "peak.comment")) result$peaks[,field] <- factor(result$peaks[,field])
	for (field in c("region.comment")) result$regions[,field] <- factor(result$regions[,field])
	
	result
}


read.tapestation <- function(xml.file, gel.image.file = NULL, fit = "spline") {
	stopifnot(fit %in% c("interpolation", "spline", "regression"))
	if (is.null(gel.image.file)) gel.image.file <- sub("\\.xml$", ".png", xml.file)
	
	parsed.data <- read.tapestation.xml(xml.file)
	samples <- parsed.data$samples
	stopifnot(length(unique(samples$batch)) == 1)
	batch = samples$batch[1]
	peaks <- parsed.data$peaks
	regions <- parsed.data$regions
	result <- read.tapestation.gel.image(gel.image.file)
	stopifnot(length(unique(result$gel.lane)) == nrow(samples))
	
	result <- cbind(samples[result$gel.lane,], result)
	
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
	peaks$lower.relative.distance <- (peaks$lower.distance - marker.distances$upper[peaks$well.number]) / marker.distances$range[peaks$well.number]
	peaks$upper.relative.distance <- (peaks$upper.distance - marker.distances$upper[peaks$well.number]) / marker.distances$range[peaks$well.number]
	
	# allocate new variables, which might or might not get populated
	result$length <- NA
	result <- cbind(result, do.call(rbind, lapply(unique(result$well.number), function(this.well) {
		result.this.well <- subset(result, well.number == this.well)
		data.frame(
			delta.distance = c(NA, diff(result.this.well$distance)),
			delta.fluorescence = c(NA, diff(result.this.well$fluorescence))
		)
	})))
	result$area <- (2 * result$fluorescence - result$delta.fluorescence) / 2 * result$delta.distance
	result$mass <- NA
	result$molarity <- NA
	
	# determine ladder scheme
	ladder.wells <- as.character(subset(samples, sample.observations == "Ladder")$well.number)
	wells.by.ladder <- list()
	rows.by.ladder <- list()
	
	# scheme: no ladder	
	if (length(ladder.wells) == 0) {
		warning("warning: no ladder specified so lengths and molarities are not calculated")
	
	# scheme: only one ladder for the whole run
	} else if (length(ladder.wells) == 1) {
		wells.by.ladder[[ladder.wells]] <- unique(result$well.number)
		rows.by.ladder[[ladder.wells]] <- 1:nrow(result) # all rows of the frame are associated with this ladder
	
	# scheme: one electronic ladder for the whole run but it's displayed more than once
	} else if (
		length(unique(subset(peaks, well.number %in% ladder.wells)$reagent.id)) == 1 &&
		unique(subset(peaks, well.number %in% ladder.wells)$sample.name) == "Electronic Ladder"
	) {
		wells.by.ladder[[ladder.wells[1]]] <- unique(result$well.number)
		rows.by.ladder[[ladder.wells[1]]] <- 1:nrow(result) # all rows of the frame are associated with this ladder
	
	# scheme: one ladder per reagent.id (ScreenTape)
	} else if (all.equal(unique(peaks$reagent.id), unique(subset(peaks, well.number %in% ladder.wells)$reagent.id))) { # 
		rows.by.ladder <- lapply(ladder.wells, function(well) which(result$reagent.id == unique(subset(peaks, well.number == well)$reagent.id)))
		names(rows.by.ladder) <- ladder.wells
		wells.by.ladder <- lapply(rows.by.ladder, function(x) unique(result[x,]$well.number))
	
	# scheme: something unexpected
	}	else {
		warning("warning: unknown ladder scheme so lengths and molarities are not calculated")
	}
	
	peaks$lower.length <- NA
	peaks$upper.length <- NA
	if (! is.null(regions)) {
		regions$lower.relative.distance <- NA
		regions$upper.relative.distance <- NA
	}
	mobility.functions <- list()
	mobility.inverses <- list()
	mass.coefficients <- list()
	for (ladder.well in names(rows.by.ladder)) { # use names(rows.by.ladder) because it only exists if we're in a recognized scheme
		peaks.ladder <- subset(peaks, well.number == ladder.well)
		which.rows <- rows.by.ladder[[ladder.well]]
		which.peaks <- which(peaks$well.number %in% wells.by.ladder[[ladder.well]])
		which.regions <- which(regions$well.number %in% wells.by.ladder[[ladder.well]])
		
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
		mobility.functions[[ladder.well]] <- standard.curve.function
		mobility.inverses[[ladder.well]] <- standard.curve.inverse
		
		# apply model to raw data
		result$length[which.rows] <- standard.curve.function(result$relative.distance[which.rows])
		
		# apply model to peaks
		peaks$lower.length[which.peaks] <- standard.curve.function(peaks$upper.relative.distance[which.peaks])
		peaks$upper.length[which.peaks] <- standard.curve.function(peaks$lower.relative.distance[which.peaks])
		
		# apply inverse model to regions
		if (! is.null(regions)) {
			regions$lower.relative.distance[which.regions] <- standard.curve.inverse(regions$upper.length[which.regions])
			regions$upper.relative.distance[which.regions] <- standard.curve.inverse(regions$lower.length[which.regions])
		}
		
		# convert to molarity
		# first solve for a coefficient that relates known masses (known molarities scaled by length) to area under the electropherogram peaks, then apply that coefficient to each individual measurement
		# this is done by fitting a one-parameter model on the non-marker peaks of the ladder, because the marker peaks in all samples are unreadable (blocked by green and purple bands)
		# preferably it would be scaled to each sample's own markers, but that's impossible because the markers are unreadable! and we can't even use their Agilent-reported molarity estimates or areas to scale relative to the ladder because those are always normalized to the upper marker! (which is often the one that's contaminated by sample anyway)
		# so all we can do is normalize to the ladder's non-marker peaks, therefore each sample will be randomly off by some constant scaling factor, but at least molarity comparisons within a sample ought to be accurate
		peaks.ladder$mass <- peaks.ladder$length * peaks.ladder$molarity
		peaks.ladder$estimated.area <- sapply(1:nrow(peaks.ladder), function(i) sum(subset(result, well.number == ladder.well & distance >= peaks.ladder$lower.distance[i] & distance <= peaks.ladder$upper.distance[i])$area))
		fluorescence.model <- lm(mass ~ estimated.area - 1, data = peaks.ladder)
		mass.coefficient <- fluorescence.model$coefficients[1]
		mass.coefficients[[ladder.well]] <- mass.coefficient
		result$mass[which.rows] <- mass.coefficient * result$area[which.rows]
		result$molarity[which.rows] <- result$mass[which.rows] / result$length[which.rows]
	}
	
	# convert inferred relative distances of regions back to raw distances
	if (! is.null(regions)) {
		regions$lower.distance <- regions$lower.relative.distance * marker.distances$range[regions$well.number] + marker.distances$upper[regions$well.number]
		regions$upper.distance <- regions$upper.relative.distance * marker.distances$range[regions$well.number] + marker.distances$upper[regions$well.number]
	}
	
	# annotate which peak each data point is in, if any
	# WARNING: if peaks overlap, this will overwrite and each point will only be mapped to the last-occuring one!
	# WARNING: if multiple tables are combined later, these will all be wrong!
	result$peak <- NA
	for (i in 1:nrow(peaks)) result$peak[result$well.number == peaks$well.number[i] & result$distance >= peaks$lower.distance[i] & result$distance <= peaks$upper.distance[i]] <- i
	
	rownames(result) <- NULL # clean up row names again
	
	# wrap lists by ladder in larger lists by batch so they'll survive being combined with other batches
	wells.by.ladder <- list(wells.by.ladder)
	names(wells.by.ladder) <- batch
	mobility.functions <- list(mobility.functions)
	names(mobility.functions) <- batch
	mass.coefficients <- list(mass.coefficients)
	names(mass.coefficients) <- batch
	
	structure(list(
		data = result,
		samples = samples,
		wells.by.ladder = wells.by.ladder,
		peaks = peaks,
		regions = regions,
		mobility.functions = mobility.functions,
		mass.coefficients = mass.coefficients
	), class = "electrophoresis")
}

