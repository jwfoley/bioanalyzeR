library(XML)
library(png)

rgb.upper.marker <- c(128, 0, 128) / 255 # upper marker is purple
rgb.lower.marker <- c(0, 128, 0) / 255 # lower marker is green


read.tapestation.gel.image <- function(gel.image.files) {

	combined.results <- do.call(rbind, lapply(gel.image.files, function(gel.image.file) {

		gel.image.rgb <- readPNG(gel.image.file) # this is in the form (y, x, channel); [1,1,] is the upper left corner

		# find marker bands
		which.pixel.upper.marker <- gel.image.rgb[,,1] == rgb.upper.marker[1] & gel.image.rgb[,,2] == rgb.upper.marker[2] & gel.image.rgb[,,3] == rgb.upper.marker[3]
		which.pixel.lower.marker <- gel.image.rgb[,,1] == rgb.lower.marker[1] & gel.image.rgb[,,2] == rgb.lower.marker[2] & gel.image.rgb[,,3] == rgb.lower.marker[3]

		# use marker bands to identify gel lanes
		which.x.upper.marker <- apply(which.pixel.upper.marker, 2, any)
		which.x.lower.marker <- apply(which.pixel.lower.marker, 2, any)
		stopifnot(which.x.upper.marker == which.x.lower.marker)

		# isolate a single representative x-value (the first one from the left) for each gel lane
		x.gel <- which(diff(c(FALSE, which.x.upper.marker)) == 1) # add this FALSE so that column 1 will test positive if necessary, and this also offsets all the indices correctly
		which.position.upper.marker <- which.pixel.upper.marker[,x.gel]
		which.position.lower.marker <- which.pixel.lower.marker[,x.gel]
		stopifnot(apply(which.position.upper.marker, 2, any))
		stopifnot(apply(which.position.lower.marker, 2, any))

		# find the center y-value of the marker in each gel lane (the lanes might not be aligned perfectly at both ends)
		center.upper.marker <- apply(which.position.upper.marker, 2, function(lane) mean(which(lane)))
		center.lower.marker <- apply(which.position.lower.marker, 2, function(lane) mean(which(lane)))

		# find which positions safely contain no residual color from the markers (the edges are blurred and we won't try to deconvolve the color from the signal)
		which.position.not.marker <- apply(gel.image.rgb[,x.gel,], 1:2, function(channels) channels[1] == channels[2] && channels[2] == channels[3])

		# now finally extract the intensities!
		position.fluorescence <- 1 - gel.image.rgb[,x.gel,1] # only get red fluorescence because all channels are equal in the places we care about; subtract from 1 because it's a negative
		results <- do.call(rbind, lapply(1:length(x.gel), function(gel.lane) {
			positions.to.use <- which(
				1:nrow(position.fluorescence) < center.lower.marker[gel.lane] &
				1:nrow(position.fluorescence) > center.upper.marker[gel.lane] &
				which.position.not.marker[,gel.lane]
			)
			relative.distance <- (center.upper.marker[gel.lane] - positions.to.use) / (center.upper.marker[gel.lane] - center.lower.marker[gel.lane])
			data.frame(gel.lane, relative.distance, fluorescence = position.fluorescence[positions.to.use, gel.lane])	
		}))
		
		cbind(batch = sub("\\.png$", "", basename(gel.image.file)), results)
	}))
	
	rownames(combined.results) <- NULL
	combined.results$batch <- factor(combined.results$batch, levels = unique(combined.results$batch)) # make batches into a factor that keeps them in the observed order
	
	combined.results
}


read.tapestation.xml <- function(xml.files) {
	combined.results <- do.call(rbind, lapply(xml.files, function(xml.file) {
		xml.root <- xmlRoot(xmlParse(xml.file))
		results <- do.call(rbind, xmlApply(xml.root[["Samples"]], function(sample.xml) {
			well.number <- xmlValue(sample.xml[["WellNumber"]])
			name <- trimws(xmlValue(sample.xml[["Comment"]]))
			sample.observations <- trimws(xmlValue(sample.xml[["Observations"]]))
			if (sample.observations == "Marker(s) not detected") {
				warning(paste(sample.observations, "for well", well.number, name))
				return(NULL)
			}
			well.row <- gsub("[^A-H]", "", well.number) # assumes all row names are in A through H!
			well.col <- as.integer(gsub("[^[:digit:]]", "", well.number)) # assumes all column names are numbers!
			
			reagent.id <- xmlValue(sample.xml[["ScreenTapeID"]])
			
			peaks <- do.call(rbind, xmlApply(sample.xml[["Peaks"]], function(peak.xml) data.frame(
				peak.observations =  trimws(xmlValue(peak.xml[["Observations"]])),
				length =             as.integer(xmlValue(peak.xml[["Size"]])),
				distance =           as.numeric(xmlValue(peak.xml[["RunDistance"]])),
				lower.bound =        as.numeric(xmlValue(peak.xml[["FromPercent"]])),
				upper.bound =        as.numeric(xmlValue(peak.xml[["ToPercent"]])), 
				area =               as.numeric(xmlValue(peak.xml[["Area"]])),
				molarity =           as.numeric(xmlValue(peak.xml[["Molarity"]])
			))))
			peaks$relative.distance <- (peaks$distance - min(peaks$distance)) / (max(peaks$distance) - min(peaks$distance))
			
			cbind(well.number, well.row, well.col, name, reagent.id, sample.observations, peaks)
		}))
		
		cbind(batch = sub("\\.xml$", "", basename(xml.file)), results)
	}))
	
	rownames(combined.results) <- NULL
	combined.results$batch <- factor(combined.results$batch, levels = unique(combined.results$batch)) # make batches into a factor that keeps them in the observed order
	
	combined.results
}


read.tapestation <- function(xml.file, gel.image.file, fit = "spline") {
	stopifnot(fit %in% c("interpolate", "spline", "regression"))
	
	peaks <- read.tapestation.xml(xml.file)
	result <- read.tapestation.gel.image(gel.image.file)
	stopifnot(length(unique(result$gel.lane)) == length(unique(peaks$well.number)))
	
	# convert relative distances to absolute
	marker.distances <- data.frame(
		lower = subset(peaks, peak.observations == "Lower Marker")$distance,
		upper = subset(peaks, peak.observations == "Upper Marker")$distance,
		row.names = unique(peaks$well.number)
	)
	result$distance <- marker.distances$upper[result$gel.lane] + result$relative.distance * (marker.distances$lower[result$gel.lane] - marker.distances$upper[result$gel.lane])
	
	# fit standard curve for molecule length
	# do this in relative-distance space so it's effectively recalibrated for each sample's markers
	which.well.is.ladder <- unique(subset(peaks, sample.observations == "Ladder")$well.number)
	stopifnot(length(which.well.is.ladder) == 1)
	peaks.ladder <- subset(peaks, well.number == which.well.is.ladder)
	peaks.ladder$relative.distance <- (peaks.ladder$distance - subset(peaks.ladder, peak.observations == "Upper Marker")$distance) / (subset(peaks.ladder, peak.observations == "Lower Marker")$distance - subset(peaks.ladder, peak.observations == "Upper Marker")$distance)
	if (fit == "interpolate") {
		warning("linear interpolation gives ugly results for molarity estimation")
		standard.curve.function <- approxfun(peaks.ladder$relative.distance, peaks.ladder$length)
	} else if (fit == "spline") {
		standard.curve.function <- splinefun(peaks.ladder$relative.distance, peaks.ladder$length, method = "natural")
	} else if (fit == "regression") {
		mobility.model <- lm(relative.distance ~ log(length), peaks.ladder)
		standard.curve.function <- function(distance) exp((distance - mobility.model$coefficients[1]) / mobility.model$coefficients[2])
	}
	result$length <- standard.curve.function(result$relative.distance)
	
	# convert to molarity
	fluorescence.coefficient <- 1 / lm(area ~ length : molarity - 1, peaks.ladder)$coefficients
	result$molarity <- unlist(lapply(unique(result$gel.lane), function(this.well) {
		this.result <- subset(result, gel.lane == this.well)
		length.derivative <- diff(this.result$length) / diff(this.result$relative.distance)
		this.result$fluorescence * fluorescence.coefficient / this.result$relative.distance / c(NA, -length.derivative)
	}))
	
	# get sample names (might not be unique but we assume well numbers are)
	sample.table <- unique(peaks[,c("batch", "name", "well.number", "well.row", "well.col")])
	stopifnot(nrow(sample.table) == length(unique(result$gel.lane)))
	
	result <- cbind(sample.table[result$gel.lane,], subset(result, select = -batch))
	rownames(result) <- NULL # clean up row names again
	
	result
}

