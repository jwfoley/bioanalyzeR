library(png)

rgb.upper.marker <- c(128, 0, 128) / 255 # upper marker is purple
rgb.lower.marker <- c(0, 128, 0) / 255 # lower marker is green


read.gel.image <- function(gel.image.files) {

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
		x.gel <- which(diff(which.x.upper.marker) == 1) + 1
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
		position.intensity <- 1 - gel.image.rgb[,x.gel,1] # only get red intensity because all channels are equal in the places we care about; subtract from 1 because it's a negative
		results <- do.call(rbind, lapply(1:length(x.gel), function(lane) {
			positions.to.use <- which(
				1:nrow(position.intensity) < center.lower.marker[lane] &
				1:nrow(position.intensity) > center.upper.marker[lane] &
				which.position.not.marker[,lane]
			)
			distance <- (positions.to.use - center.lower.marker[lane]) / (center.upper.marker[lane] - center.lower.marker[lane])
			data.frame(
				lane = lane,
				distance = distance,
				intensity = position.intensity[positions.to.use, lane]
			)	
		}))
		
		cbind(batch = sub("\\.png$", "", basename(gel.image.file)), results)
	}))
	
	rownames(combined.results) <- NULL
	combined.results$batch <- factor(combined.results$batch, levels = unique(combined.results$batch)) # make batches into a factor that keeps them in the observed order
	
	combined.results
}

