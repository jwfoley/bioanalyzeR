library(png)

color.upper.marker <- rgb(128, 0, 128, maxColorValue = 255) # upper marker is purple
color.lower.marker <- rgb(0, 128, 0, maxColorValue = 255) # lower marker is green


read.gel.image <- function(gel.image.file) {
	gel.image.rgb <- readPNG(gel.image.file) # this is in the form (y, x, channel); [1,1,] is the upper left corner
	gel.image.hex <- apply(gel.image.rgb, 1:2, function(pixel) rgb(pixel[1], pixel[2], pixel[3])) # now it's a 2D matrix of hexadecimal values; this takes a long time to convert but saves more time on other computations later

	# find marker bands
	which.pixel.upper.marker <- gel.image.hex == color.upper.marker
	which.pixel.lower.marker <- gel.image.hex == color.lower.marker

	# use marker bands to identify gel lanes
	which.x.upper.marker <- apply(which.pixel.upper.marker, 2, any)
	which.x.lower.marker <- apply(which.pixel.lower.marker, 2, any)
	stopifnot(which.x.upper.marker == which.x.lower.marker)

	# isolate a single representative x-value (the first one from the left) for each gel lane
	x.gel <- which(diff(which.x.upper.marker) == 1) + 1
	gel.lanes.raw <- gel.image.hex[,x.gel]
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
	intensities <- do.call(rbind, lapply(1:length(x.gel), function(lane) {
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
	
	intensities
}

