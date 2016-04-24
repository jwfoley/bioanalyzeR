get.ladder <- function(results.file, ...) {
	
	# assume file contains something in the form:
	#Sample Name,Ladder
	#   
	#Peak Table
	#[headers]
	#[data]
	#[data]
	#...
	# 
	# and then returns it as a data frame
	
	file.lines <- readLines(results.file)
	
	# find the ladder table
	which.ladder <- which(file.lines == "Sample Name,Ladder")
	stopifnot(
		length(which.ladder) == 1 &&
		file.lines[which.ladder + 1] == " " &&
		file.lines[which.ladder + 2] == "Peak Table"
	)
	which.breaks <- which(file.lines == " ")
	
	return(read.csv(text = 
		gsub("\"(\\d+),([\\d.]+)\"", "\\1\\2", file.lines[ 	# sanitize the ladder table (remove thousands separators and the quotation marks around them)
			(which.ladder + 3):
			(min(which.breaks[which.breaks > which.ladder + 1]) - 1)
		], perl = T),  check.names = F, ...)
	)
}

get.sample <- function(data.file, ...) {

	# assume the fluorescence vs. time data are in the last block

	file.lines = readLines(data.file)
	which.breaks <- which(file.lines == " ")
	return(read.csv(text = file.lines[(
		which.breaks[length(which.breaks) - 1] + 1):
		(which.breaks[length(which.breaks)] - 1)
	], ...))
}

fit.standard.curve <- function(ladder, ...) loess(log10(ladder$Size) ~ 1 / ladder$`Aligned Migration Time [s]`, ...) # returns an object of type "loess"; use it e.g. with predict() but make sure to take 10^predict()
# (1 / migration time) is proportional to velocity, which is what should be logarithmic with molecule size

time.to.bp <- function(standard.curve, time, ...) 10 ^ predict(standard.curve, time, ...)

plot.raw <- function(sample.data, ...) qplot(sample.data$Time, sample.data$Value) + geom_line()

plot.bp <- function(sample.data, standard.curve, ...) qplot(10 ^ predict(standard.curve, sample.data$Time), sample.data$Value) + geom_line()

normalize.fluorescence <- function(ladder) ladder$Area / ladder$`Size [bp]` / ladder$`Aligned Migration Time [s]` # gets a value that is proportional to molarity
# because area under the peak is total fluorescence during a time window, so we normalize by the size to get to molecules instead of mass, and normalize by migration time because it is proportional to velocity (to account for the time the molecule spends in the detector)

get.fluorescence.coefficient <- function(ladder, ...) lm(ladder$`Molarity [pmol/l]` ~ normalize.fluorescence(ladder) - 1)$coefficients # finds the coefficient to convert normalize.fluorescence to molarity (over an area) by fitting a linear model to the standards of known molarities

convert.to.molarity <- function(ladder, sample.data, ...) {
	std.crv <- fit.standard.curve(ladder)
	bp <- time.to.bp(std.crv, sample.data$Time)
	fluorescence.coefficient <- get.fluorescence.coefficient(ladder)
	std.crv.deriv <- diff(bp) / diff(sample.data$Time) 
	molarity <- sample.data$Value * fluorescence.coefficient / sample.data$Time / c(NA, std.crv.deriv) / bp
	# fluorescence times coefficient isn't corrected for area
	# divide by time to compensate for time spent in detector
	# divide by derivative of standard curve because we're converting from an area of the fluorescence vs. time graph to an area of the fluorescence vs. bp graph
	molarity
}

