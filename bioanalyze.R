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

fit.standard.curve <- function(ladder, ...) loess(ladder[,"Size [bp]"] ~ ladder[,"Aligned Migration Time [s]"], ...) # returns an object of type "loess"; use it e.g. with predict()


