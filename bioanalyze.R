library(XML)
library(openssl)

read.bioanalyzer <- function(xml.files, fit = "spline") {
	stopifnot(fit %in% c("linear", "spline", "regression"))

	results <- do.call(rbind, lapply(xml.files, function(xml.file) {
	
		xml.root <- xmlRoot(xmlParse(xml.file))
		
		assay.definitions <- xml.root[["Chips"]][["Chip"]][["AssayBody"]][["DAAssaySetpoints"]][["DAMAssayInfoMolecular"]]
		ladder.definition <- data.frame(t(xmlSApply(assay.definitions[["LadderPeaks"]], function(standard) sapply(c("Size", "Concentration"), function (field) as.numeric(xmlValue(standard[[field]]))))), row.names = NULL)
		marker.aligned.times <- sapply(c("VirtualLowerMarkerTime", "VirtualUpperMarkerTime"), function (field) as.numeric(xmlValue(assay.definitions[[field]])))
		
		# read raw sample data
		samples <- xml.root[["Chips"]][["Chip"]][["Files"]][["File"]][["Samples"]]
		result <- do.call(rbind, xmlApply(samples, function(this.sample) {
			if (xmlValue(this.sample[["HasData"]]) != "true") NULL else {
				signal.data <- this.sample[["DASignals"]][["DetectorChannels"]][[1]][["SignalData"]]
				n.values <- as.integer(xmlValue(signal.data[["NumberOfSamples"]]))
				raw.data <- data.frame(
					well.number = as.integer(xmlValue(this.sample[["WellNumber"]])),
					name = xmlValue(this.sample[["Name"]]),
					time = as.numeric(xmlValue(signal.data[["XStart"]])) + as.numeric(xmlValue(signal.data[["XStep"]])) * (1:n.values - 1),
					fluorescence = readBin(base64_decode(xmlValue(signal.data[["ProcessedSignal"]])), "numeric", size = 4, n = n.values, endian = "little"),
					stringsAsFactors = F
				)
				
				# align the observation times according to the markers in this sample
				peaks.sample <- this.sample[["DAResultStructures"]][["DARIntegrator"]][["Channel"]][["PeaksMolecular"]]
				which.lower.marker <- which(xmlSApply(peaks.sample, function(peak) xmlValue(peak[["Observations"]]) == "Lower Marker"))
				stopifnot(length(which.lower.marker) == 1)
				stopifnot(as.numeric(xmlValue(peaks.sample[[which.lower.marker]][["AlignedMigrationTime"]])) == marker.aligned.times[1])
				which.upper.marker <- which(xmlSApply(peaks.sample, function(peak) xmlValue(peak[["Observations"]]) == "Upper Marker"))
				stopifnot(length(which.upper.marker) == 1)
				stopifnot(as.numeric(xmlValue(peaks.sample[[which.upper.marker]][["AlignedMigrationTime"]])) == marker.aligned.times[2])
				alignment.coefficient <- diff(marker.aligned.times) / diff(sapply(c(which.lower.marker, which.upper.marker), function(i) as.numeric(xmlValue(peaks.sample[[i]][["MigrationTime"]]))))
				alignment.offset <- marker.aligned.times[1] - alignment.coefficient * as.numeric(xmlValue(peaks.sample[[which.lower.marker]][["MigrationTime"]]))
				raw.data$aligned.time <- raw.data$time * alignment.coefficient + alignment.offset
				
				raw.data
			}
		}))
		
		# analyze ladder
		which.ladder <- which(xmlApply(samples, function(this.sample) xmlValue(this.sample[["Category"]])) == "Ladder")
		stopifnot(length(which.ladder) == 1)
		sample.ladder <- samples[[which.ladder]]
		stopifnot(xmlValue(sample.ladder[["HasData"]]) == "true")
		peaks.ladder <- data.frame(t(xmlSApply(sample.ladder[["DAResultStructures"]][["DARIntegrator"]][["Channel"]][["PeaksMolecular"]], function(peak) sapply(c("AlignedMigrationTime", "Area", "FragmentSize", "Molarity", "Concentration", "StartTime", "EndTime", "AlignedStartTime", "AlignedEndTime"), function(field) as.numeric(xmlValue(peak[[field]]))))), row.names = NULL)
		stopifnot(all.equal(ladder.definition, peaks.ladder[,c("FragmentSize", "Concentration")], check.names = F))
		
		# fit standard curve for molecule length
		# do this with aligned times so it's effectively recalibrated for each sample's markers
		if (fit == "interpolate") {
			warning("linear interpolation gives ugly results for molarity estimation")
			standard.curve.function <- approxfun(peaks.ladder$AlignedMigrationTime, peaks.ladder$FragmentSize)
		} else if (fit == "spline") {
			standard.curve.function <- splinefun(peaks.ladder$AlignedMigrationTime, peaks.ladder$FragmentSize, method = "natural")
		} else if (fit == "regression") {
			mobility.model <- lm(1/AlignedMigrationTime ~ log(FragmentSize), data = peaks.ladder)
			standard.curve.function <- function(time) exp((1 / time - mobility.model$coefficients[1]) / mobility.model$coefficients[2])
		}
		result$length <- sapply(result$aligned.time, function(x) if (x < marker.aligned.times[1] || x > marker.aligned.times[2]) NA else  standard.curve.function(x)) # avoid extrapolating
		
		# convert to molarity
		peaks.ladder$mass <- peaks.ladder$FragmentSize * peaks.ladder$Molarity
		peaks.ladder$mass.times.length <- peaks.ladder$mass * peaks.ladder$FragmentSize
		result <- cbind(result, do.call(rbind, lapply(unique(result$well.number), function(this.well) {
			result.this.well <- subset(result, well.number == this.well)
			data.frame(
				delta.time = c(NA, diff(result.this.well$time)),
				delta.fluorescence = c(NA, diff(result.this.well$fluorescence))
			)
		})))
		result$delta.area <- (2 * result$fluorescence - result$delta.fluorescence) / 2 * result$delta.time
		peaks.ladder$estimated.area <- sapply(1:nrow(peaks.ladder), function(i) {
			sum(subset(result, well.number == which.ladder & time >= peaks.ladder$StartTime[i] & time <= peaks.ladder$EndTime[i])$delta.area)
		})
		# problem: these areas don't match the ones in the XML (except the second ladder band) and they're not off by a constant either
		
#		
#	
#				# fit standard curve for molecule length (done once per sample because that's how the results are provided, even though it should be the same standard curve for every sample)
#				ladder.peaks <- this.sample[["DAResultStructures"]][["DARStandardCurve"]][[1]][["LadderPeaks"]]
#				stopifnot(all.equal(xmlSApply(ladder.peaks, function(peak) as.numeric(xmlValue(peak[["FragmentSize"]]))), ladder$length, check.names = F)) # verify it matches the known ladder
#				this.ladder <- cbind(ladder, data.frame(t(xmlSApply(ladder.peaks, function(peak) c(as.numeric(xmlValue(peak[["MigrationTime"]])), as.numeric(xmlValue(peak[["Area"]]))))), row.names = NULL))
#				names(this.ladder)[ncol(ladder) + 1:2] <- c("time", "area")
#				
#				if (fit == "linear") {
#					standard.curve.function <- approxfun(this.ladder$time, this.ladder$length)		
#				} else if (fit == "spline") {
#					standard.curve.function <- splinefun(this.ladder$time, this.ladder$length, method = "natural")
#				}
#				result$length <- standard.curve.function(result$time)
#				
#				# convert to molarity
#				this.ladder$normalized.fluorescence <- this.ladder$area / this.ladder$length / this.ladder$time
#				fluorescence.coefficient <- lm(molarity ~ normalized.fluorescence - 1, this.ladder)$coefficients
#				length.derivative <- diff(result$length) / diff(result$time)
#				result$molarity <- result$fluorescence * fluorescence.coefficient / result$time / c(NA, length.derivative) / result$length
#			
#				result
#			}
#		}))
		
		cbind(batch = sub("\\.xml$", "", basename(xml.file)), result)
	}))
	
	# format output nicely
	rownames(results) <- NULL
	results$batch <- factor(results$batch, levels = unique(results$batch)) # make batches into a factor that keeps them in the observed order
	results$name <- factor(results$name, levels = unique(results$name)) # make names into a factor that keeps them in the observed order
	
	results
}

