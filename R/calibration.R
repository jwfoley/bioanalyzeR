#' Estimate molecular weight of nucleic acids by length
#'
#' This function estimates the molecular weight of a single-stranded RNA or double-stranded DNA from its length in bases. 
#'
#' @param length The length of the nucleic acid in bases (nt for RNA, bp for DNA).
#' @param type The type of nucleic acid: \code{"RNA"} for single-stranded RNA, \code{"DNA"} for double-stranded DNA.
#'
#' @return The estimated molecular weight in daltons.
#'
#' @references
#' Thermo Fisher Scientific: DNA and RNA Molecular Weights and Conversions. Accessed 13 March 2020. \url{https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-tools-and-calculators/dna-and-rna-molecular-weights-and-conversions.html}
#'
#' @export
molecular.weight <- function(length, type) switch(type,
	DNA = length * 607.4 + 157.9,
	RNA = length * 320.5 + 159.0
)


#' Electrophoresis calibrations
#'
#' Estimate the derived values from raw electrophoresis data.
#'
#' Ladder peaks of known length and possibly concentration are used to fit mobility models and coefficients of fluorescence area vs. concentration. If there are multiple ladders (e.g. one per ScreenTape), each sample is fit according to its corresponding ladder. Derived values are estimated for every row in \code{electrophoresis$data}, though they are \code{NA} outside the range of interpolation from the ladder.
#'
#' Spline fitting seems to perform reasonably well on all data, and among R's built-in methods Hyman filtering of the FMM spline (\code{method = "hyman"}) is preferred to ensure monotonicity. Agilent appears to use linear interpolation with DNA data and log-linear regression on RNA data, so you could choose those options if you want to reproduce the results of the software more precisely. However, linear interpolation creates sudden spikes in the derivative that make the concentration and molarity estimates unstable; spline fitting is basically a smoother version of that. Log-linear regression is the standard theoretical approach but does not actually fit the data very well; more sophisticated parametric models may be added in the future.
#'
#' @param electrophoresis An \code{\link{electrophoresis}} object.
#' @param method The method used to fit the mobility model of molecule length vs. migration distance, either \code{"interpolation"} (linear interpolation via \code{\link{approxfun}}), \code{"loglinear"} (log-linear regression via \code{\link{lm}} with the model \code{relative.distance ~ log(length)}), or one of the methods of \code{\link{splinefun}} (monotone methods recommended). Can be abbreviated.
#' @param ladder.concentrations The true concentrations of the ladder peaks. If provided (from the Bioanalyzer), the concentration coefficient is fit according to the non-marker ladder peaks and then adjusted for each sample according to the relative fluorescence area of its markers compared with the ladder. If \code{NULL} (TapeStation), concentrations are fit according to only the upper marker if present or the lower marker otherwise; these marker concentrations are assumed to be correct.
#'
#' @return The same \code{\link{electrophoresis}} object with the new estimated variable added to its \code{$data} member. \code{calculate.length} also adds estimated length boundaries to the \code{$peaks} member and aligned-time or relative-distance boundaries to \code{$regions}.
#'
#' @name calibrate.electrophoresis
NULL


#' @rdname calibrate.electrophoresis
#' @export
calculate.length <- function(electrophoresis, method = union(c("hyman", "interpolation", "loglinear"), eval(formals(splinefun)$method)), extrapolate = FALSE) {
	method <- match.arg(method)
	x.name <- get.x.name(electrophoresis)
	lower.name <- paste0("lower.", x.name)
	upper.name <- paste0("upper.", x.name)
		
	# prepare new fields to be filled in piecemeal from each ladder
	electrophoresis$data$length <- NA
	electrophoresis$peaks$lower.length <- NA
	electrophoresis$peaks$upper.length <- NA
	if (! is.null(electrophoresis$regions)) {
		electrophoresis$regions[[lower.name]] <- NA
		electrophoresis$regions[[upper.name]] <- NA
	}
	electrophoresis$calibration <- list()
	
	# fit a mobility model for each ladder and apply it to the appropriate samples
	for (batch in unique(electrophoresis$samples$batch)) {
		electrophoresis$calibration[[batch]] <- list()
		in.this.batch <- electrophoresis$samples$batch == batch
		for (ladder.well in unique(electrophoresis$samples$ladder.well[which(in.this.batch)])) {
			which.ladder.index <- which(in.this.batch & electrophoresis$samples$well.number == ladder.well)
			peaks.ladder <- subset(electrophoresis$peaks, sample.index == which.ladder.index)
			which.samples <- which(in.this.batch & electrophoresis$samples$ladder.well == ladder.well)
			which.rows <- which(electrophoresis$data$sample.index %in% which.samples)
			which.peaks <- which(electrophoresis$peaks$sample.index %in% which.samples)
			which.regions <- which(electrophoresis$regions$sample.index %in% which.samples)
			
			peaks.ladder$x <- peaks.ladder[[x.name]]
			
			# fit standard curve for molecule length vs. x-value
			# do this in relative x space so it's effectively recalibrated for each sample's markers
			if (method == "interpolation") {
				warning("linear interpolation gives ugly results for molarity estimation")
				standard.curve.function <- approxfun(peaks.ladder$x, peaks.ladder$length)
				standard.curve.inverse <- approxfun(peaks.ladder$length, peaks.ladder$x)
				
			} else if (method == "loglinear") {
				if (x.name == "relative.distance") {
					mobility.model <- lm(x ~ log(length), peaks.ladder)
					coefs <- coefficients(mobility.model)
					standard.curve.function <- function(relative.distance) exp((relative.distance - coefs[1]) / coefs[2])
					standard.curve.inverse <- function(length) coefs[1] + coefs[2] * log(length)
				} else if (x.name == "aligned.time") {
					# if x-variable is time, we must correct for the fact that faster-moving molecules spend less time in front of the detector
					mobility.model <- lm(1/aligned.time ~ log(length), data = peaks.ladder)
					coefs <- coefficients(mobility.model)
					standard.curve.function <- function(aligned.time) exp((1 / aligned.time - coefs[1]) / coefs[2])
					standard.curve.inverse <- function(length) 1/(coefs[1] + log(length) * coefs[2])
				}
				
			} else { # if it's not one of those then it must be one of the splinefun methods
				standard.curve.function <- splinefun(peaks.ladder$x, peaks.ladder$length, method = method)
				standard.curve.inverse <- splinefun(peaks.ladder$length, peaks.ladder$x, method = method)
			}
			electrophoresis$calibration[[batch]][[ladder.well]] <- list(
				mobility.function = standard.curve.function,
				mobility.inverse = standard.curve.inverse,
				ladder.peaks = setNames(data.frame(peaks.ladder$length, peaks.ladder$x), c("length", x.name))
			)
			
			# apply model to raw data
			electrophoresis$data$length[which.rows] <- standard.curve.function(electrophoresis$data[[x.name]][which.rows])
			if (! extrapolate) { # don't extrapolate outside ladder range
				electrophoresis$data$length[! (
					in.custom.region(electrophoresis$data, min(peaks.ladder$length), max(peaks.ladder$length)) &
					in.custom.region(electrophoresis$data, min(peaks.ladder$x), max(peaks.ladder$x), bound.variable = x.name)
				)] <- NA
			} else { # just don't extrapolate below zero length (sometimes the function will have weird non-monotonic behavior there so we have to discard everything up to the last estimate below 0)
				if (method != "loglinear") warning("mobility model was fit by interpolation so there is no expectation of accuracy outside the ladder range")
				for (i in unique(electrophoresis$data$sample.index)) {
					which.negative <- which(electrophoresis$data$sample.index == i & electrophoresis$data$length < 0)
					if (any(which.negative, na.rm = T)) electrophoresis$data$length[seq(
						which(electrophoresis$data$sample.index == i)[1],
						tail(which.negative, 1)
					)] <- NA
				}
			}
			
			# before applying model to other annotations, lower and upper depend on the x.name (inverse directions)
			if (x.name == "relative.distance") {
				lower.length.analog <- upper.name
				upper.length.analog <- lower.name
			} else if (x.name == "aligned.time") {
				lower.length.analog <- lower.name
				upper.length.analog <- upper.name
			}
			
			# apply model to peaks
			electrophoresis$peaks$lower.length[which.peaks] <- standard.curve.function(electrophoresis$peaks[[lower.length.analog]][which.peaks])
			electrophoresis$peaks$upper.length[which.peaks] <- standard.curve.function(electrophoresis$peaks[[upper.length.analog]][which.peaks])
			
			# apply inverse model to regions
			if (! is.null(electrophoresis$regions)) {
				electrophoresis$regions[[lower.length.analog]][which.regions] <- standard.curve.inverse(electrophoresis$regions$lower.length[which.regions])
				electrophoresis$regions[[upper.length.analog]][which.regions] <- standard.curve.inverse(electrophoresis$regions$upper.length[which.regions])
			}
		}
		electrophoresis$assay.info[[batch]]$method <- method
	}
	
	electrophoresis
}


#' @rdname calibrate.electrophoresis
#' @export
calculate.concentration <- function(electrophoresis, ladder.concentrations = NULL) {
	x.name <- get.x.name(electrophoresis)
	delta <- do.call(rbind, by(electrophoresis$data, electrophoresis$data$sample.index, function(data.subset) data.frame(
		fluorescence = c(NA, diff(data.subset$fluorescence)),
		x = c(NA, diff(data.subset[[x.name]])) # these will probably be constant but we'd better not assume
	), simplify = F))
	if (x.name == "relative.distance") delta$x <- -delta$x # distances are stored in decreasing order so the deltas are negative
	# estimate area under each measurement with the trapezoidal rule; to simplify math, each point's sum is for the trapezoid to the left of it
	electrophoresis$data$area <- (2 * electrophoresis$data$fluorescence - delta$fluorescence) * delta$x
	if (x.name == "aligned.time") electrophoresis$data$area <- electrophoresis$data$area / electrophoresis$data[[x.name]] # compensate for time spent in the detector (faster-moving molecules get less signal)
	
	has.upper.marker <- any(electrophoresis$peaks$peak.observations %in% UPPER.MARKER.NAMES)
	which.markers <- lapply(seq(nrow(electrophoresis$samples)), function(sample.index) which(
		electrophoresis$peaks$sample.index == sample.index & (
			(electrophoresis$peaks$peak.observations %in% LOWER.MARKER.NAMES & (
				if (is.null(ladder.concentrations)) TRUE else electrophoresis$peaks$concentration == ladder.concentrations[1] # verify this is the right peak (sometimes Bioanalyzer annotates more than one as the marker but it only assigns the hardcoded concentration to one)
				)
			) | (
				electrophoresis$peaks$peak.observations %in% UPPER.MARKER.NAMES & (
					if (is.null(ladder.concentrations)) TRUE else electrophoresis$peaks$concentration == ladder.concentrations[length(ladder.concentrations)]
				)
			)
		)
	))
	for (i in seq(nrow(electrophoresis$samples))) if (length(which.markers[[i]]) == 0) warning(paste("no markers detected in sample", i))
	
	mass.coefficients <- rep(NA, nrow(electrophoresis$samples))
	for (batch in unique(electrophoresis$samples$batch)) {
		in.this.batch <- electrophoresis$samples$batch == batch
		for (ladder.well in unique(electrophoresis$samples$ladder.well[which(in.this.batch)])) {
			ladder.index <- which(in.this.batch & electrophoresis$samples$well.number == ladder.well)
			stopifnot("conflicting ladders" = length(ladder.index) == 1)
			which.samples <- which(in.this.batch & electrophoresis$samples$ladder.well == ladder.well)
			this.ladder.concentrations <- if (! is.null(ladder.concentrations)) ladder.concentrations else subset(electrophoresis$peaks, sample.index == ladder.index)$concentration # if true concentrations aren't provided, use the reported ones in the peak table
			non.marker.concentrations <- this.ladder.concentrations[-c(1, if (has.upper.marker) length(this.ladder.concentrations) else NULL)]
			ladder.peaks <- which(
				electrophoresis$peaks$sample.index == ladder.index &
				(if (is.null(ladder.concentrations)) TRUE else electrophoresis$peaks$concentration %in% ladder.concentrations)
			)
			stopifnot("no ladder peaks recognized" = length(ladder.peaks) > 0)
			non.marker.peaks <- setdiff(ladder.peaks, which.markers[[ladder.index]])
			ladder.areas <- integrate.peak(electrophoresis, ladder.peaks, "area")
			non.marker.areas <- ladder.areas[! ladder.peaks %in% which.markers[[ladder.index]]]
			ladder.marker.areas <- ladder.areas[ladder.peaks %in% which.markers[[ladder.index]]]
			
			ladder.mass.coefficient <- mean(electrophoresis$peaks$concentration[non.marker.peaks] / non.marker.areas, na.rm = T)
			
			# then modify the mass coefficient for each sample according to the ratio of its marker area(s) to the ladder's (compensate for differential fluorescence/detection/loading)
			for (sample.index in which.samples) {
				sample.marker.areas <- integrate.peak(electrophoresis, which.markers[[sample.index]], "area")
				mass.coefficients[sample.index] <- if (length(sample.marker.areas) == 0) NA else ladder.mass.coefficient * mean(ladder.marker.areas / sample.marker.areas, na.rm = T)
			}
		}
	}
	
	electrophoresis$data$concentration <- mass.coefficients[electrophoresis$data$sample.index] * electrophoresis$data$area	
	electrophoresis$data$area <- NULL # don't need this anymore so keep the structure clean
	electrophoresis
}


#' @rdname calibrate.electrophoresis
#' @export
calculate.molarity <- function(electrophoresis) {
	electrophoresis$data$molarity <- NA
	for (batch in unique(electrophoresis$samples$batch)) {
		which.rows <- which(electrophoresis$data$sample.index %in% which(electrophoresis$samples$batch == batch))
		electrophoresis$data$molarity[which.rows] <- electrophoresis$data$concentration[which.rows] / molecular.weight(electrophoresis$data$length[which.rows], electrophoresis$assay.info[[batch]]$assay.type) * 1E6 # we're converting ng/uL to nmol/L or pg/uL to pmol/L so we need to scale by 1E6
	}
	
	electrophoresis
}

