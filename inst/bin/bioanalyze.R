#! /usr/bin/Rscript

library(parallel)
library(argparse)
library(bioanalyzeR)

parser <- ArgumentParser(description = "Simple automation of bioanalyzeR functions.")
parser$add_argument("--version", "-v",
	action = "version",
	version = as.character(packageVersion("bioanalyzeR"))
)

files <- parser$add_argument_group("files", "Filenames and settings for input/output.")
files$add_argument("xml_files",
	nargs = "+",
	help = "one or more XML files exported from the Bioanalyzer or TapeStation software",
	metavar = "INPUT_FILE.xml[.gz]"
)
files$add_argument("--output_file", "-o",
	help = "path of file to write sample table (TSV)",
	metavar = "OUTPUT_FILE.tsv"
)
files$add_argument("--plot_file", "-p",
	help = "path of file to write plot(s) (PDF)",
	metavar = "PLOT_FILE.pdf"
)
files$add_argument("--dimensions", "-d",
	help = "dimensions of plot(s) (in.)",
	metavar = "WxH",
	default = "8x6"
)
files$add_argument("--method",
	help = "mobility standard curve method",
	default = "hyman",
	choices = eval(formals(calculate.length)$method)
)
files$add_argument("--mc_cores",
	help = "maximum CPU cores",
	default = if (.Platform$OS.type == "windows") 1 else detectCores(),
	type = "integer",
	metavar = "N"
)

annotations <- parser$add_argument_group("annotations", "Add sample annotations from a table. See help(annotate.electrophoresis).")
annotations$add_argument("--annotation_file", "-a",
	help = "annotation file",
	metavar = "ANNOTATIONS.csv"
)
annotations$add_argument("--annot_args",
	help = "string of arguments for read.table, e.g. 'sep = \",\"' for comma-separated values",
	metavar = "ARG_STRING"
)

plotting <- parser$add_argument_group("plotting", "Settings for the electropherogram plot. See help(qplot.electrophoresis).")
plotting$add_argument("-x", help = "x-variable",
	default = "length",
	metavar = "VARIABLE"
)
plotting$add_argument("-y", help = "y-variable",
	default = "molarity",
	metavar = "VARIABLE"
)
plotting$add_argument("--log", "-l",
	help = "axes to log-scale",
	default = "",
	choices = c("x", "y", "xy")
)
plotting$add_argument("--normalize", "-n",
	help = "normalize y-axis per sample",
	action = "store_true"
)
plotting$add_argument("--facets", "-f",
	help = "faceting formula, or 'none'",
	default = "~ sample.index",
	metavar = "FORMULA"
)
plotting$add_argument("--margins",
	help = "display marginal facets",
	action = "store_true"
)
plotting$add_argument("--scales", "-s",
	help = "scaling rules for facets",
	default = "fixed",
	choices = c("fixed", "free_x", "free_y", "free")
)
plotting$add_argument("--geom", "-g",
	help = "geom to draw",
	default = "line",
	choices = c("line", "area")
)
plotting$add_argument("--include_ladder",
	help = "include ladder(s)",
	action = "store_true"
)
plotting$add_argument("--include_markers",
	help = "extend range to include marker(s)",
	action = "store_true"
)
plotting$add_argument("--lower_marker_spread",
	help = "amount to extend lower marker",
	default = 5,
	type = "double",
	metavar = "PROPORTION"
)
plotting$add_argument("--xlim",
	help = "limits of x-axis",
	type = "double",
	nargs = 2,
	metavar = c("MIN", "MAX")
)
plotting$add_argument("--ylim",
	help = "limits of y-axis",
	type = "double",
	nargs = 2,
	metavar = c("MIN", "MAX")
)
plotting$add_argument("--hide_peaks",
	help = "don't fill the area under peaks",
	action = "store_true"
)
plotting$add_argument("--region_alpha",
	help = "alpha of region highlights",
	default = 0.2,
	type = "double",
	metavar = "ALPHA"
)
plotting$add_argument("--area_alpha",
	help = "alpha of filled areas in unfaceted area geom",
	default = 0.2,
	type = "double",
	metavar = "ALPHA"
)
plotting$add_argument("--title",
	help = "plot title",
	metavar = "LABEL"
)
plotting$add_argument("--xlab",
	help = "x-axis label",
	metavar = "LABEL"
)
plotting$add_argument("--ylab",
	help = "y-axis label",
	metavar = "LABEL"
)

integration <- parser$add_argument_group("integration", "Integrate a variable under the curve. See help(integrate.custom) and help(region.ratio).")
integration$add_argument("--integrate_region", "-i",
	help = "integrate within the selected boundaries",
	nargs = "*",
	metavar = "MIN-MAX"
)
integration$add_argument("--region_ratio", "-r",
	help = "compare region sums",
	nargs = "*",
	metavar = "MIN-MAX"
)
integration$add_argument("--dv200",
	help = "compute DV200",
	action = "store_true"
)
integration$add_argument("--illumina",
	help = "compute Illumina library ratio",
	action = "store_true"
)
integration$add_argument("--bound_variable",
	help = "variable of bounds",
	default = "length",
	metavar = "VARIABLE"
)
integration$add_argument("--sum_variable",
	help = "variable to integrate",
	default = "molarity",
	metavar = "VARIABLE"
)


qc <- parser$add_argument_group("QC", "Quality-control plots. See help(qc.electrophoresis).")
qc$add_argument("--stdcrv",
	help = "draw mobility standard curve(s)",
	action = "store_true"
)
qc$add_argument("--qc",
	nargs = "*",
	help = "draw QC plot for the selected variable(s)",
	metavar = "VARIABLE"
)

args <- parser$parse_args()

data <- do.call(read.electrophoresis, c(as.list(args$xml_files), method = args$method, mc.cores = args$mc_cores))

# annotations
if (! is.null(args$annotation_file)) {
	annot.args <- if (is.null(args$annot_args)) NULL else eval(parse(textConnection(paste0("list(", args$annot_args, ")")))) # put annotation arguments into a list expression, then parse and evaluate it
	data <- do.call(annotate.electrophoresis, c(list(data, args$annotation_file), annot.args)) 
}

# integration and table generation
result <- data$samples
if (! is.null(args$integrate_region)) {
	bounds.list <- lapply(strsplit(args$integrate_region, "-"), as.numeric)
	result <- do.call(data.frame, c(result, lapply(bounds.list, function(bounds.pair) integrate.custom(data, lower.bound = bounds.pair[1], upper.bound = bounds.pair[2], bound.variable = args$bound_variable, sum.variable = args$sum_variable)), check.names = F, stringsAsFactors = F))
	colnames(result)[(ncol(result) - length(bounds.list) + 1):ncol(result)] <- sapply(args$integrate_region, function(bounds) paste(args$sum_variable, "in", bounds))
}
if (! is.null(args$region_ratio)) {
	bounds.list <- lapply(strsplit(args$region_ratio, "-"), as.numeric)
	result <- data.frame(result, do.call(region.ratio, c(list(data), bounds.list, list(bound.variable = args$bound_variable, sum.variable = args$sum_variable))), check.names = F, stringsAsFactors = F)
}
if (args$dv200) result$DV200 <- dv200(data)
if (args$illumina) result$`Illumina library ratio` <- illumina.library.ratio(data)
write.table(result, file = if (is.null(args$output_file)) stdout() else args$output_file, quote = F, sep = "\t", row.names = F)

if (! is.null(args$plot_file)) {
	# parse dimensions
	dimensions <- as.numeric(strsplit(args$dimensions, "x")[[1]])
	stopifnot(length(dimensions) == 2)
	
	# parse faceting
	facets <- if (args$facets == "none") NULL else as.formula(args$facets)
	
	write(paste("writing", args$plot_file), stderr())
	pdf(args$plot_file, width = dimensions[1], height = dimensions[2])
		print(qplot.electrophoresis(
			data,
			x = args$x,
			y = args$y,
			log = args$log,
			normalize = args$normalize,
			facets = facets,
			margins = args$margins,
			scales = args$scales,
			geom = args$geom,
			include.ladder = args$include_ladder,
			include.markers = args$include_markers,
			lower.marker.spread = args$lower_marker_spread,
			xlim = if (! is.null(args$xlim)) args$xlim else c(NA, NA),
			ylim = if (! is.null(args$ylim)) args$ylim else c(NA, NA),
			show.peaks = ! args$hide_peaks,
			region.alpha = args$region_alpha,
			area.alpha = args$area_alpha,
			title = args$title,
			xlab = args$xlab,
			ylab = args$ylab
		))
		if (args$stdcrv) print(stdcrv.mobility(data))
		for (variable in args$qc) print(qc.electrophoresis(data, variable))
	dev.off()
	write("done", stderr())
}

