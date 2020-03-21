#! /usr/bin/Rscript

library(bioanalyzeR)
library(argparse)

parser <- ArgumentParser(description = "Simple automation of bioanalyzeR functions.")

files <- parser$add_argument_group("files", "Filenames and settings for input/output.")
files$add_argument("xml_files",
	nargs = "+",
	help = "one or more XML files exported from the Bioanalyzer or TapeStation software",
	metavar = "INPUT_FILE.xml"
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
files$add_argument("--fit",
	help = "mobility standard curve method",
	default = "spline",
	choices = c("spline", "interpolation", "regression")
)

plotting <- parser$add_argument_group("plotting", "Settings for the electropherogram plot. See help(qplot.electrophoresis) for more information.")
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
	choices = c("", "x", "y", "xy")
)
plotting$add_argument("--facets", "-f",
	help = "faceting formula, or 'none'",
	default = "~ sample.index",
	metavar = "FORMULA"
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
plotting$add_argument("--peak_fill",
	help = "color of peaks",
	default = "darkred",
	metavar = "COLOR"
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

integration <- parser$add_argument_group("integration", "Integrate a variable under the curve. See help(integrate.custom) and help(region.ratio) for more information.")
integration$add_argument("--integrate_region", "-i",
	help = "integrate within the selected boundaries",
	nargs = "*",
	metavar = c("MIN-MAX")
)
integration$add_argument("--region_ratio", "-r",
	help = "compare region sums",
	nargs = "*",
	metavar = c("MIN-MAX")
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


qc <- parser$add_argument_group("QC", "Quality-control plots.")
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

data <- do.call(read.electrophoresis, c(as.list(args$xml_files), fit = args$fit))

# integration and table generation
result <- data$samples
if (! is.null(args$integrate_region)) {
	bounds.list <- lapply(strsplit(args$integrate_region, "-"), as.numeric)
	result <- do.call(data.frame, c(result, lapply(bounds.list, function(bounds.pair) integrate.custom(data, lower.bound = bounds.pair[1], upper.bound = bounds.pair[2], bound.variable = args$bound_variable, sum.variable = args$sum_variable)), check.names = F, stringsAsFactors = F))
}
if (! is.null(args$region_ratio)) {
	bounds.list <- lapply(strsplit(args$region_ratio, "-"), as.numeric)
	result <- data.frame(result, region.ratio(data, bounds = bounds.list, bound.variable = args$bound_variable, sum.variable = args$sum_variable), check.names = F, stringsAsFactors = F)
}
if (args$illumina) result <- data.frame(result, illumina.library.ratio(data), check.names = F, stringsAsFactors = F)
write.table(result, file = if (is.null(args$output_file)) stdout() else args$output_file, quote = F, sep = "\t", row.names = F)

if (! is.null(args$plot_file)) {
	# parse dimensions
	dimensions <- as.numeric(strsplit(args$dimensions, "x")[[1]])
	stopifnot(length(dimensions) == 2)
	
	# parse faceting
	facets <- if (args$facets == "none") NULL else as.formula(args$facets)
	
	write(paste("writing", args$plot_file), stderr())
	pdf(args$plot_file, width = dimensions[1], height = dimensions[2])
		print(qplot.electrophoresis(data,
			x = args$x,
			y = args$y,
			log = args$log,
			facets = facets,
			scales = args$scales,
			geom = args$geom,
			include.ladder = args$include_ladder,
			between.markers = ! args$include_markers,
			peak.fill = args$peak_fill,
			region.alpha = args$region_alpha,
			area.alpha = args$area_alpha
		))
		if (args$stdcrv) print(stdcrv.mobility(data))
		for (variable in args$qc) print(qc.electrophoresis(data, variable))
	dev.off()
	write("done", stderr())
}

