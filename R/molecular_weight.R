# molecular weight as a function of length, with relevant scaling factor, i.e. for DNA and RNA we're converting ng/uL to nmol/L (or pg/uL to pmol/L) so we need to scale Daltons by 1E6
# source: https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-tools-and-calculators/dna-and-rna-molecular-weights-and-conversions.html
#' @export
molecular.weight <- list(
	DNA = function(length) (length * 607.4 + 157.9)/1E6,
	RNA = function(length) (length * 320.5 + 159.0)/1E6
)
