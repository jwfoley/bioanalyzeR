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

