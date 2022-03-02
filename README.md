Simple R functions for importing and analyzing electrophoresis data from Agilent systems (Bioanalyzer, TapeStation, Fragment Analyzer, ZAG DNA Analyzer, Femto Pulse).

# Installation

Install the required dependencies:

    > install.packages(c("XML", "base64enc", "png", "plyr", "ggplot2"))

Install the newest release (includes 22 MB of demo data):

    > install.packages("https://github.com/jwfoley/bioanalyzeR/releases/download/v0.9.1/bioanalyzeR_0.9.1.tar.gz", repos = NULL)

or, install the newest release without the demo data:

    > install.packages("https://github.com/jwfoley/bioanalyzeR/releases/download/v0.9.1/bioanalyzeR_0.9.1-no_data.tar.gz", repos = NULL)


# Documentation

See the vignette [online](https://stanford.edu/~jwfoley/bioanalyzeR.html) or in R:

    > vignette("bioanalyzeR")


# Support

This package is a work in progress and likely to change significantly from version to version. Watch [the repository on GitHub](https://github.com/jwfoley/bioanalyzeR) and be sure to keep up to date with the latest release. **R will not automatically install updates** because the package is not yet hosted on an official R repository.

If you discover a problem with the package (even if it's just confusing documentation), first verify you have the [latest release](https://github.com/jwfoley/bioanalyzeR/releases) installed. Then the best way to get help is by submitting a [GitHub issue](https://github.com/jwfoley/bioanalyzeR/issues). Otherwise email the author. Please include minimal code and data required to reproduce the problem; you can anonymize the samples if appropriate.

