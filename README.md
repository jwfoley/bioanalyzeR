Simple R functions for importing and analyzing electrophoresis data from Agilent systems (Bioanalyzer, TapeStation, Fragment Analyzer, ZAG DNA Analyzer, Femto Pulse).

# Installation

Install required dependencies:

    > install.packages(c("XML", "base64enc", "png", "plyr", "ggplot2"))

Install the newest release (includes 22 MB of demo data):

    > install.packages("https://github.com/jwfoley/bioanalyzeR/releases/download/v0.8.1/bioanalyzeR_0.8.1.tar.gz", repos = NULL)

or, install the newest release without the demo data:

    > install.packages("https://github.com/jwfoley/bioanalyzeR/releases/download/v0.8.1/bioanalyzeR_0.8.1-no_data.tar.gz", repos = NULL)


# Documentation

See the vignette [online](https://stanford.edu/~jwfoley/bioanalyzeR.html) or in R:

    > vignette("bioanalyzeR")

