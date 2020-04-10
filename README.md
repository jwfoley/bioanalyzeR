Simple R functions for importing and analyzing electrophoresis data from an Agilent Bioanalyzer or TapeStation.

# Installation

Easy way to get the newest release (includes 30 MB of demo data):

    > install.packages("https://github.com/jwfoley/bioanalyzeR/releases/download/v0.3.0/bioanalyzeR_0.3.0.tar.gz")

Easy way without the demo data:

    > install.packages("https://github.com/jwfoley/bioanalyzeR/releases/download/v0.3.0/bioanalyzeR_0.3.0-no_data.tar.gz")

For hackers (it may take a minute to build the vignette):

    > library(devtools)
    > install_github("jwfoley/bioanalyzeR", build_vignettes = TRUE)

# Documentation

See the vignette [online](https://stanford.edu/~jwfoley/bioanalyzeR.html) or in R:

    > vignette("bioanalyzeR")

