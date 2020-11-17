Simple R functions for importing and analyzing electrophoresis data from an Agilent Bioanalyzer or TapeStation.

# Installation

Install required dependencies:

    > install.packages(c("XML", "base64enc", "png", "ggplot2"))

Easy way to get the newest release (includes 30 MB of demo data):

    > install.packages("https://github.com/jwfoley/bioanalyzeR/releases/download/v0.5.0/bioanalyzeR_0.5.0.tar.gz", repos = NULL)

or, easy way without the demo data:

    > install.packages("https://github.com/jwfoley/bioanalyzeR/releases/download/v0.5.0/bioanalyzeR_0.5.0-no_data.tar.gz", repos = NULL)

or, for hackers (it may take a minute to build the vignette):

    > library(devtools)
    > install_github("jwfoley/bioanalyzeR", build_vignettes = TRUE)

# Documentation

See the vignette [online](https://stanford.edu/~jwfoley/bioanalyzeR.html) or in R:

    > vignette("bioanalyzeR")

