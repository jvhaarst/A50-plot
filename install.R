# Script to install dependancies.

# Set CRAN mirror to download from
options("repos" = c(CRAN = "http://cran.r-mirror.de/"))

# Install packagess
install.packages("gdata")
install.packages("plyr")
install.packages("sitools")

require("gdata")
require("plyr")
require("sitools")

source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("IRanges")

require("Biostrings")
require("IRanges")


