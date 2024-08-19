library(MALDIcellassay)
dataDir <- system.file("extdata", package="MALDIcellassay")
unzip(file.path(dataDir, "example-raw-spectra.zip"))
spec <- loadSpectra("example-raw-spectra/")
unlink("example-raw-spectra/", recursive = TRUE)

MALDIquantForeign::exportMzMl(spec, file = "inst/extdata/0.0001.mzML", force = TRUE)