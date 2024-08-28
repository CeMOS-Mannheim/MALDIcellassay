library(MALDIquant)
data("Blank2022spec")
conc <- names(Blank2022spec)
Blank2022peaks <- MALDIquant::detectPeaks(Blank2022spec, 
                                 method = "SuperSmoother",
                                 SNR = 5)

names(Blank2022peaks) <- conc

usethis::use_data(Blank2022peaks, overwrite = TRUE)
