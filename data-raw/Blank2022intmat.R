library(MALDIquant)
data("Blank2022spec")
data("Blank2022peaks")

Blank2022intmat <- intensityMatrix(Blank2022peaks, Blank2022spec)
rownames(Blank2022intmat) <- names(Blank2022spec)

usethis::use_data(Blank2022intmat, overwrite = TRUE, compress = "xz")
