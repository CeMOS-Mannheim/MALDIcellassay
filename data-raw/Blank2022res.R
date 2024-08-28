data("Blank2022spec")
Blank2022res <- fitCurve(spec = Blank2022spec,
                         SinglePointRecal = TRUE, 
                         normMz = 760.585, 
                         alignTol = 0.1, 
                         normTol = 0.1, 
                         varFilterMethod = "mean")

usethis::use_data(Blank2022res, overwrite = TRUE, compress = "xz")
