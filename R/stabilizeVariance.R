stabilizeVariance <- function(spec,
                              method = c("sqrt", "log", "log2", "log10"),
                              correctBaseline = c(NA, "SNIP", "TopHat", "ConvexHull", "median")) {
  method <- match.arg(method)
  correctBaseline <- match.arg(correctBaseline)

  spec <- transformIntensity(spec, method = method)
  if(!is.na(correctBaseline)) {
    spec <- removeBaseline(spec, method = correctBaseline)
  }
  return(spec)
}
