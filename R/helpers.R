#### helper functions ######

#' Current time
#'

#' @return current time
#' @noRd
timeNow <- function() {
  format(Sys.time(), "%H:%M")
}

#' get direction of curve
#'
#' @param model nplr model
#'
#' @return Numeric direction (positive or negative)
#'
#' @noRd
getDirection <- function(model) {
  y <- getYcurve(model)

  y0 <- y[1]
  y1 <- y[length(y)]

  return(y1-y0)
}

#' go up x folders
#'
#' @param path Character, Path
#' @param x    Numeric
#'
#' @return     new path
#' @noRd
goUpXFolders <- function(path, x) {
  for(i in 1:x) {
    resPath <- dirname(path)
    path <- resPath
  }
  return(resPath)
}

#' Convert spectrum to data.frame
#'
#' @param specs MALDIquant::MassSpectrum or list there of
#'
#' @return      data.frame of spectra
#' @noRd
spec2df <- function(specs) {
  df_l <- lapply(1:length(specs), function(i) {
    mz <- mass(specs[[i]])
    int <- intensity(specs[[i]])
    name <- names(specs)[i]
    return(tibble(name = name,
                  mz = mz,
                  int = int,
                  plotIdx = i))
  })
  return(bind_rows(df_l))
}
