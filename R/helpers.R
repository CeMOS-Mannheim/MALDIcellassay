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

#' Check if object if of class MALDIassay
#'
#' @param object object
#'
#' @return
#' logical
#' @export
isMALDIassay <- function(object) {
  if(!class(object) == "MALDIassay") {
    return(FALSE)
  }
  return(TRUE)
}

#' Stop and throw an error if object is not of class MALDIassay
#'
#' @param object object
#'
#' @export
stopIfNotIsMALDIassay <- function(object) {
  if(!isMALDIassay(object)) {
    stop("object needs to be of class MALDIassay.")
  }
}

#' Write ggplot objects to disk as png
#'
#' @param object object of class MALDIassay
#' 
#' @importFrom ggplot2 ggsave
#'
savePlots <- function(object, dir = NULL) {
  p_list <- plotCurves(object)
  normMeth <- getNormMethod(object)
  if(is.null(dir)) {
    # if no directory is given, take it from the MALDIassay object.
    dir <- getDirectory(object)
  }
  
  for(i in 1:length(p_list)) {
    p <- p_list[[i]]
    mz <- as.numeric(names(p_list[i]))
    ggsave(filename = file.path(dir, paste0(as.character(Sys.Date()),"_plotR2_", 
                                            normMeth, "norm_", round(mz,2),".png")), 
           plot = p)
  }
}
