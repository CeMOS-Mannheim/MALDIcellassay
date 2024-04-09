#' load bruker MALDI target plate spectra
#'
#' @param Dir         Character, parent directory of spectra.
#' @param filter      Character vector, filter out spectra which match the given vector.
#' @param nameSpectra Logical, if TRUE the spectra in the resulting list will be named according to the dirname.
#'
#' @return List of MALDIquant::MassSpectra
#' @export
#'
#' @importFrom svMisc progress
#' @importFrom MALDIquantForeign importBrukerFlex

loadSpectra <- function(Dir, filter = NA, nameSpectra = TRUE) {
  # get names of all anaylses
  analyses <- basename(list.dirs(Dir, recursive = F))

  # filter Calibration spectra
  analyses <- analyses[which(!analyses %in% filter)]



  # determine number of measured spots
  total_n <- vector("numeric", 1)
  for (i in 1:length(analyses)) {
    n <- length(list.dirs(file.path(Dir, analyses[i]), recursive = F))
    total_n <- total_n + n
  }

  spectra <- vector("list", length = total_n)
  counter <- 0
  cat(timeNow(), "Loading spectra...\n\n")
  for (i in analyses) {
    path <- file.path(Dir, i)
    spots <- list.files(path, recursive = T)[grepl(pattern = "fid", x = list.files(path, recursive = T))]
    for (spot in spots) {
      counter <- counter + 1
      spec <- MALDIquantForeign::importBrukerFlex(path = file.path(path, spot), verbose = F)
      # metaData(spec[[1]]) <- list(name = i)
      spectra[[counter]] <- spec[[1]]

      if (nameSpectra) {
        names(spectra)[counter] <- i
      }

      svMisc::progress(counter / total_n * 100)
    }
  }
  cat("\n")

  return(spectra)
}

#' load mzML spectra
#'
#' @param Dir         Character, parent directory of spectra.
#' @param filter      Character vector, filter out spectra which match the given vector.
#' @param nameSpectra Logical, if TRUE the spectra in the resulting list will be named according to the dirname.
#'
#' @return List of MALDIquant::MassSpectra
#' @export
#'
#' @importFrom svMisc progress
#' @importFrom MALDIquantForeign importMzMl
loadSpectraMzML <- function(Dir, filter = NA, nameSpectra = TRUE) {
  # get names of all anaylses
  analyses <- list.files(Dir, recursive = F)

  # filter Calibration spectra
  analyses <- analyses[which(!analyses %in% filter)]

  spectra <- vector("list", length = length(analyses))

  cat(MALDIcellassay:::timeNow(), "Loading spectra...\n\n")
  counter <- 0
  for (i in analyses) {
    counter <- counter + 1
    path <- file.path(Dir, i)

    spec <- MALDIquantForeign::importMzMl(path = file.path(path), verbose = F)


    if (nameSpectra) {
      names(spec) <- rep(tools::file_path_sans_ext(i), length(spec))
    }

    spectra[[counter]] <- spec

    svMisc::progress(counter / length(analyses) * 100)
  }

  cat("\n")

  return(unlist(spectra))
}
