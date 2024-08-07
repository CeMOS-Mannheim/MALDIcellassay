#' Get mass shift for target mz - internal
#'
#' @param peaksdf       data.frame with peaks information as generated by peaks2df()
#' @param tol           Numeric, tolerance around targetMz
#' @param targetMz      Numeric, target mass
#' @param tolppm        Logical, is the tolerance provided in ppm (TRUE) or Dalton (FALSE)
#' @param allowNoMatch  Logical, stop if targetMz is not found in single spectrum?
#'                      If TRUE spectra without targetMz match will be excluded.
#'
#' @return             List with two entries:
#'                                       MzShift The mass shift for each spectrum
#'                                       specIdx The index of the spectra with a match for targetMz
#'
#' @importFrom dplyr %>% mutate filter group_by arrange pull .data

.getMzShift <- function(peaksdf, tol, targetMz, tolppm = TRUE, allowNoMatch = TRUE) {
  plot_Idx <- sort(unique(peaksdf$plotIdx))

  if(tolppm) {
    tol <- tol / 1e6
  }


  f_resdf <- peaksdf %>%
    mutate(match = .data$mz > targetMz - tol & .data$mz < targetMz + tol) %>%
    filter(match) %>%
    mutate(mz.diff = round(targetMz - .data$mz, 4)) %>%
    group_by(.data$plotIdx) %>%
    filter(abs(.data$mz.diff) == min(abs(.data$mz.diff))) %>%
    arrange(.data$plotIdx)



  if (!all(plot_Idx %in% (f_resdf %>% pull(.data$plotIdx)))) {
    if (!allowNoMatch) {
      stop("Could not find ", targetMz, " for all spectra! Consider adjusting tol.\n")
    }
    warning("Could not find ", targetMz, " in spectrum ", paste(which(!(plot_Idx %in% (f_resdf %>% pull(.data$plotIdx)))), collapse = ", "), ".\n")
    specIdx <- sort(which(plot_Idx %in% (f_resdf %>% pull(.data$plotIdx))))
  } else {
    specIdx <- plot_Idx
  }
  if (length(specIdx) < 1) {
    stop("Could not find targetMz in any spectrum! Consider adjusting tol.\n")
  }

  return(list(
    mzshift = pull(f_resdf, .data$mz.diff),
    specIdx = specIdx
  ))
}


#' Get mass shift for target mz
#'
#' @param peaks    List of MALDIquant::MassPeak
#' @param targetMz Numeric, target mass
#' @param tol      Numeric, tolerance around targetMz
#' @param tolppm   Logical, tolerance supplied in ppm
#'
#' @return
#' List with two entries:
#' `MzShift` The mass shift for each spectrum
#' `specIdx` The index of the spectra with a match for targetMz
#' @export
getMzShift <- function(peaks,
                       targetMz,
                       tol,
                       tolppm = FALSE) {
  stopifnot(isMassPeaksList(peaks))
  nm <- names(peaks)
  stopifnot(!is.null(nm))
  stopifnot(is.numeric(as.numeric(nm)))


  # perform single point mass recalibration
  mzShift <- .getMzShift(
    peaksdf = peaks2df(peaks),
    tol = tol,
    targetMz = targetMz,
    tolppm = tolppm,
    allowNoMatch = TRUE
  )
  cat("found mz", targetMz, "in", length(mzShift$specIdx), "/",
      length(peaks), "spectra\n")
  cat(timeNow(), "mzshift was", mean(mzShift$mzshift),
      "in mean and", max(abs(mzShift$mzshift)), " abs. max.\n")

  if(length(unique(nm)) != length(unique(nm[mzShift$specIdx]))) {
    # stop if a single condition got filtered completely
    u_nm <- unique(nm)
    u_fil <- unique(nm[mzShift$specIdx])
    label_removed <- u_nm[which(!(u_nm %in% u_fil))]

    stop("Could not find ", targetMz, " in any spectrum with label ",
         paste0(label_removed, collapse = ", "),
         ".\n Consider increasing tol.\n")
  }
  return(mzShift)
}

