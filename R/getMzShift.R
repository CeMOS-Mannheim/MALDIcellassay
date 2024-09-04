#' Get mass shift for target mz
#'
#' @param peaks    List of MALDIquant::MassPeak
#' @param targetMz Numeric, target mass
#' @param tol      Numeric, tolerance around targetMz
#' @param tolppm   Logical, tolerance supplied in ppm
#' @param verbose  Logical, print logs to the console.
#'
#' @return
#' List with two entries:
#' `MzShift` The mass shift for each spectrum
#' `specIdx` The index of the spectra with a match for targetMz
#' @export
#' @importFrom MALDIquant match.closest
#' @importFrom purrr map_vec
#' 
#' @examples
#' data(Blank2022peaks)
#' getMzShift(Blank2022peaks, targetMz = 760.585, tol = 0.1, tolppm = FALSE)
getMzShift <- function(peaks,
                       targetMz,
                       tol,
                       tolppm = FALSE,
                       verbose = TRUE) {
  stopifnot(isMassPeaksList(peaks))
  nm <- names(peaks)
  stopifnot(!is.null(nm))
  stopifnot(all(!is.na((as.numeric(nm)))))
  
  allIdx <- 1:length(peaks)
  
  if(tolppm) {
    tol <- tol * targetMz / 1e6
  }
  
  mzShift <- map_vec(peaks, 
          function(x) {
            mz <- mass(x)
            
            idx <- match.closest(targetMz, 
                          table = mz, 
                          tolerance = tol, 
                          nomatch = NA_integer_)
            mzdiff <- targetMz-mz[idx]
            return(mzdiff)
          })
  
  specIdx <- which(!is.na(mzShift))
  mzShift <- na.omit(mzShift)
  
  if (!all(allIdx %in% specIdx)) {
    if (!allowNoMatch) {
      stop("Could not find ", targetMz, " for all spectra! Consider adjusting tol.\n")
    }
    warning("Could not find ", targetMz, " in spectrum ", paste(which(!(allIdx %in% specIdx)), collapse = ", "), ".\n")
    specIdx <- sort(which(allIdx %in% specIdx))
  } else {
    specIdx <- allIdx
  }
  if (length(specIdx) < 1) {
    stop("Could not find targetMz in any spectrum! Consider adjusting tol.\n")
  }
  
  if(verbose) {
    cat("found mz", targetMz, "in", length(specIdx), "/",
        length(peaks), "spectra\n")
    cat(timeNow(), "mzshift was", mean(mzShift),
        "in mean and", max(abs(mzShift)), " abs. max.\n")  
  }

  if(length(unique(nm)) != length(unique(nm[specIdx]))) {
    # stop if a single condition got filtered completely
    u_nm <- unique(nm)
    u_fil <- unique(nm[specIdx])
    label_removed <- u_nm[which(!(u_nm %in% u_fil))]

    stop("Could not find ", targetMz, " in any spectrum with label ",
         paste0(label_removed, collapse = ", "),
         ".\n Consider increasing tol.\n")
  }
  return(list(mzshift = mzShift,
              specIdx = specIdx))
}

