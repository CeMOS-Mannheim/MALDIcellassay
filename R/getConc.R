#' Extract the concentrations used in a MALDIassay
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' Concentrations used in a MALDIassay
#' @export

getConc <- function(object) {
  if(!class(object) == "MALDIassay") {
    stop("object needs to be of class MALDIassay.")
  }
  return(object@settings$Conc)
}


