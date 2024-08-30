#' Current time
#'
#' @return current time
#' @noRd
timeNow <- function() {
  format(Sys.time(), "%H:%M")
}

#' Convert concentration to log10 and replace zero's
#'
#' @param conc numeric, concentrations.
#'
#' @return
#' numeric, log10 transformed concentrations
#' @export
#' @examples
#' transformConc2Log(c(0.1, 0.01,0.001))
transformConc2Log <- function(conc) {
  concLog <- log10(conc)
  if (any(concLog == -Inf)) {
    concLog[which(concLog == -Inf)] <- (min(concLog[which(!concLog == -Inf)]) - 1)
  }
  return(concLog)
}


#' Check if object is of class MALDIassay
#'
#' @param object object to text
#'
#' @return
#' logical, TRUE if object is of class MALDIassay
#' @export
#' @examples 
#' x <- 1
#' # FALSE
#' isMALDIassay(x)
#' # TRUE
#' isMALDIassay(Blank2022res)
isMALDIassay <- function(object) {
  if (!is(object, "MALDIassay")) {
    return(FALSE)
  }
  return(TRUE)
}

#' Stop and throw an error if object is not of class MALDIassay
#'
#' @param object object
#' @noRd
stopIfNotIsMALDIassay <- function(object) {
  if (!isMALDIassay(object)) {
    stop("object needs to be of class MALDIassay.")
  }
}

#' Calculate the fold-change
#'
#' @param model nplr model object
#'
#' @return
#' Numeric, fold-change
#'
#' @importFrom tibble is_tibble
#' @noRd
calculateFC <- function(model) {
  par <- getPar(model)$params

  if(par[["top"]] <= 0 | par[["bottom"]] <= 0) {
    return(NA)
  }
  
  # Define FC as largeValue/smallValue
  # Is this correct?
  if(par[["scal"]] < 0) {
    FC <- par[["bottom"]]/par[["top"]]
    return(FC)
  }

  FC <- par[["top"]]/par[["bottom"]]

  return(FC)
}
