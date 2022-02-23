#' generate ggplot objects for each of the curve fits in a MALDIassay object
#'
#' @param object    object of class MALDIassay
#' @param fc_thresh numeric, max/min fold-change above which plots should be generated
#' @param markValue numeric, a value to add as a reference to the plots
#' @param R2_thresh numeric, min. R-squared (goodness of curve fit) to plot curve
#'
#' @return
#' list of ggplot objects
#' @export
plotCurves <- function(object, fc_thresh = 1, R2_tresh = 0, markValue = NA) {
  stopIfNotIsMALDIassay(object)
  res_list <- getCurveFits(object)
  len <- length(names(res_list))
  mz_vals <- as.numeric(names(res_list))
  
  p_list <- vector("list", length = len)
  for(i in 1:len) {
    mz <- mz_vals[i]
    model <- res_list[[as.character(mz)]]$model
    df <- res_list[[as.character(mz)]]$df
    
    ic50 <- 10^getEstimates(model, targets = 0.5)[,3]
    min <- min(df$value)
    max <- max(df$value)
    fc_window <- max/min
    R2 <- getGoodness(model)[[1]]
    
    if(fc_window >= fc_thresh & R2 >= R2_tresh) {
      df_C <- tibble(xC = getXcurve(model), yC = getYcurve(model))
      df_P <- tibble(x = getX(model), y = getY(model))
      
      
      ggplot(data = df_P, aes(x = x, y = y)) +
        geom_line(data = df_C, aes(x = xC, y = yC)) +
        geom_point() +
        scale_x_continuous(labels = c(0, 10^df_P$x[-1]), breaks = df_P$x) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Conc.",
             y = "relative Int. [% of max Int.]",
             title = paste0("mz ", round(mz,2), " Da, R\u00B2=", round(R2,3), "\n",
                            "logIC50=", round(log(ic50),3), " min=", round(min, 3), " max=", round(max, 3), " FC=", round(fc_window, 2))) -> p
      
      
      if(!is.na(markValue)) {
        p <- p + geom_vline(aes(xintercept = markValue), linetype = "dashed")
      }
      p_list[[i]] <- p
      names(p_list) <- mz_vals
    }
  }
  # check for empty entries (result of filtering for FC or R2) and remove them
  idx <- vapply(p_list, function(x) {length(x)>0}, FUN.VALUE = TRUE)
  if(sum(!idx) == len) {
    stop("Nothing to plot. Condsider decreasing fc_thresh or R2_tresh.\n")
  }
  return(p_list[idx])
}
