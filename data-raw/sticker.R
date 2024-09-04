library(MALDIcellassay)
library(MALDIquant)
library(ggplot2)
library(hexSticker)
library(tidyverse)

set.seed(42)

# helper
plotCurves2 <- function(object,
                        mzIdx = NULL,
                        errorbars = c("none", "sd", "sem"),
                        col = "black",
                        linetype = "solid",
                        shape = 21,
                        limits = NULL) {
  MALDIcellassay:::stopIfNotIsMALDIassay(object)
  errorbars <- match.arg(errorbars)
  
  if (!is.null(mzIdx)) {
    res_list <- getCurveFits(object)[mzIdx]
  }
  len <- length(names(res_list))
  mz_vals <- as.numeric(names(res_list))
  
  p_list <- vector("list", length = len)
  for (i in 1:len) {
    mz <- mz_vals[i]
    model <- res_list[[as.character(mz)]]$model
    df <- res_list[[as.character(mz)]]$df
    
    ic50 <- MALDIcellassay:::.getEstimates(model, 0.5)
    min <- min(df$value)
    max <- max(df$value)
    
    fc <- MALDIcellassay:::calculateFC(model)
    R2 <- nplr::getGoodness(model)[[1]]
    
    df_C <- tibble(xC = nplr::getXcurve(model),
                   yC = nplr::getYcurve(model))
    df_P <- tibble(x = nplr::getX(model),
                   y = nplr::getY(model))
    
    int <- vapply(getSinglePeaks(object), function(x) {
      targetmass <- mz
      mass <- mass(x)
      idx <- match.closest(targetmass, mass, tolerance = 0.01)
      int <- intensity(x)
      return(int[idx])
    }, numeric(1))
    
    df_singlePeaks <- tibble(
      x = getConc(object),
      int = int) %>%
      group_by(x) %>%
      summarise(
        sd = sd(int),
        sem = sd/sqrt(n())
      )
    
    df_P <- df_P %>%
      mutate(sem = pull(df_singlePeaks, sem),
             sd = pull(df_singlePeaks, sd))
    
    p <- ggplot(data = df_P, aes(x = -x, y = y)) +
      geom_line(data = df_C,
                aes(x = -xC, y = yC),
                linewidth = 0.5,
                alpha = 1,
                col = col,
                linetype = linetype) +
      geom_point(size = 1,
                 shape = shape,
                 col = col,
                 stroke = 1) +
      scale_x_continuous(labels = scales::comma,
                         limits = limits,
                         n.breaks = 5,
                         trans = "reverse")
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        x = "Conc.",
        y = "Intensity",
        title = paste0(
          "mz ", round(mz, 2), " Da, R\u00B2=", round(R2, 3), "\n",
          "pIC50=", round(ic50, 3),
          " min=", round(min, 3),
          " max=", round(max, 3),
          " FC=", round(fc, 2)
        )
      )
    
    if (errorbars == "sem") {
      p <- p +
        geom_errorbar(
          aes(
            y = y,
            ymin = y - sem,
            ymax = y + sem
          ),
          alpha = 1,
          width = 0.2,
          size = 1,
          col = col
        )
    }
    
    if (errorbars == "sd") {
      p <- p +
        geom_errorbar(
          aes(
            y = y,
            ymin = y - sd,
            ymax = y + sd
          ),
          alpha = 1,
          width = 0.2,
          size = 1,
          col = col
        )
    }
    
    p_list[[i]] <- p
    names(p_list) <- mz_vals
    
  }
  # check for empty entries (result of filtering for FC or R2) and remove them
  idx <- vapply(p_list, function(x) {
    length(x) > 0
  }, FUN.VALUE = TRUE)
  if (sum(!idx) == len) {
    stop("Nothing to plot. Condsider decreasing fc_thresh or R2_tresh.\n")
  }
  
  if (length(p_list) == 1) {
    return(p_list[[1]])
  }
  
  return(p_list[idx])
}

# download data and unzip
curl::curl_download("https://figshare.com/ndownloader/files/46156788", "testdata_mzML.zip")
if(!file.exists("testdata_mzML.zip")) {
  stop("Could not download testdata_mzML.zip.\n")
}

unzip("testdata_mzML.zip")

fs::file_delete("testdata_mzML.zip")

# overwrite settings
if(!file.exists("data-raw/settings_mzML_data.csv")) {
  stop("Could not find settings_mzML_data.csv.\n")
}
fs::file_copy(path = "data-raw/settings_mzML_data.csv",
              new_path = "settings.csv",
              overwrite = TRUE)

spec <- loadSpectraMzML("mzML")
settings <- read.csv("settings.csv")
fs::file_delete("settings.csv")
fs::dir_delete("mzML")

spec_prc <- MALDIquant::removeBaseline(spec)
conc <- as.numeric(names(spec))
names(spec_prc) <- conc

res <- fitCurve(spec_prc, 
                unit = settings$concUnits, 
                varFilterMethod = settings$VarFilterMethod, 
                monoisotopicFilter = settings$monoisotopicFilter, 
                normMz = settings$normMz, 
                normTol = settings$normTol, 
                averageMethod = settings$avgMethod, 
                binTol = settings$binTol, 
                alignTol = settings$alignTol, 
                SNR = settings$SNR, 
                SinglePointRecal = settings$SinglePointRecal, 
                normMeth = settings$normMeth, 
                halfWindowSize = settings$halfWindowSize)

idx <- MALDIquant::match.closest(table = getAllMz(res), 
                                 x = 349.11, 
                                 tolerance = 0.1)

val <- map(unique(conc), function(x) {abs(rnorm(mean = x, n = 27, sd = x/2))}) %>% unlist() %>% .[-c(1:3)]
well <- paste0(LETTERS[2:13], rep(2:21, each = 12))

pm <- 
  platetools::raw_map(data = log10(val), well = well, plate = "384", shape = 21, size = 1.5, na_alpha = 0) + 
  theme_void() + 
  theme_transparent() +
  theme(legend.position = "none") +
  scale_fill_viridis_c(end = 0.85)

p <- 
  plotCurves2(object = res, 
           mzIdx = idx, 
           errorbars = "none") + 
  theme_void() + 
  theme_transparent() +
  theme(text = element_blank()) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  scale_x_continuous(expand = expansion(mult = 0.1))


pall <- pm + annotation_custom(grob = ggplotGrob(p))
pall

sticker(package = "MALDIcellassay", 
        subplot = pall, 
        filename = "man/figures/MALDIcellassay_sticker.png", 
        dpi = 300, 
        p_size=16, p_fontface = "bold", 
        s_x=1.05, 
        s_y=0.8, 
        s_width=1.7, 
        s_height=1.6, 
        h_fill = "white", 
        h_color = "midnightblue", 
        p_color = "black")

usethis::use_logo("man/figures/MALDIcellassay_sticker.png")