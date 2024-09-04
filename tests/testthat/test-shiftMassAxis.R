test_that("shiftMassAxis stops correctly", {
  data(Blank2022spec)
  
  # length(spec) != length(mzdiff)
  expect_error(shiftMassAxis(spec = Blank2022spec, mzdiff = 1))
})

test_that("shiftMassAxis works for single spectra/peaks", {
  data(Blank2022spec)
  data(Blank2022peaks)
  
  # spectra
  shifted <- shiftMassAxis(spec = Blank2022spec[[1]], mzdiff = 1)
  expect_true(length(MALDIquant::mass(shifted)) == length(MALDIquant::mass(Blank2022spec[[1]])))
  # round because it seems that 1 some times is not 1...
  expect_true(all(round(MALDIquant::mass(shifted) - MALDIquant::mass(Blank2022spec[[1]])) == 1))
  
  
  # peaks
  shifted <- shiftMassAxis(spec = Blank2022peaks[[1]], mzdiff = 1)
  expect_true(length(MALDIquant::mass(shifted)) == length(MALDIquant::mass(Blank2022peaks[[1]])))
  expect_true(all(round(MALDIquant::mass(shifted) - MALDIquant::mass(Blank2022peaks[[1]])) == 1))
})

test_that("shiftMassAxis works for spectra/peak-lists", {
  data(Blank2022spec)
  data(Blank2022peaks)
  
  # spectra
  len <- length(Blank2022spec)
  shifted <- shiftMassAxis(spec = Blank2022spec, mzdiff = rep(1, len))
  expect_true(length(shifted) == length(Blank2022spec))
  
  lengthCorrect <- vapply(1:len, 
                       function(i) {
                         length(MALDIquant::mass(shifted[[i]])) == length(MALDIquant::mass(Blank2022spec[[i]]))
                       }, 
                       FUN.VALUE = logical(1))
  expect_true(all(lengthCorrect))
  
  shiftCorrect <- vapply(1:len, 
                          function(i) {
                            all(round(MALDIquant::mass(shifted[[i]]) - MALDIquant::mass(Blank2022spec[[i]])) == 1)
                          }, 
                          FUN.VALUE = logical(1))
  expect_true(all(shiftCorrect))
  
  # peaks
  len <- length(Blank2022peaks)
  shifted <- shiftMassAxis(spec = Blank2022peaks, mzdiff = rep(1, len))
  expect_true(length(shifted) == length(Blank2022peaks))
  
  lengthCorrect <- vapply(1:len, 
                          function(i) {
                            length(MALDIquant::mass(shifted[[i]])) == length(MALDIquant::mass(Blank2022peaks[[i]]))
                          }, 
                          FUN.VALUE = logical(1))
  expect_true(all(lengthCorrect))
  
  shiftCorrect <- vapply(1:len, 
                         function(i) {
                           all(round(MALDIquant::mass(shifted[[i]]) - MALDIquant::mass(Blank2022peaks[[i]])) == 1)
                         }, 
                         FUN.VALUE = logical(1))
  expect_true(all(shiftCorrect))
})