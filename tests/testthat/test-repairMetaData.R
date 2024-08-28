# helper 
spectraListConstructor <- function(n = 3, metaData = list(test = "test")) {
  l <- lapply(1:n, 
              function(x) 
              {MALDIquant::createMassSpectrum(mass = 1:10, 
                                              intensity = abs(rnorm(10, 0, 1)), 
                                              metaData = metaData)
              })
  return(l)
} 

test_that(".repairMetaData removes NA entries", {
  s <- spectraListConstructor(metaData = list(test = NA))
  s <- suppressWarnings(.repairMetaData(s))
  expect_true(
    all(
      unlist(
        lapply(s, 
               function(x) {
                 MALDIquant::metaData(x)$test == ""
               })
      )
    )
  ) 
})

test_that(".repairMetaData keeps non NA entries", {
  s <- spectraListConstructor(metaData = list(test = "test"))
  s <- suppressWarnings(.repairMetaData(s))
  expect_true(
    all(
      unlist(
        lapply(s, 
               function(x) {
                 MALDIquant::metaData(x)$test == "test"
               })
      )
    )
  ) 
  
  s <- spectraListConstructor(metaData = list(test = 1))
  s <- suppressWarnings(.repairMetaData(s))
  expect_true(
    all(
      unlist(
        lapply(s, 
               function(x) {
                 MALDIquant::metaData(x)$test == 1
               })
      )
    )
  ) 
  
  s <- spectraListConstructor(metaData = list(test = list(x = 3, y = 2)))
  s <- suppressWarnings(.repairMetaData(s))
  expect_true(
    all(
      unlist(
        lapply(s, 
               function(x) {
                 md <- MALDIquant::metaData(x)
                 md$test$x == 3 & md$test$y == 2
                 
               })
      )
    )
  ) 
})