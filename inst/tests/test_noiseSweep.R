context("Noise Sweep")

test_that("noise sweep works properly", {
  noiseGenes <- as.character(seq(1, 5000))
  nSamples <- 2
  sizeNoise <- seq(0, 500, 20)
  fracShared <- seq(0, 1, .1)
  
  expLenNoise <- length(sizeNoise)
  expLenShared <- length(fracShared)
  expNumNoise <- sizeNoise
  expNumShared <- round(rep(sizeNoise, each=expLenShared) * rep(fracShared, expLenNoise) )
  
  outNoise <- sweepNoiseSample(noiseGenes, nSamples, sizeNoise, fracShared)
  
  hasLenNoise <- length(outNoise)
  expect_equal(hasLenNoise, expLenNoise)
  
  hasNumNoise <- sapply(outNoise, function(x){length(x[[1]][[1]])})
  expect_equal(hasNumNoise, expNumNoise)
  
  hasNumShared <- lapply(outNoise, function(x){
    lapply(x, function(y){
      length(intersect(y[[1]], y[[2]]))
    })
  })
  
  expect_equal(unlist(hasNumShared), expNumShared)
})