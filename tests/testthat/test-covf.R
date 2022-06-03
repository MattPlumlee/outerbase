testcovgrad <- function(covname) {
  
  ss = 10
  covobj = new(get(paste("covf_",covname,sep="")))
  range = covobj$uppbnd-covobj$lowbnd
  xo = seq(covobj$lowbnd+range/ss/2,covobj$uppbnd-range/ss/2,
           length.out=ss)
  
  hyp0 = covobj$hyp0 + (covobj$hypub-covobj$hyplb)*
    (runif(length(covobj$hyp))-0.5)/2
  
  covobj$hyp = hyp0
  A = covobj$cov(xo, xo)
  A_gh = covobj$cov_gradhyp(xo, xo)
  
  eps = 10^(-4) 
  totdiff = rep(0,length(covobj$hyp))
  for(k in 1:length(covobj$hyp)){
    hyph = hyp0
    hyph[k] = hyp0[k]+eps
    covobj$hyp = hyph
    
    Aalt = covobj$cov(xo, xo)
    totdiff[k] = sum(abs((Aalt-A)/eps-A_gh[,,k]))
  }
  sum(abs(totdiff))
}


testdiagcov <- function(covname) {
  
  ss = 10
  covobj = new(get(paste("covf_",covname,sep="")))
  
  set.seed(42)
  range = covobj$uppbnd-covobj$lowbnd
  xo = seq(covobj$lowbnd+range/ss/2,covobj$uppbnd-range/ss/2,
           length.out=ss)
  
  A = covobj$cov(xo, xo)
  dA = covobj$covdiag(xo)
  
  sum(abs(dA-diag(A)))
}


test_that("mat25", {
  expect_equal(testdiagcov('mat25'), 0, tolerance = 0.01,
               label="grad ok")
  expect_equal(testcovgrad('mat25'), 0, tolerance = 0.01,
               label="grad ok")
})

test_that("mat25ang", {
  expect_equal(testdiagcov('mat25ang'), 0, tolerance = 0.01,
               label="grad ok")
  expect_equal(testcovgrad('mat25ang'), 0, tolerance = 0.01,
               label="grad ok")
})

test_that("mat25pow", {
  expect_equal(testdiagcov('mat25pow'), 0, tolerance = 0.01,
                label="grad ok")
  expect_equal(testcovgrad('mat25pow'), 0, tolerance = 0.01,
               label="grad ok")
})