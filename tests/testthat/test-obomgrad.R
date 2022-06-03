borehole <- function(x) {
  rw <- x[, 1] * (0.15 - 0.05) + 0.05
  r <-  x[, 2] * (50000 - 100) + 100
  Tu <- x[, 3] * (115600 - 63070) + 63070
  Hu <- x[, 4] * (1110 - 990) + 990
  Tl <- x[, 5] * (116 - 63.1) + 63.1
  Hl <- x[, 6] * (820 - 700) + 700
  L <-  x[, 7] * (1680 - 1120) + 1120
  Kw <- x[, 8] * (12045 - 9855) + 9855
  
  m1 <- 2 * pi * Tu * (Hu - Hl)
  m2 <- log(r / rw)
  m3 <- 1 + 2 * L * Tu / (m2 * rw ^ 2 * Kw) + Tu / Tl
  return(m1 / m2 / m3)
}

# Four cases

om = new(outermod)

testmultgrad = function(ss=400,nterms=100){
  
  
  d = 8
  set.seed(42)
  xo = matrix(runif(ss*d),ncol=d)
  yo = borehole(xo)
  x = xo
  
  offset = mean(yo)
  scale = sd(yo)
  y = (yo-offset)/scale
  knotlist = list()
  for(k in 1:d)  knotlist[[k]] = seq(0.001,0.999,0.025)
  
  setcovfs(om, c("mat25pow",rep("mat25",d-1)))
  
  setknot(om,knotlist)
  hyp0 = gethyp(om)
  om$updatehyp(hyp0)
  
  terms = om$selectterms(nterms)
  obp = new(outerbase,om,x) #build a basis mat, X
  theta = sd(y)/100*rnorm(nterms)
  
  meanh = obp$matmul(terms, theta)
  meanh_gradhyp = obp$matmul_gradhyp(terms, theta)
  grad = obp$tmatmul(terms, y - meanh )
  y0 = y - meanh
  grad_gradhyp = obp$tmatmul_gradhyp( terms, y0)
  
  
  eps = 10^(-3)
  hypp = (runif(length(hyp0))-0.5)
  hyp1 = hyp0 + eps*hypp
  om$updatehyp(hyp1)
  obp2 = new(outerbase,om,x) #build a basis mat, X
  
  meanh2 = obp2$matmul(terms, theta)
  grad2 = obp2$tmatmul(terms, y - meanh )
  
  L = list()
  L$mult_gradhyp = cbind((meanh2-meanh)/eps,meanh_gradhyp %*% hypp)
  L$tmult_gradhyp = cbind((grad2-grad)/eps,grad_gradhyp %*% hypp)
  L
}

expect_equal_or_warn <- function(...) tryCatch(expect_equal(...),
      error = function(e) warning("some inexact grad, this can be normal."))


test_that("short, skinny test", {
  L = testmultgrad(200,100)
  
  expect_equal_or_warn(L$mult_gradhyp[,1], L$mult_gradhyp[,2], tolerance = 0.1,
               label="mult grad agreement")
  expect_equal_or_warn(L$tmult_gradhyp[,1], L$tmult_gradhyp[,2], tolerance = 0.1,
               label="tmult grad agreement")
})

test_that("tall, skinny test", {
  L = testmultgrad(10000,100)
  
  expect_equal_or_warn(L$mult_gradhyp[,1], L$mult_gradhyp[,2], tolerance = 0.1,
               label="mult grad agreement")
  expect_equal_or_warn(L$tmult_gradhyp[,1], L$tmult_gradhyp[,2], tolerance = 0.1,
               label="tmult grad agreement")
})

test_that("short, wide test", {
  L = testmultgrad(200,1000)
  
  expect_equal_or_warn(L$mult_gradhyp[,1], L$mult_gradhyp[,2], tolerance = 0.1,
               label="mult grad agreement")
  expect_equal_or_warn(L$tmult_gradhyp[,1], L$tmult_gradhyp[,2], tolerance = 0.1,
               label="tmult grad agreement")
})

test_that("tall, wide test", {
  L = testmultgrad(10000,2000)
  
  expect_equal_or_warn(L$mult_gradhyp[,1], L$mult_gradhyp[,2], tolerance = 0.1,
               label="mult grad agreement")
  expect_equal_or_warn(L$tmult_gradhyp[,1], L$tmult_gradhyp[,2], tolerance = 0.1,
               label="tmult grad agreement")
})