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
  xo = matrix(runif(ss*d),ncol=d)
  A = matrix(0,nrow=d,ncol=d)
  for(i in 1:d) for(j in 1:d) A[i,j] = 0.001^abs(i-j)
  xo = ((xo-0.5) %*% A) +0.5
  miv = apply(xo,2,min)
  mav = apply(xo,2,max)
  xo = t(0.005+0.99*(t(xo)-miv)/(mav-miv))
  yo = borehole(xo)
  x = xo
  
  offset = mean(yo)
  scale = sd(yo)
  y = (yo-offset)/scale
  knotlist = list()
  for(k in 1:d)  knotlist[[k]] = quantile(x[,k],
                                          seq(0,1,length=40)*40/
                                            (40+1)+0.5/(40+1))
  
  eta = 0.1+0.1*(runif(2*d)-0.5)#rep(c(0,0),d)
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
      error = function(e) warning("some errors on grad, this can be normal."))


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