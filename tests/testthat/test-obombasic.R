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

testsimple = function(ss=400){
  
  nterms = 20
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
  ob = new(outerbase,om,x) #build a basis mat, X
  theta = sqrt(om$getvar(terms)/20)*rnorm(nterms)
  
  basemat_getbase = matrix(1,ncol=nterms,nrow=ss)
  for(k in 1:8){
    basemat_getbase = basemat_getbase*ob$getbase(k)[,terms[,k]+1]
  }
  basemat_getmat = ob$getmat(terms)
  
  getmatgetbasediff = sum(abs(basemat_getmat-basemat_getbase))
  vec1_getmat = basemat_getmat %*% theta
  vec1_matmul = ob$matmul(terms,theta)
  getmatmatmuldiff = sum(abs(vec1_getmat-vec1_matmul))
  vec2_getmat = t(basemat_getmat) %*% vec1_getmat
  vec2_tmatmul = ob$tmatmul(terms,vec1_getmat)
  getmattmatmuldiff = sum(abs(vec1_getmat-vec1_matmul))
  
  Lr = list()
  Lr$getbase = getmatgetbasediff 
  Lr$getmat = getmatgetbasediff 
  Lr$matmul = getmatmatmuldiff
  Lr$tmatmul = getmattmatmuldiff
  Lr
}

test_that("basic test", {
  L = testsimple(500)
  expect_equal(L$getbase, 0, tolerance = 0.01,
               label="base agreement")
  expect_equal(L$getmat, 0, tolerance = 0.01,
               label="mat grad agreement")
  expect_equal(L$matmul, 0, tolerance = 0.01,
               label="getmat grad agreement")
  expect_equal(L$tmatmul, 0, tolerance = 0.01,
               label="tmult grad agreement")
})
