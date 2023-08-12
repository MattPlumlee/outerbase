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
d = 8

om = new(outermod)

getvals = function(obj, coeff, hyp, para) {
  obj$compute_gradhyp =TRUE
  obj$compute_gradpara = TRUE
  om$updatehyp(hyp)
  obj$updateom()
  obj$updatepara(para)
  obj$update(coeff)
  L = list()
  L$val = obj$val
  L$grad = obj$grad
  L$gradhyp = obj$gradhyp
  L$gradpara = obj$gradpara
  L$diaghess = obj$diaghess()
  L$diaghessgradhyp = obj$diaghessgradhyp()
  L$diaghessgradpara = obj$diaghessgradpara()
  L
}

perturbvals = function(obj, coeff, coeffp,
                       hyp, hypp,
                       para, parap, ep = 10^(-4)) {
  
  L = getvals(obj, coeff, hyp, para)
  L_coeff = getvals(obj, coeff+ep*coeffp, 
                    hyp, para) 
  L_para = getvals(obj, coeff, 
                  hyp, para+ep*parap) 
  L_hyp = getvals(obj, coeff, 
                  hyp+ep*hypp, para) 
  L2 = getvals(obj, coeff, hyp, para)
  
  L$grad = (L$grad+L2$grad+L_coeff$grad)/3
  L$gradhyp = (L$gradhyp+L2$gradhyp+L_para$gradhyp)/3
  L$gradpara = (L$gradpara+L2$gradpara+L_para$gradpara)/3
  
  gradcheck = list()
  gradcheck$coeff = c((L_coeff$val - L$val),ep*sum((L$grad*coeffp)))
  gradcheck$hyp = c((L_hyp$val - L$val),ep*sum((L$gradhyp*hypp)))
  gradcheck$para =c( (L_para$val - L$val),ep*sum((L$gradpara*parap)))
  gradcheck$rep = c(L$val,L2$val)
  gradcheck$dh = c(L$diaghess,L2$diaghess)
  rvec = rnorm(length(L$diaghess))
  gradcheck$dhcoeff = cbind((L_hyp$diaghess - L$diaghess),
                            ep*L$diaghessgradhyp %*% hypp)
  gradcheck$dhpara = cbind((L_para$diaghess - L$diaghess),
                          ep*L$diaghessgradpara %*% parap)
  gradcheck
}



fulltest = function(ss, nterms, corr=0, ep = 10^(-2)){
  
  d = 8
  
  set.seed(42)
  xo = matrix(runif(ss*d),ncol=d)
  yo = borehole(xo)
  x = xo
  
  offset = mean(yo)
  scale = sd(yo)
  d = dim(x)[2]
  y = (yo-offset)/scale
  setcovfs(om, rep("mat25",d))
  knotlist = list()
  for(k in 1:d)  knotlist[[k]] = seq(0.001,0.999,0.05)
  
  hyp = gethyp(om)
  prpara = log(1)
  likpara = c(log(0.1))
  
  
  setknot(om,knotlist)
  om$updatehyp(hyp)
  terms = om$selectterms(nterms)
  
  logpr = new(logpr_gauss,om,terms)
  loglik = new(loglik_gauss,om,terms,y,x)
  
  coeff = sd(y)/100*rnorm(nterms)
  coeffp =sd(y)/100*rnorm(nterms)
  
  hypp = (runif(length(hyp))-0.5)
  likparap = rnorm(length(likpara))
  prparap = rnorm(length(prpara))
  
  pinfo = list()
  pinfo$prinfo = perturbvals(logpr, coeff, coeffp,
                             hyp, hypp,
                             prpara, prparap,ep=ep)
  
  pinfo$likinfo = perturbvals(loglik, coeff, coeffp,
                              hyp, hypp,
                              likpara, likparap,ep=ep)
  
  logpdf = new(lpdfvec,loglik,logpr)
  vecpara = getpara(logpdf)
  vecparap = rnorm(length(vecpara))
  
  
  pinfo$vecinfo = perturbvals(logpdf, coeff, coeffp,
                              hyp, hypp,
                              vecpara, vecparap,ep=ep)
  
  pinfo
}

expect_equal_or_warnh <- function(...) tryCatch(expect_equal(...),
        error = function(e) warning("some inexact grad, this can be normal."))


test_that("short, skinny test", {
  ep = 10^(-4)
  pinfo = fulltest(200,100, ep = ep)
  
  expect_equal_or_warnh(pinfo$prinfo$coeff[1], pinfo$prinfo$coeff[2], tolerance = 1,
               label="pr grad agreement")
  expect_equal_or_warnh(pinfo$prinfo$hyp[1], pinfo$prinfo$hyp[2], tolerance = 1,
               label="pr gradhyp agreement")
  expect_equal_or_warnh(pinfo$prinfo$para[1], pinfo$prinfo$para[2], tolerance = 1,
               label="pr gradpara agreement")
  
  expect_equal_or_warnh(pinfo$likinfo$coeff[1], pinfo$likinfo$coeff[2], tolerance = 1,
               label="lik grad agreement")
  expect_equal_or_warnh(pinfo$likinfo$hyp[1], pinfo$likinfo$hyp[2], tolerance = 1,
               label="lik gradhyp agreement")
  expect_equal_or_warnh(pinfo$likinfo$para[1], pinfo$likinfo$para[2], tolerance = 1,
               label="lik gradpara agreement")
  
  expect_equal_or_warnh(pinfo$vecinfo$coeff[1], pinfo$vecinfo$coeff[2], tolerance = 1,
                        label="lik grad agreement")
  expect_equal_or_warnh(pinfo$vecinfo$hyp[1], pinfo$vecinfo$hyp[2], tolerance = 1,
                        label="lik gradhyp agreement")
  expect_equal_or_warnh(pinfo$vecinfo$para[1], pinfo$vecinfo$para[2], tolerance = 1,
                        label="lik gradpara agreement")
  
  expect_equal_or_warnh(pinfo$likinfo$dhcoeff[,1], 
               pinfo$likinfo$dhcoeff[,2], tolerance = 1,
               label="lik dhgradcoeff agreement")
  expect_equal_or_warnh(pinfo$likinfo$dhpara[,1], 
               pinfo$likinfo$dhpara[,2], tolerance = 1,
               label="lik dhgradpara agreement")
  
  expect_equal_or_warnh(pinfo$likinfo$rep[1], 
               pinfo$likinfo$rep[2], tolerance = 10^(-5),
               label="rep agreement")
})



test_that("tall, skinny test", {
  ep = 10^(-4)
  pinfo = fulltest(10000,100, ep = ep)
  
  expect_equal_or_warnh(pinfo$prinfo$coeff[1], pinfo$prinfo$coeff[2], tolerance = 1,
               label="pr grad agreement")
  expect_equal_or_warnh(pinfo$prinfo$hyp[1], pinfo$prinfo$hyp[2], tolerance = 1,
               label="pr gradhyp agreement")
  expect_equal_or_warnh(pinfo$prinfo$para[1], pinfo$prinfo$para[2], tolerance = 1,
               label="pr gradpara agreement")
  
  expect_equal_or_warnh(pinfo$likinfo$coeff[1], pinfo$likinfo$coeff[2], tolerance = 1,
               label="lik grad agreement")
  expect_equal_or_warnh(pinfo$likinfo$hyp[1], pinfo$likinfo$hyp[2], tolerance = 1,
               label="lik gradhyp agreement")
  expect_equal_or_warnh(pinfo$likinfo$para[1], pinfo$likinfo$para[2], tolerance = 1,
               label="lik gradpara agreement")
  
  expect_equal_or_warnh(pinfo$vecinfo$coeff[1], pinfo$vecinfo$coeff[2], tolerance = 1,
                        label="lik grad agreement")
  ###
  expect_equal_or_warnh(pinfo$vecinfo$hyp[1], pinfo$vecinfo$hyp[2], tolerance = 1,
                        label="lik gradhyp agreement")
  expect_equal_or_warnh(pinfo$vecinfo$para[1], pinfo$vecinfo$para[2], tolerance = 1,
                        label="lik gradpara agreement")
  
  expect_equal_or_warnh(pinfo$likinfo$dhcoeff[,1], 
               pinfo$likinfo$dhcoeff[,2], tolerance = 1,
               label="lik dhgradcoeff agreement")
  expect_equal_or_warnh(pinfo$likinfo$dhpara[,1], 
               pinfo$likinfo$dhpara[,2], tolerance = 1,
               label="lik dhgradpara agreement")
  
  expect_equal_or_warnh(pinfo$likinfo$rep[1], pinfo$likinfo$rep[2], tolerance = 10^(-4),
               label="rep agreement")
})


test_that("short, wide test", {
  ep = 10^(-4)
  pinfo = fulltest(200,1000, ep = ep)
  
  
  expect_equal_or_warnh(pinfo$prinfo$coeff[1], pinfo$prinfo$coeff[2], tolerance = 1,
               label="pr grad agreement")
  expect_equal_or_warnh(pinfo$prinfo$hyp[1], pinfo$prinfo$hyp[2], tolerance = 1,
               label="pr gradhyp agreement")
  expect_equal_or_warnh(pinfo$prinfo$para[1], pinfo$prinfo$para[2], tolerance = 1,
               label="pr gradpara agreement")
  
  expect_equal_or_warnh(pinfo$likinfo$coeff[1], pinfo$likinfo$coeff[2], tolerance = 1,
               label="lik grad agreement")
  expect_equal_or_warnh(pinfo$likinfo$hyp[1], pinfo$likinfo$hyp[2], tolerance = 1,
               label="lik gradhyp agreement")
  expect_equal_or_warnh(pinfo$likinfo$para[1], pinfo$likinfo$para[2], tolerance = 1,
               label="lik gradpara agreement")
  
  
  expect_equal_or_warnh(pinfo$vecinfo$coeff[1], pinfo$vecinfo$coeff[2], tolerance = 1,
                        label="lik grad agreement")
  expect_equal_or_warnh(pinfo$vecinfo$hyp[1], pinfo$vecinfo$hyp[2], tolerance = 1,
                        label="lik gradhyp agreement")
  expect_equal_or_warnh(pinfo$vecinfo$para[1], pinfo$vecinfo$para[2], tolerance = 1,
                        label="lik gradpara agreement")
  
  expect_equal_or_warnh(pinfo$likinfo$dhcoeff[,1], 
               pinfo$likinfo$dhcoeff[,2], tolerance = 1,
               label="lik dhgradcoeff agreement")
  expect_equal_or_warnh(pinfo$likinfo$dhpara[,1], 
               pinfo$likinfo$dhpara[,2], tolerance = 1,
               label="lik dhgradpara agreement")
  
  expect_equal_or_warnh(pinfo$likinfo$rep[1], pinfo$likinfo$rep[2], tolerance = 10^(-4),
               label="rep agreement")
})


test_that("tall, wide test", {
  ep = 10^(-4)
  pinfo = fulltest(10000,2500, ep = ep)
  
  
  expect_equal_or_warnh(pinfo$prinfo$coeff[1], pinfo$prinfo$coeff[2], tolerance = 1,
               label="pr grad agreement")
  expect_equal_or_warnh(pinfo$prinfo$hyp[1], pinfo$prinfo$hyp[2], tolerance = 1,
               label="pr gradhyp agreement")
  expect_equal_or_warnh(pinfo$prinfo$para[1], pinfo$prinfo$para[2], tolerance = 1,
               label="pr gradpara agreement")
  
  expect_equal_or_warnh(pinfo$likinfo$coeff[1], pinfo$likinfo$coeff[2], tolerance = 1,
               label="lik grad agreement")
  expect_equal_or_warnh(pinfo$likinfo$hyp[1], pinfo$likinfo$hyp[2], tolerance = 10,
               label="lik gradhyp agreement")
  expect_equal_or_warnh(pinfo$likinfo$para[1], pinfo$likinfo$para[2], tolerance = 1,
               label="lik gradpara agreement")
  
  expect_equal_or_warnh(pinfo$vecinfo$coeff[1], pinfo$vecinfo$coeff[2], tolerance = 1,
                        label="lik grad agreement")
  expect_equal_or_warnh(pinfo$vecinfo$hyp[1], pinfo$vecinfo$hyp[2], tolerance = 1,
                        label="lik gradhyp agreement")
  expect_equal_or_warnh(pinfo$vecinfo$para[1], pinfo$vecinfo$para[2], tolerance = 1,
                        label="lik gradpara agreement")
  
  expect_equal_or_warnh(pinfo$likinfo$dhcoeff[,1], 
               pinfo$likinfo$dhcoeff[,2], tolerance = 1,
               label="lik dhgradcoeff agreement")
  expect_equal_or_warnh(pinfo$likinfo$dhpara[,1], 
               pinfo$likinfo$dhpara[,2], tolerance = 1,
               label="lik dhgradpara agreement")
  
  expect_equal_or_warnh(pinfo$likinfo$rep[1], pinfo$likinfo$rep[2], tolerance = 10^(-4),
               label="rep agreement")
})