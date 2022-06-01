#' Three demensional Borehole
#' 
#' @param x a n by 3 vector of inputs
#' @return a length n vector of outputs
#' @export
obtest_borehole3d <- function(x) {
  rw <- x[, 1] * (0.15 - 0.05) + 0.05
  Hl <- x[, 2] * (820 - 700) + 700
  L <-  x[, 3] * (1680 - 1120) + 1120
  
  r <-  0.5 * (50000 - 100) + 100
  Tu <- 0.5 * (115600 - 63070) + 63070
  Hu <- 0.5 * (1110 - 990) + 990
  Tl <- 0.5 * (116 - 63.1) + 63.1
  Kw <- 0.5 * (12045 - 9855) + 9855
  
  m1 <- 2 * pi * Tu * (Hu - Hl)
  m2 <- log(r / rw)
  m3 <- 1 + 2 * L * Tu / (m2 * rw ^ 2 * Kw) + Tu / Tl
  return(m1 / m2 / m3 - 77) #approximately centering
}

#' Eight demensional Borehole
#' 
#' @param x a n by 8 vector of inputs
#' @return a length n vector of outputs
#' @export
obtest_borehole8d <- function(x) {
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
  return((m1 / m2 / m3-77))
}