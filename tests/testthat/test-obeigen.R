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
library(outerbase)

om = new(outermod)

ss=2000
nterms = 500
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
ob = new(outerbase,om,x) #build a basis mat, X
theta = sqrt(om$getvar(terms)/20)*rnorm(nterms)

basemat_getbase = matrix(1,ncol=nterms,nrow=ss)
for(k in 1:8){
  Bh = ob$getbase(k)
  basemat_getbase = basemat_getbase * Bh[,terms[,k]+1]
}
basemat_getmat = ob$getmat(terms)

getmatgetbasediff = sum(abs(basemat_getmat-basemat_getbase))

vec1_getmat = basemat_getmat %*% theta
vec1_matmul = ob$matmul(terms,theta)

plot(vec1_getmat,vec1_matmul)

vec1_matmulEigen = ob$matmulEigen(terms,theta)
plot(vec1_getmat,vec1_matmulEigen)

library(microbenchmark)

microbenchmark(ob$matmul(terms,theta), ob$matmulEigen(terms,theta))