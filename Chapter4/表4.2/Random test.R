source("../function.R")

set.seed(1234)
NN_random = c()
QNN_random = c()
N = c(200, 500, 1000, 5000)
xr = 1000
yr = 500
for (n in N){
  r = (-(xr+yr)-sqrt((xr+yr)^2-4*(1-n)*xr*yr))/(2*(1-n))
  maplst <- list()
  for (i in 1:1000){
    df = data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 500))
    maplst[[i]] = df
  }
  
  NN_p <- c()
  for (i in 1:1000){
    meuse = sf::st_as_sf(maplst[[i]], coords = c("x", "y"),  
                         crs = 28992, agr = "constant")
    p = nni(meuse)$p
    NN_p = c(NN_p, p)
  }
  NN_p[which(NN_p < 0.05)] = 'aggregation'
  NN_p[which(NN_p != 'aggregation')] = 'random'
  NN_random = c(NN_random, length(NN_p[NN_p == 'random'])/length(NN_p))
  
  
  AR_p <- c()
  for (i in 1:1000){
    data = maplst[[i]]
    t1 = QNN(data, 1000, 500, r)
    p = 2*stats::pnorm(-abs(t1$c))
    AR_p = c(AR_p, p)
  }
  AR_p[which(AR_p < 0.05)] = 'aggregation'
  AR_p[which(AR_p != 'aggregation')] = 'random'
  QNN_random = c(QNN_random, length(AR_p[AR_p == 'random'])/length(AR_p))
  
}

result = data.frame(N = N, QNN = QNN_random, NN = NN_random)

result