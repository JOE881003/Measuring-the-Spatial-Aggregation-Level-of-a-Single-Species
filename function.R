library(ggplot2)
library(MASS)
library(doParallel)
library(vegan)
library(ggpattern)
library(spatialEco)
library(magrittr)
library(spatstat)
library(reshape2)
library(scales)

for(i in 1:8){
  load(paste('./doi_10.15146_5xcp-0d46__v2/bci.tree/bci.tree', i,'.rdata', sep = ''))
}

BCI_List <- list()
BCI_List[[1]] <- bci.tree1
BCI_List[[2]] <- bci.tree2
BCI_List[[3]] <- bci.tree3
BCI_List[[4]] <- bci.tree4
BCI_List[[5]] <- bci.tree5
BCI_List[[6]] <- bci.tree6
BCI_List[[7]] <- bci.tree7
BCI_List[[8]] <- bci.tree8





get_elevation <- function(x, y, elevation_matrix, elve_data) {
  # 確保座標在有效範圍內
  if (x >= 0 & x <= max(elve_data$x) & y >= 0 & y <= max(elve_data$y)) {
    return(elevation_matrix[round(y/5)+1, round(x/5)+1])
  } else {
    return(NA)  # 若座標超出範圍則返回 NA
  }
}
#get_elevation(24, 10, elve_mat)


#BCI_data = BCI_List[[8]][which(BCI_List[[8]]$status=="A" & BCI_List[[8]]$dbh >= 10), 5:6]
#sp = unique(BCI_List[[8]]$sp)
species.distribution<-function(delta,lambda,prob,cnum,rnum)
{
  pnum<-rpois(1,lambda)+1 # avoiding zeros, Poission distribution,lambda =  parameter intensity ??
  offnum<-rgeom(pnum,prob)+1 # avoiding zeros,The Geometric Distribution, prob = 0.01
  # generate the locations of parental points
  x<-sample(cnum,pnum,replace=FALSE)
  y<-sample(rnum,pnum,replace=FALSE)
  #
  #now generating offspring locations
  sig<-matrix(c(delta,0,0,delta),2,2)
  pp<-vector()
  for(i in 1:pnum)
  {
    this<-mvrnorm(offnum[i],c(x[i],y[i]),sig)#bivariate normal distribution,delta link to ??2 
    this<-rbind(this,c(x[i],y[i])) #including parental point is tricky, otherwise sometimes this is NULL
    #exluding out-of-bound points
    idd<-which(this[,1]>cnum | this[,2]>rnum | this[,1]<1 | this[,2]<1)
    if(length(idd)!=0)
    {
      this<-this[-idd,]
    }
    pp<-rbind(pp,this)
  }
  #
  return(pp)
}#end

##################################################
#simulation of multinomial-Dirichlet distribution by avoiding NAs
#N is the total number of organisms, can be a single value or a vector
rMDir1<-function(n,alpha,N)
{
  p=rdirichlet1(n,alpha)
  #
  if(length(N)==1)
  {
    N=rep(N,1,n)
  }#
  #
  mat=vector()
  for(i in 1:n)
  {
    v=rmultinom(1,size=N[i],prob=p[i,])
    mat=rbind(mat,as.vector(v))
  }#
  #
  return(mat)
}#
#
#################
#avoid NAs:
rdirichlet1<-function (n = 1, alpha) 
{
  rv=rdirichlet(n,alpha)
  if(n==1)
  {
    rv=t(as.matrix(rv))
  }
  #
  ids=which(is.na(rv[,1]))
  if(length(ids)>0)
  {
    for(i in 1:length(ids))
    {
      rv[ids[i],]=rep(0,1,length(alpha))
      rv[ids[i],sample(1:length(alpha),1)]=1
    }#i
  }#
  #
  return(rv)
}
#################
##############
#random variate generation
rdirichlet<-function (n = 1, alpha) 
{
  Gam <- matrix(0, n, length(alpha))
  for (i in 1:length(alpha)) Gam[, i] <- rgamma(n, shape = alpha[i])
  Gam/rowSums(Gam)
}


#################
quad.num <- function(x, y, diste, xr = 1000, yr = 500){
  rlt <- c()
  xn <- xr/diste#50
  yn <- yr/diste#25
  
  
  for (i in 0:(yn-2)){
    for (j in 0:(xn-2)){
      rlt <- c(rlt, length(which(x >= j*diste & x < (j+1)*diste & y >= i*diste & y < (i+1)*diste)))
    }
    rlt <- c(rlt, length(which(x >= (xn-1)*diste & x <= xn*diste & y >= i*diste & y < (i+1)*diste)))
  }
  
  for (j in 0:(xn-2)){
    rlt <- c(rlt, length(which(x >= j*diste & x < (j+1)*diste & y >= (yn-1)*diste & y <= yn*diste)))
  }
  rlt <- c(rlt, length(which(x >= (xn-1)*diste & x <= (xn)*diste & y >= (yn-1)*diste & y <= yn*diste)))
  
  return(rlt)
}



#####Shannon of evenness index
D.SH <- function(data){
  p <- data/sum(data)
  s <- length(data)
  d1 <- sum(p*log(p))
  d2 <- log(s)
  return(list(shannon.evennes = -d1/d2, shannon = -d1))
}
#####Gini-simpsom of evenness index
D.GS <- function(data){
  p <- data/sum(data)
  s <- length(data)
  d1 <- s*(1-sum(p^2))
  d2 <- s-1
  return(list(gini.simpson.evennes = d1/d2, gini.simpson = d1/s))
}
#####reciprocal simpson of evenness index
D.IS <- function(data){
  p <- data/sum(data)
  s <- length(data)
  d1 <- 1
  d2 <- s*sum(p^2)
  return(d1/d2)
}
#####CCV
ccv <- function(xr, yr, r, type){
  lst <- c()
  ccv.size <- c(50, 40, 25, 20, 10, 5)
  for (i in ccv.size){
    dat = as.ppp(r, c(0, xr, 0, yr))
    dat <- quadratcount(dat, xr/i, yr/i)
    if (type == 'shannon'){
      lst <- c(lst, D.SH(dat[dat != 0])$shannon.evennes)
    }else if (type == 'gini-simpson'){
      lst <- c(lst, D.GS(dat[dat != 0])$gini.simpson.evennes)
    }else if (type == 'reciprocal'){
      lst <- c(lst, D.IS(dat[dat != 0]))
    }
  }
  return(sd(lst)/mean(lst))
}


bubble<-function(ssmat,factor=1)
{
  spsn=dim(ssmat)[1]
  siten=dim(ssmat)[2]
  mat=ssmat/sum(ssmat)
  mat=mat/max(mat)
  ##########################
  xy=expand.grid(1:spsn,1:siten)
  # windows()
  par(cex.axis=1.5,cex.lab=1.5)
  plot(xy,xlim=c(1,(spsn+1)),ylim=c(1,(siten+1)),axes=FALSE,xlab=expression(alpha),ylab="Sites",type="n")
  #axis(1,c(1:spsn)+.5,1:spsn,tck=FALSE,col=NA)
  axis(2,c(1:siten)+.5,1:siten,tck=FALSE,col=NA)
  ##########################
  for(i in 1:siten)
  {
    id=which(xy[,2]==i)
    z=mat[,i]
    points(xy[id,]+.5,cex=z*factor*1.5,col=adjustcolor("red",alpha.f=.5),pch=16)
  }#i
  ##########################
  for(i in 1:(siten+1))
  {
    lines(1:(spsn+1),rep(i,1,spsn+1),col="gray")
  }#
  #
  for(i in 1:(spsn+1))
  {
    lines(rep(i,1,siten+1),1:(siten+1),col="gray")
  }#
  ##########################
  #
}#end

### The Cluster-Random-Regular Continuum
#  Clark, P. J., & Evans, F. C. (1954). Distance to nearest neighbor as a measure of spatial relationships in populations. Ecology, 35, 445-453.
#  R is the Clark and Evans competition index
#  the input mat contains x and y coordinates 

CECI <- function(mat, area, method = "P")
{
  # if(length(x) != length(y)){
  #   stop("x and y should have the same number")	
  # }
  n = dim(mat)[1] # numbers of points
  area = area   # study area (square meter)
  Rho = n/area  # density of points per square meter
  # minimal distance between a target point to its neighbours
  if(method == "NP"){
    dd <- as.matrix(dist(mat, p=2))
    diag(dd) <- NA
    dmin = apply(dd, 1, min, na.rm = TRUE)
  } else if (method == "P"){
    
    t1 <- Sys.time()
    library(doParallel)
    detectCores()
    cl <- makeCluster(10)
    registerDoParallel(cl)
    Res <- foreach(x = 1:n,.combine = c) %dopar% {
      d <- sqrt((mat[x,1] - mat[,1])^2 + (mat[x,2] - mat[,2])^2)
      dmin = sort(d)[2]
      dmin
    }
    t2 <- Sys.time()-t1
    t2
    stopCluster(cl)
    dmin = Res
  }
  # R is the Clark and Evans competition index
  R = (1/n) * sum(dmin) * 2 * sqrt(Rho)
  names(R) = "Clark & Evans Competition Index"
  
  return(list(R = R, t1 = (1/n) * sum(dmin), t2 = 1/(2*sqrt(Rho))))
}



#########Find nearest neighbor distance in line
roll <- function(data, d){
  data <- as.data.frame(data)
  n <- dim(data)[1]
  if ((d <= n) & (d > 0)){
    re = rbind(data[-(1:d),,drop = FALSE], data[1:d,])
  }else if((d <= n) & (d < 0)){
    re = rbind(data[(n+d+1):n,], data[-((n+d+1):n),,drop = FALSE])
  }
  
  return(re)
}



#########Other way to find nearest neighbor distance in line
NN_1d <- function(data){
  n = length(data)
  if(n > 1){
  data <- sort(data)
  da <- data.frame(data, roll(data, 1), roll(data, -1))
  d <- apply(da, 1, function(x) min(abs(x[1]-x[-1])))
  t1 <- mean(d)
    return(t1)
  }else {
    return(0)
  }
}

#########Other way to find nearest neighbor distance in line
NNL <- function(data, L){
  n <- length(data)
  data <- sort(data)
  ra <- NN_1d(data)
  lam <- (n-1)/L
  re <- 1/(2*lam)
  sig <- 1/(2*lam*sqrt(n))
  c <- (ra-re)/sig
  return(list(NNL = ra/re, ra = ra, re = re, p = 2*pt(-abs(c/2), df = n-1)))
}

#########NN index use C
nni <- function(x, win = c("hull", "extent"), method = 'C', A = 1000*500) {
  if(!inherits(x, "sf"))		
    stop(deparse(substitute(x)), " must be an sf POINT object")	
  if(unique(as.character(sf::st_geometry_type(x))) != "POINT")
    stop(deparse(substitute(x)), " must be an sf POINT object")		
  if (win[1] == "hull") {
    w <- spatstat.geom::convexhull.xy( sf::st_coordinates(x)[,1:2] )
  }
  if (win[1] == "extent") {
    e <- as.vector(sf::st_bbox(x))
    w <- spatstat.geom::as.owin(c(e[1], e[3], e[2], e[4]))
  }
  x <- spatstat.geom::as.ppp(sf::st_coordinates(x)[,1:2], w)
  A <- A
  obsMeanDist <- sum(spatstat.geom::nndist(x, method = method))/x$n
  expMeanDist <- 0.5 * sqrt(A / x$n)
  se <- 0.26136 / ((x$n**2.0 / A)**0.5)
  nni <- obsMeanDist / expMeanDist
  z <- (obsMeanDist - expMeanDist) / se
  return(list(NNI = nni, z.score = z, p = 2*stats::pnorm(-abs(z)),  
              expected.mean.distance = expMeanDist,
              observed.mean.distance = obsMeanDist))
}




#########QNN use C 
QNN_C <- function(data, xr, yr, vr, hr = vr){
  mat_t1 = c()
  mat_t2 = c()
  da <- data.frame(x = data[,1], y = data[,2], xcut = cut(data[,1], breaks = seq(0, xr, vr), labels = 1:(xr/vr)), ycut = cut(data[,2], breaks = seq(0, yr, hr), labels = 1:((yr/hr))))
  for(i in 1:(xr/vr)){
    dat = sort(da[which(da$xcut == i),2])
    n = length(dat)
    if (n <= 1) {
      mat_t1 <- c(mat_t1, 0)
      mat_t2 <- c(mat_t2, 0)
    }else{
      mat_t1 <- c(mat_t1, mean(nndist(dat))*n)
      mat_t2 <- c(mat_t2, (yr*n)/(2*(n-1)))
    }
  }
  
  for(i in 1:((yr/hr))){
    dat = sort(da[which(da$ycut == i),1])
    n = length(dat)
    if(n <= 1){
      mat_t1 <- c(mat_t1, 0)
      mat_t2 <- c(mat_t2, 0)
    }else{
      mat_t1 <- c(mat_t1, mean(nndist(dat))*n)
      mat_t2 <- c(mat_t2, (xr*n)/(2*(n-1)))
    }
  }
  
  return(list(ar = sum(mat_t1)/sum(mat_t2), t1 = sum(mat_t1), t2 = sum(mat_t2)))
}


###########QNN use R
QNN <- function(data, xr, yr, vr, hr = vr){
  da <- data.frame(x = data[,1], y = data[,2], xcut = cut(data[,1], breaks = seq(0, xr, vr), labels = 1:(xr/vr)), ycut = cut(data[,2], breaks = seq(0, yr, hr), labels = 1:((yr/hr))))
  N <- dim(da)[1]
  x_t1_lst <- c()
  x_t2_lst <- c()
  y_t1_lst <- c()
  y_t2_lst <- c()
  sigre <- c()
  clst <- c()
  k = 0
  for(i in 1:(xr/vr)){
    dat = da[which(da$xcut == i),2]
    n <- length(dat)
    if(n <= 1){
      x_t1_lst[i] = 0
      x_t2_lst[i] = 0
    }else{
      k = k + 1
      lam = (n-1)/yr
      sig2 <- n*(1/(2*lam))^2
      t2 = (1)/(2*lam)
      NN1d = NN_1d(dat)
      x_t1_lst[i] = NN1d*(n)
      x_t2_lst[i]= t2*n
      sigre = c(sigre, sig2)
    }
  }
  
  for(i in 1:((yr/hr))){
    dat = da[which(da$ycut == i),1]
    n <- length(dat)
    if(n <= 1){
      y_t1_lst[i] = 0
      y_t2_lst[i] = 0
    }else{
      k = k + 1
      lam = (n-1)/xr
      sig2 <- n*(1/(2*lam))^2
      t2 = (1)/(2*lam)
      NN1d = NN_1d(dat)
      y_t1_lst[i] = NN1d*(n)
      y_t2_lst[i]= t2*n
      sigre = c(sigre, sig2)
    }
  }
  t1_lst <- c(x_t1_lst, y_t1_lst)
  t2_lst <- c(x_t2_lst, y_t2_lst)
  return(list(ar = sum(t1_lst)/sum(t2_lst), t1 = sum(t1_lst), t2 = sum(t2_lst), c = ((sum(t1_lst)/(N))-(sum(t2_lst)/(N)))/(sqrt(sum(sigre))/(N)), deg = k))
}



########Calculate traditional NN index, but use quadrat data
NNQ <- function(data, xr, yr, r){
  da <- data.frame(x = data[,1], y = data[,2], xcut = cut(data[,1], breaks = seq(0, xr, r), labels = 1:(xr/r)), ycut = cut(data[,2], breaks = seq(0, yr, r), labels = 1:((yr/r))))
  N <- dim(da)[1]
  mat <- c()
  
  for (i in 1:(xr/r)){
    for (j in 1:(yr/r)){
      dat <- da[which((da$xcut == i) & (da$ycut == j)),1:2]
      n <- dim(dat)[1]
      
      if (n <= 1){
        mat <- c(mat, 0)
      }else{
        
        mat <- c(mat, mean(nndist(dat))*2*sqrt(n/r^2)*n)
      }
      
    }
  }
  
  return(list(R = sum(mat)/N))
}




##############################################################################
#probability mass (or density) function of the multinomial-Dirichlet distribution
dMDir<-function(x,alpha,log=FALSE)
{
  x=as.vector(x)
  if(length(x)==length(alpha))
  {
    num=length(x)
    n=sum(x)
    #
    p1=lgamma(n+1)
    p2=lgamma(sum(alpha))
    p3=lgamma(n+sum(alpha))
    #
    p4=sum(lgamma(x+alpha))
    p5=sum(lgamma(x+1))
    p6=sum(lgamma(alpha))
    #
    f=p1+p2-p3+p4-p5-p6
    #
    if(log==FALSE)
    {
      return(exp(f))
    }else
    {
      return(f)
    }
    #
  }else
  {
    return(NA)
  }
}#
##################################################
######################################################
#fitting of the SDM model
#mat is a species-site matrix
#if it is a vector, convert it into a matrix
fit<-function(mat)
{
  
  #mat = mat[which(mat > 0)]
  
  if(is.vector(mat))
  {
    mat=t(as.matrix(mat))
  }#
  ##################
  n=dim(mat)[2] #number of sites
  likelihood<-function(pars)
  {
    alpha=pars[1]
    alpha=rep(alpha,1,n)
    #
    one=function(v)
    {
      return(dMDir(v,alpha,log=TRUE))
    }#
    #
    all=apply(mat,1,one)
    return(-sum(all))
  }#
  ##################
  res=nlminb(runif(1),likelihood,lower=1e-20,upper=1e+20)
  #
  return(res)
}#
#
###################Timing

timing <- function(x){
  t1 = Sys.time()
  x
  t2 = Sys.time()
  return(difftime(t2, t1, units = "secs")[[1]])
}


####
# DC means deviation coefficient or diffusion coefficient 
# G. E. BLACKMAN, Statistical and Ecological Studies in the Distribution of Species in Plant Communities: I. Dispersion as a Factor in the Study of Changes in Plant Populations, Annals of Botany, Volume 6, Issue 2, April 1942, Pages 351?C370, https://doi.org/10.1093/oxfordjournals.aob.a088411
DC <- function(dat)
{ 
  N  = length(dat)
  dat <- array(dat)
  dat <- data.frame(table(dat))
  dat$dat <- as.numeric(as.character(dat$dat))
  #variance
  V = (sum(dat$dat^2 * dat$Freq) - (((sum(dat$dat * dat$Freq))^2)/N))/(N-1)
  #mean
  X = sum(dat$dat * dat$Freq)/N
  # deviation coefficient
  dc = V/X
  # names(dc) = "Diffusion coefficient"
  # S is the standard error
  # t value
  t = (dc - 1) / sqrt(2/(N - 1))
  # p value 
  p_value = round(2 * pt(-abs(t),df = N-1, lower.tail = T),3)
  # p_value = 2 * pt(abs(t),df = N-1, lower.tail = F)
  return(list("dc" = dc,"t" = t, "p" = p_value))
}
######################################################
#fitting of the independent NBD model
#mat is a species-site matrix
#if it is a vector, convert it into a matrix
fitNBD<-function(mat)
{
  x=as.vector(mat) #use vector data for fitting
  #x = x[which(x > 0)]
  likelihood<-function(pars)
  {
    k=pars[1]
    mu=pars[2]
    #
    all=dnbinom(x,size=k,mu=mu,log=TRUE)
    return(-sum(all,na.rm=TRUE))
  }#
  ##################
  res=nlminb(runif(2),likelihood,lower=rep(1e-20,1,2),upper=rep(1e+20,1,2))
  #
  return(res)
}#
#
