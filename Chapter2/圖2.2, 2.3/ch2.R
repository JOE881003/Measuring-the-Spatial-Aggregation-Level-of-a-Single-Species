source('../function.R')

##############The maximum value of NN index
oh = ((5-sqrt(3))/2)
x <- c(seq(2, 3, 1), seq(1.5, 3.5, 1), seq(2, 3, 1))
y <- c(rep(oh+sqrt(3), 2), rep(oh+sqrt(3)/2, 3), rep(oh, 2))
#x11(10, 10)
plot(x, y, xlim = c(0, 5), ylim = c(0, 5), pch = 19, xlab = '', ylab = '', xaxt='n', yaxt = 'n')


polygon(c(1, 2, 1.5, 1, 2), c(oh, oh, oh+sqrt(3)/2, oh+sqrt(3), oh+sqrt(3)),
        col = "blue",
        density = 5, angle = 45, lwd = 3)

polygon(c(1.5, 2.5, 2, 1.5, 2.5), c(oh-sqrt(3)/2, oh-sqrt(3)/2, oh, oh+sqrt(3)/2, oh+sqrt(3)/2),
        col = "peru",
        density = 5, angle = 45, lwd = 3)

polygon(c(2.5, 3.5, 3, 2.5, 3.5), c(oh-sqrt(3)/2, oh-sqrt(3)/2, oh, oh+sqrt(3)/2, oh+sqrt(3)/2),
        col = "violetred",
        density = 5, angle = 45, lwd = 3)

polygon(c(2, 3, 2.5, 2, 3), c(oh, oh, oh+sqrt(3)/2, oh+sqrt(3), oh+sqrt(3)),
        col = "green4",
        density = 5, angle = 45, lwd = 3)

polygon(c(3, 4, 3.5, 3, 4), c(oh, oh, oh+sqrt(3)/2, oh+sqrt(3), oh+sqrt(3)),
        col = "cyan3",
        density = 5, angle = 45, lwd = 3)

polygon(c(1.5, 2.5, 2, 1.5, 2.5), c(oh+sqrt(3)/2, oh+sqrt(3)/2, oh+sqrt(3), oh+1.5*sqrt(3), oh+1.5*sqrt(3)),
        col = "gray64",
        density = 5, angle = 45, lwd = 3)

polygon(c(2.5, 3.5, 3, 2.5, 3.5), c(oh+sqrt(3)/2, oh+sqrt(3)/2, oh+sqrt(3), oh+1.5*sqrt(3), oh+1.5*sqrt(3)),
        col = "orange",
        density = 5, angle = 45, lwd = 3)

points(x, y, cex = 3, pch=19)



################Different alpha in eight quadrat
set.seed(3)
dat1=rMDir1(8,rep(.5,1,10),sample(500:1000,8))#When alpha is 0.5
set.seed(4)
dat2=rMDir1(8,rep(1,1,10),sample(500:1000,8))#When alpha is 1
set.seed(1123)
dat3=rMDir1(8,rep(10,1,10),sample(500:1000,8))#When alpha is 10

dat = cbind(dat1[,1], dat2[,1], dat3[,1])

#x11(15, 10)
bubble(t(dat),factor=5)
axis(1,c(1:3)+.5, c(0.5, 1, 10),tck=FALSE,col=NA)

