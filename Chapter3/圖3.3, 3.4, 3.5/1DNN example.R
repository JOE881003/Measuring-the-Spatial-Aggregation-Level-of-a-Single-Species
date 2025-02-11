library(ggplot2)


## Random
set.seed(123)
x <- runif(10, 0, 500)
df_ra = data.frame(x = x, y = rep(0, length(x)))
#x11(20, 10)
ggplot(df_ra, aes(x = x, y = y))+
  geom_point(size = I(4))+
  geom_segment(aes(x = 0, y = 0, xend = 500, yend = 0), linewidth = I(1))+
  xlim(0, 500)+
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())



## Regular
x = seq(0, 500, 50)
df_re = data.frame(x = x, y = rep(0, length(x)))
#x11(20, 10)
ggplot(df_re, aes(x = x, y = y))+
  geom_point(size = I(4))+
  geom_line(linewidth = I(1))+
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())



## Aggregation
set.seed(123)
sig<-matrix(c(200,0,0,200),2,2)
x = mvrnorm(5,c(150,400),sig) 
x = c(x[,1], x[,2])

df_agg = data.frame(x = x, y = rep(0, length(x)))
#x11(20, 10)
ggplot(df_agg, aes(x = x, y = y))+
  geom_point(size = I(4))+
  geom_segment(aes(x = 0, y = 0, xend = 500, yend = 0), linewidth = I(1))+
  xlim(0, 500)+
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

