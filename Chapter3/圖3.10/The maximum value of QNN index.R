library(ggplot2)

xr = 1000
yr = 500

r = 1000/10
h = 500/5

x = rep(seq(0, xr , r), 10)
y = t = rep(yr, xr/r+1)
for (i in 1:9){
  t = t-h
  y = c(y, t)
}
df_re <- data.frame(x = x, y = y)
df_re <- df_re[which(df_re$y >= 0),]
x11(20, 10)
ggplot(data = df_re, aes(x = x, y = y))+
  geom_point(size = I(3))+
  xlim(0, 1000)+
  scale_x_continuous(limits=c(0,1000), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,500), expand=c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(.4, .7, , , "cm"))



