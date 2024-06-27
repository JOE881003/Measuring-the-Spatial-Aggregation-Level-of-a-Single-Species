library(ggplot2)
set.seed(123)
df = data.frame(x = runif(10, 0, 1000), y = runif(10, 0, 500))

######Vertical
ggplot(df, aes(x = x, y = y))+
  geom_point(size = I(4))+
  ylim(0, 500)+
  xlim(0, 1000)+
  geom_vline(xintercept = seq(0, 1000, 100), linewidth = I(1))+
  scale_x_continuous(limits=c(0,1000), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,500), expand=c(0,0)) +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.minor = element_blank())



######Horizontal
ggplot(data = df, aes(x = x, y = y))+
  geom_point(size = I(4))+
  ylim(0, 500)+
  xlim(0, 1000)+
  geom_hline(yintercept = seq(0, 500, 100), linewidth = I(1))+
  scale_x_continuous(limits=c(0,1000), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,500), expand=c(0,0)) +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.minor = element_blank())
