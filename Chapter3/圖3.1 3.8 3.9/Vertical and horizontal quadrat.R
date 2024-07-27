library(ggplot2)
set.seed(1003)
x = runif(20, 0, 1000)
y = runif(20, 0, 500)
df = data.frame(x = x, y = y, ycut = as.numeric(as.character(cut(y, breaks = seq(0, 500, 100), labels = seq(50, 450, 100)))), xcut = as.numeric(as.character(cut(x, breaks = seq(0, 1000, 100), labels = seq(50, 950, 100)))))


######Vertical
ggplot(df, aes(x = x, y = y))+
  geom_point(size = I(4))+
  ylim(0, 500)+
  xlim(0, 1000)+
  geom_vline(xintercept = seq(0, 1000, 100), linewidth = I(2))+
  geom_vline(xintercept = seq(50, 950, 100), linewidth = I(1), col = 'firebrick1')+
  geom_point(data = df, aes(x = xcut, y = y), col = 'firebrick1', size = I(4))+
  scale_x_continuous(limits=c(-1,1001), expand=c(0,0)) +
  scale_y_continuous(limits=c(-1,501), expand=c(0,0)) +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.minor = element_blank())


######Horizontal
ggplot(df)+
  geom_point(aes(x = x, y = y), size = I(4))+
  ylim(0, 500)+
  xlim(0, 1000)+
  geom_hline(yintercept = seq(0, 500, 100), linewidth = I(2))+
  geom_hline(yintercept = seq(50, 450, 100), linewidth = I(1), col = 'firebrick1')+
  geom_point(data = df, aes(x = x, y = ycut), col = 'firebrick1', size = I(4))+
  scale_x_continuous(limits=c(-1,1001), expand=c(0,0)) +
  scale_y_continuous(limits=c(-1,501), expand=c(0,0)) +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.minor = element_blank())

######Vertical and Horizontal
ggplot(df, aes(x = x, y = y))+
  geom_point(size = I(4))+
  ylim(0, 500)+
  xlim(0, 1000)+
  geom_vline(xintercept = seq(0, 1000, 100), linewidth = I(2), col = 'red')+
  geom_hline(yintercept = seq(0, 500, 100), linewidth = I(2), col = 'blue')+
  scale_x_continuous(limits=c(-1,1001), expand=c(0,0)) +
  scale_y_continuous(limits=c(-1,501), expand=c(0,0)) +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.minor = element_blank())

