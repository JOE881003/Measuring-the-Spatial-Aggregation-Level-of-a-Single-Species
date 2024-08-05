source('function.R')


################ Random(guarsp)
sp = unique(bci.tree8$sp)
s = 30
BCI_data = bci.tree8[which(bci.tree8$status=="A" & bci.tree8$dbh >= 10 & bci.tree8$sp == sp[s]), 5:6]
#x11(20, 10)

ggplot(data = BCI_data, aes(x = gx, y = gy))+
  geom_point(size = I(2))+
  xlim(0, 1000)+
  scale_x_continuous(limits=c(0,1000), expand=c(0,4.5)) +
  scale_y_continuous(limits=c(0,500), expand=c(0,4.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = margin(.4, .7, , , "cm"))


xr = 1000
yr = 500
BCI_data = bci.tree8[which(bci.tree8$status=="A" & bci.tree8$dbh >= 10 & bci.tree8$sp == sp[30]), 5:6]
N = dim(BCI_data)[1]
r = (-(xr+yr)-sqrt((xr+yr)^2-4*(1-N)*xr*yr))/(2*(1-N))
t = QNN(BCI_data, xr, yr, r)$c



result = data.frame(QNN = QNN(BCI_data, xr, yr, r)$ar, pvalue = 2*stats::pnorm(-abs(t)))
result


result2 = data.frame(index = unlist(lapply(c(1:floor(round(r, 2)), round(r, 2)), function(x) QNN(BCI_data, xr, yr, x)$ar)), size = c(as.character(1:floor(round(r, 2))), as.character(round(r, 2))))
x11(20, 10)
ggplot(data = result2, aes(x = factor(size, levels = unique(size)), y = index, group = 1))+
  geom_line(linewidth = 1)+
  geom_point(size = 4)+
  ylim(0, 1.1)+
  xlab('帶寬(公尺)')+
  ylab('QNN')+
  #annotate('text', x = result2$size, y = result2$index-0.05, label = round(result2$index, 2), col = 'red', size = 7)+
  scale_y_continuous(limits=c(0,1.1), expand=c(0,0)) +
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x=element_text(size = 20),
        axis.title.y=element_text(size = 20),
        plot.margin = margin(.4, .4, .2, .2, "cm"))
############# Aggregation(crotbi)

xr = 1000
yr = 500
s = 8
sp = unique(bci.tree8$sp)
BCI_data = bci.tree8[which(bci.tree8$status=="A" & bci.tree8$dbh >= 10 & bci.tree8$sp == sp[s]), 5:6]
N = dim(BCI_data)[1]
#x11(20, 10)

ggplot(data = BCI_data, aes(x = gx, y = gy))+
  geom_point(size = I(2))+
  xlim(0, 1000)+
  scale_x_continuous(limits=c(0,1000), expand=c(0,4.5)) +
  scale_y_continuous(limits=c(0,500), expand=c(0,4.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = margin(.4, .7, , , "cm"))




r = (-(xr+yr)-sqrt((xr+yr)^2-4*(1-N)*xr*yr))/(2*(1-N))
t = QNN(BCI_data, xr, yr, r)$c
result = data.frame(QNN = QNN(BCI_data, xr, yr, r)$ar, pvalue = 2*stats::pnorm(-abs(t)))
result

result2 = data.frame(index = unlist(lapply(c(1:floor(round(r, 2)), round(r, 2)), function(x) QNN(BCI_data, xr, yr, x)$ar)), size = c(as.character(1:floor(round(r, 2))), as.character(round(r, 2))))
x11(20, 10)

ggplot(data = result2, aes(x = factor(size, levels = unique(size)), y = index, group = 1))+
  geom_line(linewidth = 1)+
  geom_point(size = 4)+
  ylim(0, 1.1)+
  xlab('帶寬(公尺)')+
  ylab('QNN')+
  #annotate('text', x = result2$size, y = result2$index-0.05, label = round(result2$index, 2), col = 'red', size = 7)+
  scale_y_continuous(limits=c(0,1.1), expand=c(0,0)) +
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x=element_text(size = 20),
        axis.title.y=element_text(size = 20),
        plot.margin = margin(.4, .4, .2, .2, "cm"))


ggplot(data = result2, aes(x = factor(size, levels = unique(size)), y = index, group = 1))+
  geom_line(linewidth = 1)+
  geom_point(size = 4)+
  ylim(0, 1)+
  xlab('帶寬(公尺)')+
  ylab('QNN')+
  annotate('text', x = 1, y = result2$index[1]-0.03, label = round(result2$index[1], 2), col = 'red', size = 7)+
  annotate('text', x = 2, y = result2$index[2]-0.03, label = round(result2$index[2], 2), col = 'red', size = 7)+
  annotate('text', x = 3, y = result2$index[3]-0.03, label = round(result2$index[3], 2), col = 'red', size = 7)+
  annotate('text', x = 4, y = result2$index[4]-0.03, label = round(result2$index[4], 2), col = 'red', size = 7)+
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x=element_text(size = 20),
        axis.title.y=element_text(size = 20),
        plot.margin = margin(.4, .4, .2, .2, "cm"))
############## Random or aggregation(psycho)
xr = 1000
yr = 500
s = 200
sp = unique(bci.tree8$sp)
BCI_data = bci.tree8[which(bci.tree8$status=="A" & bci.tree8$dbh >= 10 & bci.tree8$sp == sp[s]), 5:6]
N = dim(BCI_data)[1]
#x11(20, 10)

ggplot(data = BCI_data, aes(x = gx, y = gy))+
  geom_point(size = I(2))+
  xlim(0, 1000)+
  scale_x_continuous(limits=c(0,1000), expand=c(0,4.5)) +
  scale_y_continuous(limits=c(0,500), expand=c(0,4.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = margin(.4, .7, , , "cm"))




r = (-(xr+yr)-sqrt((xr+yr)^2-4*(1-N)*xr*yr))/(2*(1-N))
t = QNN(BCI_data, xr, yr, r)$c
result = data.frame(QNN = QNN(BCI_data, xr, yr, r)$ar, pvalue = 2*stats::pnorm(-abs(t)))
result

result2 = data.frame(index = unlist(lapply(c(1:floor(round(r, 2)), round(r, 2)), function(x) QNN(BCI_data, xr, yr, x)$ar)), size = c(as.character(1:floor(round(r, 2))), as.character(round(r, 2))))
x11(20, 10)

ggplot(data = result2, aes(x = factor(size, levels = unique(size)), y = index, group = 1))+
  geom_line(linewidth = 1)+
  geom_point(size = 4)+
  ylim(0, 1.1)+
  xlab('帶寬(公尺)')+
  ylab('QNN')+
  #annotate('text', x = result2$size, y = result2$index-0.05, label = round(result2$index, 2), col = 'red', size = 7)+
  scale_y_continuous(limits=c(0,1.1), expand=c(0,0)) +
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x=element_text(size = 20),
        axis.title.y=element_text(size = 20),
        plot.margin = margin(.4, .4, .2, .2, "cm"))



ggplot(data = result2, aes(x = factor(size, levels = unique(size)), y = index, group = 1))+
  geom_line(linewidth = 1)+
  geom_point(size = 4)+
  ylim(0, 1)+
  xlab('帶寬(公尺)')+
  ylab('QNN')+
  annotate('text', x = 1, y = result2$index[1]-0.03, label = round(result2$index[1], 2), col = 'red', size = 7)+
  annotate('text', x = 2, y = result2$index[2]-0.03, label = round(result2$index[2], 2), col = 'red', size = 7)+
  annotate('text', x = 3, y = result2$index[3]-0.03, label = round(result2$index[3], 2), col = 'red', size = 7)+
  annotate('text', x = 4, y = result2$index[4]-0.03, label = round(result2$index[4], 2), col = 'red', size = 7)+
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x=element_text(size = 20),
        axis.title.y=element_text(size = 20),
        plot.margin = margin(.4, .4, .2, .2, "cm"))
####################### elevation
BCI_data = bci.tree8[which(bci.tree8$status=="A" & bci.tree8$dbh >= 10), c(3, 5, 6)]
sp = unique(bci.tree8$sp)

elve_data <- read.csv("BCIelev.csv")
elve_mat <- acast(elve_data, y ~ x, value.var = "elev")
BCI_data$elevation <- mapply(get_elevation, BCI_data$gx, BCI_data$gy, MoreArgs = list(elve_mat, elve_data))


######################poutre
s = 48
ggplot(elve_data, aes(x = x, y = y, fill = elev)) +
  geom_tile() +
  scale_fill_gradient(low = "#0000FF", high ="#FF0000",limits = ceiling(c(min(elve_data$elev), max(elve_data$elev))), breaks=pretty_breaks(n=4)(min(elve_data$elev):max(elve_data$elev))) +
  labs(title = "海拔高度熱力圖與樹的位置", x = "X 座標", y = "Y 座標", fill = "海拔高度") +
  theme_minimal() +
  geom_point(data = BCI_data[which(BCI_data$sp == sp[s]), ], aes(x = gx, y = gy, fill = elevation), color = "black", size = 2)+ # 添加樹的位置
  theme(plot.title = element_blank(), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        legend.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.key.height= unit(.8, 'cm'),
        legend.key.width= unit(1.5, 'cm'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))

L  = max(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation)-min(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation)
NNL(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation, L)$NNL
NNL(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation, L)$p

######################beilpe
s = 4
ggplot(elve_data, aes(x = x, y = y, fill = elev)) +
  geom_tile() +
  scale_fill_gradient(low = "#0000FF", high ="#FF0000",limits = ceiling(c(min(elve_data$elev), max(elve_data$elev))), breaks=pretty_breaks(n=4)(min(elve_data$elev):max(elve_data$elev))) +
  labs(title = "海拔高度熱力圖與樹的位置", x = "X 座標", y = "Y 座標", fill = "海拔高度") +
  theme_minimal() +
  geom_point(data = BCI_data[which(BCI_data$sp == sp[s]), ], aes(x = gx, y = gy, fill = elevation), color = "black", size = 2)+ # 添加樹的位置
  theme(plot.title = element_blank(), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        legend.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.key.height= unit(.8, 'cm'),
        legend.key.width= unit(1.5, 'cm'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))

L  = max(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation)-min(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation)
NNL(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation, L)$NNL
NNL(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation, L)$p

######################maquco
s = 115
ggplot(elve_data, aes(x = x, y = y, fill = elev)) +
  geom_tile() +
  scale_fill_gradient(low = "#0000FF", high ="#FF0000",limits = ceiling(c(min(elve_data$elev), max(elve_data$elev))), breaks=pretty_breaks(n=4)(min(elve_data$elev):max(elve_data$elev))) +
  labs(title = "海拔高度熱力圖與樹的位置", x = "X 座標", y = "Y 座標", fill = "海拔高度") +
  theme_minimal() +
  geom_point(data = BCI_data[which(BCI_data$sp == sp[s]), ], aes(x = gx, y = gy, fill = elevation), color = "black", size = 2)+ # 添加樹的位置
  theme(plot.title = element_blank(), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        legend.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.key.height= unit(.8, 'cm'),
        legend.key.width= unit(1.5, 'cm'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))


L  = max(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation)-min(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation)
NNL(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation, L)$NNL
NNL(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation, L)$p

