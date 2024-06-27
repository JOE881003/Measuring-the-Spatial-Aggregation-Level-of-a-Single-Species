source('./function.R')

################ Random(guarsp)
sp = unique(BCI_List[[8]]$sp)
s = 30
BCI_data = BCI_List[[8]][which(BCI_List[[8]]$status=="A" & BCI_List[[8]]$dbh >= 10 & BCI_List[[8]]$sp == sp[s]), 5:6]
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
BCI_data = BCI_List[[8]][which(BCI_List[[8]]$status=="A" & BCI_List[[8]]$dbh >= 10 & BCI_List[[8]]$sp == sp[30]), 5:6]
N = dim(BCI_data)[1]
r = (-(xr+yr)-sqrt((xr+yr)^2-4*(1-N)*xr*yr))/(2*(1-N))
t = QNN(BCI_data, x, y, r)$c
result = data.frame(QNN = QNN(BCI_data, xr, yr, r)$ar, pvalue = 2*stats::pnorm(-abs(t)))
result


############# Aggregation(crotbi)

xr = 1000
yr = 500
s = 8
sp = unique(BCI_List[[8]]$sp)
BCI_data = BCI_List[[8]][which(BCI_List[[8]]$status=="A" & BCI_List[[8]]$dbh >= 10 & BCI_List[[8]]$sp == sp[s]), 5:6]
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

############## Random or aggregation(psycho)
xr = 1000
yr = 500
s = 200
sp = unique(BCI_List[[8]]$sp)
BCI_data = BCI_List[[8]][which(BCI_List[[8]]$status=="A" & BCI_List[[8]]$dbh >= 10 & BCI_List[[8]]$sp == sp[s]), 5:6]
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


####################### elevation
BCI_data = BCI_List[[8]][which(BCI_List[[8]]$status=="A" & BCI_List[[8]]$dbh >= 10), c(3, 5, 6)]
sp = unique(BCI_List[[8]]$sp)
elve_data <- read.csv("./BCIelev.csv")
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


NNL(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation, 39.67)

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


NNL(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation, 39.67)

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


NNL(BCI_data[which(BCI_data$sp == sp[s]), ]$elevation, 39.67)

