source('function.R')


##########QNN VS other
xr = 1000
yr = 500
N = 50000
r = (-(xr+yr)-sqrt((xr+yr)^2-4*(1-N)*xr*yr))/(2*(1-N))
set.seed(1234)
x = runif(N, 0, 1000)
y = runif(N, 0, 500)

df_ra = data.frame(x = x, y = y)
ppp_ra = as.ppp(df_ra, c(0, 1000, 0, 500))
meuse <- sf::st_as_sf(df_ra, coords = c("x", "y"),  
                      crs = 28992, agr = "constant")
qua_dat = quadratcount(ppp_ra, 1000/20, 500/20)
nn1d_time = timing(QNN_C(df_ra, xr, yr, r)$ar)



df <- data.frame(index = c('CCV-Shannon', 'CCV-Gini-Simpson', 'SDM', 'CD', 'Morisita'), time = c(timing(ccv(xr, yr, df_ra, 'shannon'))/nn1d_time,
                                                                                                 timing(ccv(xr, yr, df_ra, 'gini-simpson'))/nn1d_time,
                                                                                                 timing(fit(qua_dat)$par)/nn1d_time,
                                                                                                 timing(DC(qua_dat)$dc)/nn1d_time,
                                                                                                 timing(dispindmorisita(qua_dat)$imor)/nn1d_time))


#x11(20, 10)
p<-ggplot(data=df, aes(x=factor(index, levels = unique(index)), y=time)) +
  geom_bar(stat="identity", width=.5)+
  xlab('')+
  ylab('運算時間比值')+
  geom_text(aes(label = round(time, 2)), vjust = -0.5, colour = "red", size=9)+
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 25), 
        axis.text.y = element_text(size = 25), 
        legend.title = element_text(size = 25), 
        legend.text = element_text(size = 20), 
        legend.key.height= unit(.8, 'cm'),
        legend.key.width= unit(2, 'cm'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 25),
        plot.margin = margin(.4, .5, , .5, "cm"))
p



######################QNN VS NN
xr = 1000
yr = 500
N = 50000
r = (-(xr+yr)-sqrt((xr+yr)^2-4*(1-N)*xr*yr))/(2*(1-N))
df = data.frame(x = runif(N, 0, xr), y = runif(N, 0, yr))
meuse <- sf::st_as_sf(df, coords = c("x", "y"),  
                      crs = 28992, agr = "constant")

qnn_c = timing(QNN_C(df, 1000, 500, r)$ar)
qnn_r = timing(QNN(df, 1000, 500, r)$ar)
nn_c = timing(nni(meuse)$NNI)
nn_r = timing(CECI(df, 1000*500, method = "NP")$R)
nn_r_p = timing(CECI(df, 1000*500, method = "P")$R)

time = data.frame(qnn_c = qnn_c, qnn_r = qnn_r, nn_c = nn_c, nn_r = nn_r, nn_r_p = nn_r_p)
time