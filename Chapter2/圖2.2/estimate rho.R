source("function.R")

f <- function(N, type){
  xr = 1000
  yr = 500
  x = runif(N, 0, 1000)
  y = runif(N, 0, 500)
  df_ra = data.frame(x = x, y = y)
  if(type == 'est'){
    return(round((N-1)/sum(nndist(df_ra)**2*pi), 4))
  }else{
    return(round(N/(1000*500),4))
  }
}

df <- data.frame()
set.seed(47)
for (i in seq(1000, 50000, 1000)){
  
  df = rbind(df, c(f(i, 'est'), 'est', i))
  df = rbind(df, c(f(i, 'theor'), 'theor',i))
}

names(df) <- c("value", "type", "N")
df[which(df$type == 'est'), 2] = '不偏估計量'
df[which(df$type == 'theor'), 2] = '觀測密度'

x11(20, 10)
ggplot(df, aes(x = as.numeric(N), y = as.numeric(value), group = type, color = type))+
  geom_line(size = 1)+
  ylab(expression(hat(rho)))+
  xlab('個體數')+
  labs(color = '計算方式')+
  theme(plot.title = element_blank(), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        legend.title = element_text(size = 30), 
        legend.text = element_text(size = 20), 
        legend.key.height= unit(.8, 'cm'),
        legend.key.width= unit(1.5, 'cm'),
        axis.title.x=element_text(size = 20),
        axis.title.y=element_text(size = 20, angle = 0, vjust = 0.5))