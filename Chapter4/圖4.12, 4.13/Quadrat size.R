source('function.R')


################## Random
set.seed(1003)
xr = 1000
yr = 500
N = 500
r = (-(xr+yr)-sqrt((xr+yr)^2-4*(1-N)*xr*yr))/(2*(1-N))
size = c(5, 10, 20, round(r,2))
df <- data.frame()
for (i in size){
  for (j in 1:100){
    df_ra = data.frame(x = runif(N, 0, xr), y = runif(N, 0, yr))
    df <- rbind(df, c(i, 'NNQ', NNQ(df_ra, xr, yr, i)$R))
    df <- rbind(df, c(i, '1dNN', QNN(df_ra, xr, yr, i)$ar))
  }
}

names(df) <- c("size", "pattern", "index")


df[which(df$pattern == '1dNN'), 2] = '區塊最近鄰指標'
df[which(df$pattern == 'NNQ'), 2] = '最近鄰指標'

p2 <- ggplot(df, aes(x = factor(size, level = c('5', '10', '20', round(r,2))), y = as.numeric(index), fill = factor(pattern, levels = unique(pattern)))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = 'point', position = position_dodge(width = 0.75), shape = 20, size = 4, color = 'red') +
  
  geom_line(stat = 'summary', fun = 'mean', aes(group = factor(pattern, levels = unique(pattern)), color = factor(pattern, levels = unique(pattern))), position = position_dodge(width = 0.75), linetype = 'dashed') +
  scale_x_discrete(name = "帶狀區塊帶寬大小")+
  scale_y_continuous(name = "指標", breaks = seq(0, max(df$index), 0.1)) +
  #scale_fill_discrete(name = "指數類別")+
  labs(color = "平均數連線", fill = "指標類別")+
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 25), 
        axis.text.y = element_text(size = 25), 
        legend.title = element_text(size = 25), 
        legend.text = element_text(size = 20), 
        legend.key.height= unit(.8, 'cm'),
        legend.key.width= unit(2, 'cm'),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size =20))


#x11(15, 10)
p2



####################### Aggregation
set.seed(1003)
xr = 1000
yr = 500
N = 500
r = (-(xr+yr)-sqrt((xr+yr)^2-4*(1-N)*xr*yr))/(2*(1-N))
size = c(5, 10, 20, round(r,2))
df <- data.frame()
for (i in size){
  for (j in 1:100){
    test <- species.distribution(delta = 500, lambda = 100, prob = 0.01, cnum = xr, rnum = yr)
    test <- as.data.frame(test[sample(1:dim(test)[1], N),])
    df <- rbind(df, c(i, 'NNQ', NNQ(test, xr, yr, i)$R))
    df <- rbind(df, c(i, '1dNN', QNN(test, xr, yr, i)$ar))
  }
}

names(df) <- c("size", "pattern", "index")


df[which(df$pattern == '1dNN'), 2] = '區塊最近鄰指標'
df[which(df$pattern == 'NNQ'), 2] = '最近鄰指標'

p2 <- ggplot(df, aes(x = factor(size, level = c('5', '10', '20', round(r,2))), y = as.numeric(index), fill = factor(pattern, levels = unique(pattern)))) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = 'point', position = position_dodge(width = 0.75), shape = 20, size = 4, color = 'red') +
  
  geom_line(stat = 'summary', fun = 'mean', aes(group = factor(pattern, levels = unique(pattern)), color = factor(pattern, levels = unique(pattern))), position = position_dodge(width = 0.75), linetype = 'dashed') +
  scale_x_discrete(name = "帶狀區塊帶寬大小")+
  scale_y_continuous(name = "指標", breaks = seq(0, max(df$index), 0.1)) +
  #scale_fill_discrete(name = "指數類別")+
  labs(color = "平均數連線", fill = "指標類別")+
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 25), 
        axis.text.y = element_text(size = 25), 
        legend.title = element_text(size = 25), 
        legend.text = element_text(size = 20), 
        legend.key.height= unit(.8, 'cm'),
        legend.key.width= unit(2, 'cm'),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size =20))


#x11(15, 10)
p2




