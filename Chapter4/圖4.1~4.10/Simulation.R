
source('function.R')


################ Regular

xr = 1000
yr = 500

r = 1000/9
h = 500/9

x = rep(seq(0, xr , r), 10)

y = t = rep(yr, xr/r+1)
for (i in 1:9){
  t = t-h
  y = c(y, t)
}

df_re <- data.frame(x = x, y = y)
df_re <- df_re[which(df_re$y >= 0),]

#x11(20, 10)
ggplot(data = df_re, aes(x = x, y = y))+
  geom_point(size = I(3))+
  xlim(0, 1000)+
  scale_x_continuous(limits=c(0,1000), expand=c(0,4.5)) +
  scale_y_continuous(limits=c(0,500), expand=c(0,4.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 25), 
        axis.text.y = element_text(size = 25), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(.4, .8, , , "cm"))

################ Random
xr = 1000
yr = 500
N = 100
set.seed(1234)

df_ra = data.frame(x = runif(N, 0, 1000), y = runif(N, 0, 500))
#x11(20, 10)
ggplot(data = df_ra, aes(x = x, y = y))+
  geom_point(size = I(3))+
  xlim(0, 1000)+
  scale_x_continuous(limits=c(0,1000), expand=c(0,4.5)) +
  scale_y_continuous(limits=c(0,500), expand=c(0,4.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 25), 
        axis.text.y = element_text(size = 25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(.4, .8, , , "cm"))

################ Aggregation-1
xr = 1000
yr = 500
N = 100
set.seed(1234)

df_ag1 = species.distribution(delta = 50, lambda = 10, prob = 0.01, cnum = xr, rnum = yr)
df_ag1 <- df_ag1[sample(1:dim(df_ag1)[1], N),] %>% as.data.frame()
colnames(df_ag1) <- c('x', 'y')
#x11(20, 10)
ggplot(data = df_ag1, aes(x = x, y = y))+
  geom_point(size = I(3))+
  xlim(0, 1000)+
  scale_x_continuous(limits=c(0,1000), expand=c(0,4.5)) +
  scale_y_continuous(limits=c(0,500), expand=c(0,4.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 25), 
        axis.text.y = element_text(size = 25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(.4, .8, , , "cm"))

################ Aggregation-2
xr = 1000
yr = 500
N = 100
set.seed(1234)

df_ag2 = species.distribution(delta = 500, lambda = 10, prob = 0.01, cnum = xr, rnum = yr)
df_ag2 <- df_ag2[sample(1:dim(df_ag2)[1], N),] %>% as.data.frame()
colnames(df_ag2) <- c('x', 'y')

#x11(20, 10)

ggplot(data = df_ag2, aes(x = x, y = y))+
  geom_point(size = I(3))+
  xlim(0, 1000)+
  scale_x_continuous(limits=c(0,1000), expand=c(0,4.5)) +
  scale_y_continuous(limits=c(0,500), expand=c(0,4.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 25), 
        axis.text.y = element_text(size = 25), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(.4, .8, , , "cm"))

################ Aggregation-3
xr = 1000
yr = 500
N = 100
set.seed(1234)

df_ag3 = species.distribution(delta = 5000, lambda = 10, prob = 0.01, cnum = xr, rnum = yr)
df_ag3 <- df_ag3[sample(1:dim(df_ag3)[1], N),] %>% as.data.frame()
colnames(df_ag3) <- c('x', 'y')
#x11(20, 10)

ggplot(data = df_ag3, aes(x = x, y = y))+
  geom_point(size = I(3))+
  xlim(0, 1000)+
  scale_x_continuous(limits=c(0,1000), expand=c(0,4.5)) +
  scale_y_continuous(limits=c(0,500), expand=c(0,4.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 25), 
        axis.text.y = element_text(size = 25), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(.4, .8, , , "cm"))


################ QNN
df <- data.frame(pattern = c('regular', 'random', 'aggregation 1', 'aggregation 2', 'aggregation 3'), index = c(QNN(df_re, xr, yr, floor(1000/19), floor(500/9))$ar, 
                                                                                                                QNN(df_ra, xr, yr, r)$ar,
                                                                                                                QNN(df_ag1, xr, yr, r)$ar,
                                                                                                                QNN(df_ag2, xr, yr, r)$ar,
                                                                                                                QNN(df_ag3, xr, yr, r)$ar))

ggplot(data=df, aes(x=factor(pattern, levels = unique(pattern)), y=index)) +
  geom_bar(stat="identity", width=.5)+
  xlab('')+
  ylab('')+
  geom_text(aes(label = round(index, 2)), vjust = -0.35, colour = "red", size=9)+
  theme(plot.title = element_text(hjust = 0.5, size = 30), axis.text.x=element_text(size=30), axis.text.y=element_text(size=30))

################ NN
df <- data.frame(pattern = c('regular', 'random', 'aggregation 1', 'aggregation 2', 'aggregation 3'), index = c(CECI(df_re, xr*yr)$R, 
                                                                                                                CECI(df_ra, xr*yr)$R,
                                                                                                                CECI(df_ag1, xr*yr)$R,
                                                                                                                CECI(df_ag2, xr*yr)$R,
                                                                                                                CECI(df_ag3, xr*yr)$R))
#x11(20, 10)
p<-ggplot(data=df, aes(x=factor(pattern, levels = unique(pattern)), y=index)) +
  geom_bar(stat="identity", width=.5)+
  xlab('')+
  ylab('')+
  geom_text(aes(label = round(index, 2)), vjust = -0.35, colour = "red", size=9)+
  theme(plot.title = element_text(hjust = 0.5, size = 30), axis.text.x=element_text(size=30), axis.text.y=element_text(size=30))
p

############ CCV
df_sh <- data.frame(Pattern = c('regular', 'random', 'aggregation 1', 'aggregation 2', 'aggregation 3'), index = c(ccv(xr, yr, df_re, 'shannon'), 
                                                                                                                   ccv(xr, yr, df_ra, 'shannon'),
                                                                                                                   ccv(xr, yr, df_ag1, 'shannon'),
                                                                                                                   ccv(xr, yr, df_ag2, 'shannon'),
                                                                                                                   ccv(xr, yr, df_ag3, 'shannon')), Type = rep('shannon', 5))


df_gi <- data.frame(Pattern = c('regular', 'random', 'aggregation 1', 'aggregation 2', 'aggregation 3'), index = c(ccv(xr, yr, df_re, 'gini-simpson'), 
                                                                                                                   ccv(xr, yr, df_ra, 'gini-simpson'),
                                                                                                                   ccv(xr, yr, df_ag1, 'gini-simpson'),
                                                                                                                   ccv(xr, yr, df_ag2, 'gini-simpson'),
                                                                                                                   ccv(xr, yr, df_ag3, 'gini-simpson')), Type = rep('gini-simpson', 5))




df_ccv <- rbind(df_sh, df_gi)
#x11(20, 10)
p <- ggplot(data = df_ccv, aes(x = factor(Pattern, levels = unique(Pattern)), y = index, fill = Type, pattern = Type)) +
  geom_bar_pattern(stat='identity', position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.05,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.8)  +
  scale_pattern_manual(values = c("gini-simpson" = "stripe", "shannon" = "none")) +
  xlab('') +
  ylab('') +
  geom_text(aes(label = round(index, 5)), vjust = -0.5, colour = "red", size = 5, position = position_dodge(width = .9)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 30), 
        axis.text.y = element_text(size = 30), 
        legend.title = element_text(size = 30), 
        legend.text = element_text(size = 25), 
        legend.key.height= unit(.8, 'cm'),
        legend.key.width= unit(2, 'cm'))
p

################ SDM alpha and NBD K
set.seed(1003)
qua_re <- df_re %>% as.ppp(c(0, 1000, 0, 500)) %>% quadratcount(1000/20, 500/20)
qua_ra <- df_ra %>% as.ppp(c(0, 1000, 0, 500)) %>% quadratcount(1000/20, 500/20)
qua_ag1 <- df_ag1 %>% as.ppp(c(0, 1000, 0, 500)) %>% quadratcount(1000/20, 500/20)
qua_ag2 <- df_ag2 %>% as.ppp(c(0, 1000, 0, 500)) %>% quadratcount(1000/20, 500/20)
qua_ag3 <- df_ag3 %>% as.ppp(c(0, 1000, 0, 500)) %>% quadratcount(1000/20, 500/20)

df1_alpha <- data.frame(pattern = c('regular', 'random'), index = c(fit(qua_re)$par, 
                                                                    fit(qua_ra )$par), type = rep('alpha', 2))

df2_alpha <- data.frame(pattern = c('aggregation 1', 'aggregation 2', 'aggregation 3'), index = c(fit(qua_ag1)$par,
                                                                                                  fit(qua_ag2)$par,
                                                                                                  fit(qua_ag3)$par), type = rep('alpha', 3))


df1_k <- data.frame(pattern = c('regular', 'random'), index = c(fitNBD(qua_re)$par[1],
                                                                fitNBD(qua_ra)$par[1]), type = rep('k', 2))

df2_k <- data.frame(pattern = c('aggregation 1', 'aggregation 2', 'aggregation 3'), index = c(fitNBD(qua_ag1)$par[1],
                                                                                              fitNBD(qua_ag2)$par[1],
                                                                                              fitNBD(qua_ag3)$par[1]), type = rep('k', 3))



df1 <- rbind(df1_alpha, df1_k)
df2 <- rbind(df2_alpha, df2_k)

#x11(20, 10)
p1<-ggplot(data=df1, aes(x = factor(pattern, levels = unique(pattern)), y = index, fill = type, pattern = type)) +
  geom_bar_pattern(stat='identity', position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.05,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.8)  +
  scale_pattern_manual(values = c('alpha' = "stripe", "k" = "none")) +
  
  xlab('') +
  ylab('') +
  geom_text(aes(label = round(index, 1)), vjust = -0.5, colour = "red", size = 7, position = position_dodge(width = .9)) +
  
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 30), 
        axis.text.y = element_text(size = 30), 
        legend.title = element_text(size = 30), 
        legend.text = element_text(size = 25), 
        legend.key.height= unit(.8, 'cm'),
        legend.key.width= unit(2, 'cm'))
p1

p2<-ggplot(data=df2, aes(x = factor(pattern, levels = unique(pattern)), y = index, fill = type, pattern = type)) +
  geom_bar_pattern(stat='identity', position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.05,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.8)  +
  scale_pattern_manual(values = c('alpha' = "stripe", "k" = "none")) +
  xlab('') +
  ylab('') +
  geom_text(aes(label = round(index, 5)), vjust = -0.5, colour = "red", size = 7, position = position_dodge(width = .9)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), 
        axis.text.x = element_text(size = 30), 
        axis.text.y = element_text(size = 30), 
        legend.title = element_text(size = 30), 
        legend.text = element_text(size = 25), 
        legend.key.height= unit(.8, 'cm'),
        legend.key.width= unit(2, 'cm'))
p2

