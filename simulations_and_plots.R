# VG simulation function
VG=function(sigma, nu, mu, T, N) {
  a=1/nu
  b=1/nu
  h=T/N
  t=(0:T)/N
  X=rep(0, N+1)
  I=rep(0,N)
  X[1]=0
  for(i in 1:N) {
    I[i]=rgamma(1,a*h,b)
    X[i+1]=X[i] + mu*I[i]+sigma*sqrt(I[i])*rnorm(1)
  }
  return(X)
}


VG_increments=function(sigma, nu, mu, T, N) {
  a=1/nu
  b=1/nu
  h=T/N
  t=(0:T)/N
  inc=rep(0, N)
  I=rep(0,N)
  for(i in 1:N) {
    I[i]=rgamma(1,a*h,b)
    inc[i] <- mu*I[i]+sigma*sqrt(I[i])*rnorm(1)
  }
  return(inc)
}


BM_increments <- function(sigma, T, N, theta =0){
  drift <- theta*seq(from = 0, to=T, length.out = N)
  B_inc <- sigma*rnorm(N,0,1)*sqrt((T/N)) + drift
  return(B_inc)
}

# ---------------------------
# downloading data

#install.packages("quantmod") # if needed
library("quantmod")
ftse100 <- new.env()
getSymbols("^FTSE", env = ftse100, src = "yahoo",
           from = as.Date("1984-03-01"), to = as.Date("2017-01-30"))

FTSE <- ftse100$FTSE
head(FTSE)

r = diff(log(FTSE$FTSE.Close))
r = r[2:length(r)]
names(r)=c("FTSE logReturn")
head(r)
plot(as.vector(r), type= 'line')

rownames(r)

x_observations <- as.vector(r) - mean(as.vector(r))
mean(x_observations)
var(x_observations)

# estimation of parameters:
sigma_sq_hat  <- sum(x_observations^2)/length(x_observations)
nu_hat <- (sum(x_observations^4)/length(x_observations))*(1/(3*sigma_sq_hat^2)) -1
sigma_hat <- sqrt(sigma_sq_hat)

# ---------------------------------------------
# plot comparing to normal
plot(density(x_observations))
lines(density(rnorm(length(x_observations), m=0, sd = sd(x_observations))), col='red', add=TRUE)
lines(density(x_vg_inc), col='green')

#----------------------------------------------
# generating from the VG distribution
N <- 3*length(x_observations)
x_vg_inc <- VG_increments(sigma_hat, nu=nu_hat, mu=0, T=length(x_observations), N)


length(x_vg_inc)
length(x_observations)

# plot of increments
plot(1:length(x_observations), x_observations)
lines(seq(from =1, to=length(x_observations), length.out = N), x_vg_inc, col = 'blue')

# plot of cumsum
plot(1:length(x_observations), cumsum(x_observations), ylim = c(-2,2))
lines(seq(from =1, to=length(x_observations), length.out = N), cumsum(x_vg_inc), col = 'blue')

#-----------------------------------
# comparison to the Brownian Motion

bm_simulation <- BM_increments(sigma=sigma_hat, T=length(x_observations), N=N, theta =0)

plot(1:length(x_observations), cumsum(x_observations), ylim = c(-2,2))
lines(seq(from =1, to=length(x_observations), length.out = N), cumsum(bm_simulation), col = 'blue')


#------------------------------------
# comparing MSE of BM and VG
MSE <- function(bm, vg, real_data){
  mse_bm <- (sum(bm - real_data)^2)/length(real_data)
  mse_vg <- (sum(vg - real_data)^2)/length(real_data)
  return(c(mse_bm, mse_vg))
}

#----------------------------------
# moments
#install.packages('moments')
library(moments)

skewness(rnorm(10000000,0,1))
kurtosis(rnorm(10000000,0,100))


#---------------------------------
# full comparison
bm_simulation <- BM_increments(sigma=sigma_hat, T=length(x_observations), N=length(x_observations), theta =0)
vg_simulation <- VG_increments(sigma = sigma_hat, nu = nu_hat, mu=0, T=length(x_observations), N=length(x_observations))

MSE(bm_simulation, vg_simulation, x_observations)
skewness(x_observations)
kurtosis(x_observations)

kurtosis(bm_simulation)
skewness(bm_simulation)

kurtosis(vg_simulation)
skewness(vg_simulation)

#---------------------------------
n_iter <- 1000
bm_data <- t(sapply(1:n_iter, function(k){
  bm_simulation <- BM_increments(sigma=sigma_hat, T=length(x_observations), N=length(x_observations), theta =0)
  return(c(kurtosis(bm_simulation), skewness(bm_simulation)))
}))

bm_data<- as.data.frame(bm_data)
colnames(bm_data) <- c("kurtosis", "skewness")

vg_data <- t(sapply(1:n_iter, function(k){
  vg_simulation <- VG_increments(sigma = sigma_hat, nu = nu_hat, mu=0, T=length(x_observations), N=length(x_observations))
  return(c(kurtosis(vg_simulation), skewness(vg_simulation)))
}))

vg_data<- as.data.frame(vg_data)
colnames(vg_data) <- c("kurtosis", "skewness")

comparison_data = rbind(bm_data, vg_data)
comparison_data$process <- rep(c('Brownian Motion', 'Variance-Gamma'), each = n_iter)


library(ggplot2)
p1 <- ggplot(vg_data, aes(x=kurtosis)) + geom_histogram(bins=50) +  
  geom_vline(xintercept = kurtosis(x_observations), col='red', size =1) +
  geom_vline(xintercept = 3, col='blue', linetype=2, size =1)+theme(legend.title = element_blank(),
                                                                    axis.title.x = element_text(size = 16),
                                                                    axis.title.y = element_text(size = 16),
                                                                    axis.text = element_text(size = 12))
p1
ggsave('kurtosis16.pdf', p1, height=6, width=8)


p2 <- ggplot(vg_data, aes(x=skewness)) + geom_histogram(bins=50) +  
  geom_vline(xintercept = skewness(x_observations), col='red', size =1) +
  geom_vline(xintercept = 0, col='blue', linetype=2, size =1) +theme(legend.title = element_blank(),
                                                                     axis.title.x = element_text(size = 16),
                                                                     axis.title.y = element_text(size = 16),
                                                                     axis.text = element_text(size = 12))
p2
ggsave('skewness16.pdf', p2, height=6, width=8)


#--------------------------
#density plot
density_data <- data.frame(value = c(x_observations,
                                    rnorm(length(x_observations), m=0, sd = sd(x_observations)),
                                    VG_increments(sigma = sigma_hat, nu = nu_hat, mu=0, T=length(x_observations), N=length(x_observations))),
                           type = rep(c('observed values', 'Gaussian', 'simulated Variance-Gamma'), each = length(x_observations))
                           )
p3 <- ggplot(density_data, aes(x=value, colour=type)) + geom_density(size =1) + xlim(c(-0.1, 0.1)) +
  scale_color_manual(values = c('blue', 'black', 'red')) +theme(legend.title = element_blank(),
                                                                axis.title.x = element_text(size = 16),
                                                                axis.title.y = element_text(size = 16),
                                                                axis.text = element_text(size = 12),
                                                                legend.text= element_text(size=16))
p3
ggsave('density16.pdf', p3, height=6, width=8)

#--------------------------
# path plots
d1 <- data.frame(value = cumsum(x_observations))
d1$Date <- index(FTSE)[-1]
d1$path <- c('observations')

d2 <- data.frame(value = c(cumsum(VG_increments(sigma = sigma_hat, nu = nu_hat, mu=0, T=length(x_observations), N=length(x_observations)/4)),
                           cumsum(VG_increments(sigma = sigma_hat, nu = nu_hat, mu=0, T=length(x_observations), N=length(x_observations)/4)),
                           cumsum(VG_increments(sigma = sigma_hat, nu = nu_hat, mu=0, T=length(x_observations), N=length(x_observations)/4)),
                           cumsum(VG_increments(sigma = sigma_hat, nu = nu_hat, mu=0, T=length(x_observations), N=length(x_observations)/4)),
                           cumsum(VG_increments(sigma = sigma_hat, nu = nu_hat, mu=0, T=length(x_observations), N=length(x_observations)/4)),
                           cumsum(VG_increments(sigma = sigma_hat, nu = nu_hat, mu=0, T=length(x_observations), N=length(x_observations)/4))))
                     
d2$path <- as.factor(rep(1:(nrow(d2)/length(x_observations)/4), each = length(x_observations)/4))
d2$Date <- rep(index(FTSE)[-1], times = (1/4)*(nrow(d2)/length(x_observations)))

#scale_x_date(format = "%b-%Y")

p4 <- ggplot(d2, aes(x= Date, y = value)) + 
  geom_segment(aes(xend = Date+1, yend = value, col=path), size =0.75) +  
  geom_line(data = d1, aes(Date, value, col = path), col = 'black') +  theme(legend.position = "none",
                                                                             axis.title.x = element_text(size = 16),
                                                                           axis.title.y = element_text(size = 16),
                                                                           axis.text = element_text(size = 12)) +
  scale_colour_manual(values = c('orange', 'red', 'green', 'blue', 'pink', 'brown'))

p4
ggsave('path16.pdf',p4, height=6, width=10)  

#--------------------------
# path plots
# MSE test

n_iter <- 1000
mse_data <- t(sapply(1:n_iter, function(k){
  bm_simulation <- BM_increments(sigma=sigma_hat, T=length(x_observations), N=length(x_observations), theta =0)
  vg_simulation <- VG_increments(sigma = sigma_hat, nu = nu_hat, mu=0, T=length(x_observations), N=length(x_observations))
  
  mse <- MSE(cumsum(bm_simulation), cumsum(vg_simulation), cumsum(x_observations))
  ks_bm <- ks.test(x_observations, bm_simulation)$p.value
  ks_vg <- ks.test(x_observations, vg_simulation)$p.value
  return(c(mse, ks_bm, ks_vg))
}))



mse_data <- as.data.frame(mse_data)
colnames(mse_data) <- c('mse_bm', 'mse_vg', 'ks_bm', 'ks_vg')
head(mse_data)
class(mse_data)
is.data.frame(mse_data)
hist(mse_data$ks_vg, 50)


shapiro.test(sample(size = 5000,x_observations))

hist(mse_data$mse_bm)
d3 <- data.frame(mse = c(mse_data$mse_bm, mse_data$mse_vg), process = rep(c('Brownian Motion', 'Variance-Gamma')))

p5 <- ggplot(d3, aes(x=mse, colour=process)) + geom_density(size =1) + 
  scale_color_manual(values = c('blue', 'red')) +theme(legend.title = element_blank(),
                                                       axis.title.x = element_text(size = 16),
                                                       axis.title.y = element_text(size = 16),
                                                       legend.text = element_text(size = 16),
                                                       axis.text = element_text(size = 12))
p5
ggsave('mse16.pdf',p5, height=6, width=10)  


