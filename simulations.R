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

CGMY=function(sigma,nu,mu,Y,T,N)
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


CGMY=function(sigma,nu,mu,Y,T,N) {
  h=T/N
  t=(0:T)/N
  eta_p = sqrt(((mu^2*nu^2)/4)+((sigma^2*nu)/2)) + (mu*nu)/2
  eta_n = sqrt(((mu^2*nu^2)/4)+((sigma^2*nu)/2)) - (mu*nu)/2
  C = 1/nu
  G = 1/eta_n
  M = 1/eta_p
  A=(G-M)/2
  B=(G+M)/2
  f=function(y){exp(-(B^2-A^2)*y/2)*gamma(Y)*CGMY_f(y,Y,B)/(gamma(Y/2)*2^(Y/2-1))}
  X=rep(0, N+1)
  I=rep(0,N)
  probJ=rep(0,N)
  J=(0:(N-1))*h + h
  for (j in 1:N) {probJ[[j]]<-f(J[[j]])}
  randf=sample(J,N,probJ,replace=TRUE)
  X[1]=0
  for(i in 1:N) {
    I[i]=randf[[i]]
    X[i+1]=X[i] + A*I[i]+sqrt(I[i])*rnorm(1)
  }
  return(X)
}

CGMY_increments=function(sigma,nu,mu,Y,T,N) {
  h=T/N
  t=(0:T)/N
  eta_p = sqrt(((mu^2*nu^2)/4)+((sigma^2*nu)/2)) + (mu*nu)/2
  eta_n = sqrt(((mu^2*nu^2)/4)+((sigma^2*nu)/2)) - (mu*nu)/2
  C = 1/nu
  G = 1/eta_n
  M = 1/eta_p
  A=(G-M)/2
  B=(G+M)/2
  f=function(y){exp(-(B^2-A^2)*y/2)*cgamma(Y)*CGMY_f(y,Y,B)/(cgamma(Y/2)*2^(Y/2-1))}
  inc=rep(0, N)
  I=rep(0,N)
  probJ=rep(0,N)
  J=(0:(N-1))*h + h
  for (j in 1:N) {probJ[[j]]<-f(J[[j]])}
  randf=sample(J,N,probJ,replace=TRUE)
  X[1]=0
  for(i in 1:N) {
    I[i]=randf[[i]]
    inc[i]= A*I[i]+sqrt(I[i])*rnorm(1)
  }
  return(X)
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

#-----------------------------------
# comparison to the CGMY Model

CGMY_simulation <- CGMY_increments(sigma=sigma_hat, nu=nu_hat, mu=0, Y=1 , T=length(x_observations), N=N)

plot 1:length(x_observations), cumsum(x_observations), ylim = c(-2,2))
lines(seq(from =1, to=length(x_observations), length.out = N), cumsum(CGMY_simulation), col = 'blue')



#------------------------------------
# comparing MSE of BM and VG
MSE <- function(bm, vg, real_data){
  mse_bm <- (bm - real_data)^2/length(real_data)
  mse_vg <- (vg - real_data)^2/length(real_data)
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
n_iter <- 100
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
p1 <- ggplot(comparison_data, aes(x=kurtosis)) + geom_histogram() + facet_grid(.~process)
p1
