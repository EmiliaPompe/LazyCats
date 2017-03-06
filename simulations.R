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
# generating from the 
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




# curtosis, skewness


BM_inc <- function(sigma, T, N, theta =0){
  drift <- theta*seq(from = 0, to=T, length.out = N)
  B <- sigma*rnorm(N,0,1)*sqrt((T/N)) + drift
}
