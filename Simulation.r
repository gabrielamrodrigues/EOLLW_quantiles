rm(list = ls())

##Gamlss Model
source('EOLLW_025.R') 
##Use EOLLE_050 for quantile 0.50 and EOLLW_075 for quantile 0.75
library(gamlss.cens)


# Function to generate regression model with 1 covariate
# (Use rEOLLE_050 for quantile 0.50 and rEOLLW_075 for quantile 0.75)

reg_eollw <- function(n,beta10,beta11,beta20,beta21,beta30,beta40){
  x1 <- rbinom(n,1,0.5)
  mu <- exp(beta10+beta11*x1)
  sigma <- exp(beta20+beta21*x1)
  nu <- exp(beta30)
  tau <- exp(beta40)
  y <- rEOLLW_025(n=n,mu=mu,sigma=sigma,nu=nu,tau=tau)
  return(list(y=y,x1=x1,mu=mu,sigma=sigma,nu=nu,tau=tau))
}

##Initial Values
beta10 <- 1.50;beta11 <- -1.32;beta20 <- 0.5;beta21 <- 0.2;
beta30 <-1.1; beta40 <- 1.4
initial <- c(beta10,beta11,beta20,beta21,beta30,beta40)


# Generate censored data
estimates_cens <- function(n,beta10,beta11,beta20,beta21,beta30,beta40,porc.cens,gama){
  reg <- reg_eollw(n,beta10,beta11,beta20,beta21,beta30,beta40)
  tempos <- reg$y
  x1 <- reg$x1
  mu0 <- reg$mu
  sigma0 <- reg$sigma
  nu0 <- reg$nu
  tau0 <- reg$tau
  tempo.cens <- runif(n,0,gama)
  delta <- c()
  tempo <- c()
  k <- 0
  for(j in 1:n){
    if(tempo.cens[j] <= tempos[j] && k < round(n*porc.cens/100))
    {delta[j] <- 0
    tempo[j] <- tempo.cens[j]
    k <- k+1}
    else
    {delta[j] <- 1
    tempo[j] <- tempos[j]}
  }
  dados <- data.frame(tempo,delta,x1,mu0,sigma0,tau0,nu0)
  return(dados)
}


# Likelihood function
# Use dEOLLE_050 and pEOLLE_050 for quantile 0.50 
# Use dEOLLW_075 and pEOLLE_075 for quantile 0.75

eollw_loglik <- function(par){
  n <- length(t)
  beta10 <- par[1]
  beta11 <- par[2]
  beta20 <- par[3]
  beta21 <- par[4]
  beta30 <- par[5]
  beta40 <- par[6]
  mu0 <- exp(beta10+beta11*x1)
  sigma0 <- exp(beta20+beta21*x1)
  nu0 <- exp(beta30)
  tau0 <- exp(beta40)
  f <- dEOLLW_025(t, mu=mu0, sigma=sigma0, nu=nu0,tau=tau0,log=F)
  S <- (1-pEOLLW_025(t,mu=mu0, sigma=sigma0, nu=nu0,tau=tau0))
  loglik <- sum(censura*log(f)+(1-censura)*(log(S)))
  return(loglik)
}



# Sample size n=500 and Censorship percentage: 0%

n <- 500
gama <- 10
porc.cens <- 0

# for 10% censorship use porc.cens=10
# for 10% censorship use gamma=10

# for 50% censorship use porc.cens=50
# for 50% censorship use gamma=5

npar <- 6
r <- 1000
theta <- matrix(0,r,npar)
se <- matrix(0,r,npar)
i <- 1
cens.vector <- c()
set.seed(1311)
while(i<=r){
  data1 <- estimates_cens(n,beta10,beta11,beta20,beta21,beta30,beta40,porc.cens,gama)
  x1 <- data1$x1
  t <- data1$tempo
  censura <- data1$delta
  fit1=try(optim(initial,fn=eollw_loglik,hessian=T,control = list(fnscale = -1),method = 'BFGS'))
  if((class(fit1) != "try-error")==T){
    if(fit1$counts[1]<=70){
      inf <- try(solve(fit1$hessian))   
      if(((class(inf) != "try-error")==T)){
        if((inf[1,1]<0)&(inf[2,2]<0)&(inf[3,3]<0)&(inf[4,4]<0)&(inf[5,5]<0)&(inf[6,6]<0)){
          se[i,] <- sqrt(-diag(inf))  
          theta[i,] <- fit1$par
          cens.vector[i] <- summary(as.factor(data1$delta))[1]
          i <- i+1
        }else{i <- i}
      }    
      else{i <- i}
    }
    else{i <- i}
  }
  else{i <- i}
}


theta500_cens0 <- theta
se500_cens0 <- se
cens.vector500.0 <- cens.vector


# Sample size n=300 and Censorship percentage: 0%

n <- 300
gama <- 10
porc.cens <- 0

# for 10% censorship use porc.cens=10
# for 10% censorship use gamma=10

# for 50% censorship use porc.cens=50
# for 50% censorship use gamma=5


theta300_cens0 <- theta
se300_cens0 <- se
cens.vector300.0 <- cens.vector


# Sample size n=100 and Censorship percentage: 0%

n <- 100
gama <- 10
porc.cens <- 0
# for 10% censorship use porc.cens=10
# for 10% censorship use gamma=10

# for 50% censorship use porc.cens=50
# for 50% censorship use gamma=4


theta100_cens0 <- theta
se100_cens0 <- se
cens.vector100.0 <- cens.vector


# Calculate results

# Sample size n=100 and Censorship percentage: 0%

estimates <- theta100_cens0
n <- 100

mean.estimates <- apply(estimates, 2, mean) 
sd <- apply(estimates, 2, sd)
errop <- sd/sqrt(n) 
var <- apply(estimates, 2,var) 
eqm <- c( var[1]+(mean.estimates[1]-initial[1])^2,
          var[2]+(mean.estimates[2]-initial[2])^2,
          var[3]+(mean.estimates[3]-initial[3])^2,
          var[4]+(mean.estimates[4]-initial[4])^2,
          var[5]+(mean.estimates[5]-initial[5])^2,
          var[6]+(mean.estimates[6]-initial[6])^2)

vies <- c( (mean.estimates[1]-initial[1]),
           (mean.estimates[2]-initial[2]),
           (mean.estimates[3]-initial[3]),
           (mean.estimates[4]-initial[4]),
           (mean.estimates[5]-initial[5]),
           (mean.estimates[6]-initial[6]))
results100 <- data.frame(mean.estimates, vies,eqm)


##  Sample size n=300 and Censorship percentage: 0%

estimates <- theta300_cens0
n <- 300

mean.estimates <- apply(estimates, 2, mean) 
sd <- apply(estimates, 2, sd)
errop <- sd/sqrt(n) 
var <- apply(estimates, 2,var) 
eqm <- c( var[1]+(mean.estimates[1]-initial[1])^2,
          var[2]+(mean.estimates[2]-initial[2])^2,
          var[3]+(mean.estimates[3]-initial[3])^2,
          var[4]+(mean.estimates[4]-initial[4])^2,
          var[5]+(mean.estimates[5]-initial[5])^2,
          var[6]+(mean.estimates[6]-initial[6])^2)

vies <- c( (mean.estimates[1]-initial[1]),
           (mean.estimates[2]-initial[2]),
           (mean.estimates[3]-initial[3]),
           (mean.estimates[4]-initial[4]),
           (mean.estimates[5]-initial[5]),
           (mean.estimates[6]-initial[6]))
results300 <- data.frame(mean.estimates, vies,eqm)


# Sample size n=500 and Censorship percentage: 0%

estimates <- theta500_cens0
n <- 500

mean.estimates <- apply(estimates, 2, mean) 
sd <- apply(estimates, 2, sd)
errop <- sd/sqrt(n) 
var <- apply(estimates, 2,var) 
eqm <- c( var[1]+(mean.estimates[1]-initial[1])^2,
          var[2]+(mean.estimates[2]-initial[2])^2,
          var[3]+(mean.estimates[3]-initial[3])^2,
          var[4]+(mean.estimates[4]-initial[4])^2,
          var[5]+(mean.estimates[5]-initial[5])^2,
          var[6]+(mean.estimates[6]-initial[6])^2)

vies <- c( (mean.estimates[1]-initial[1]),
           (mean.estimates[2]-initial[2]),
           (mean.estimates[3]-initial[3]),
           (mean.estimates[4]-initial[4]),
           (mean.estimates[5]-initial[5]),
           (mean.estimates[6]-initial[6]))
results500 <- data.frame(mean.estimates, vies,eqm)
