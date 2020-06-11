library(tidyverse)
library(rjags)
library(lubridate)
library(coda)
library(ecoforecastR)
library(imputeTS)

DO1 <- read.csv('CRAM_dissolved_oxygen_1hr.csv') %>% 
  mutate(site = 1)

DO1$Datetime <- ymd_hms(DO1$startDateTime)


DO2 <- read.csv('SUGG_dissolved_oxygen_1hr.csv') %>% 
  select(-siteID) %>% 
  mutate(site = 2)

DO2$Datetime <- ymd_hms(DO2$startDateTime)

ggplot(DO2,aes(Datetime,dissolvedOxygen))+
  geom_point()+
  geom_point(data=DO1,aes(Datetime,dissolvedOxygen))+
  scale_x_datetime(date_breaks = '1 month')

date_start <- as.Date('2019-9-15')
date_end <- as.Date('2019-10-15')

DO <- rbind(DO1,DO2) %>% 
  filter(Datetime >= date_start & Datetime <= date_end)


ggplot(DO,aes(Datetime,thermoDepth,color=site))+
  geom_point()

DO_tau <- DO %>% 
  mutate(tau = 1/dissolvedOxygenExpUncert^2)

DO_use <- DO_tau
for (k in 2:nrow(DO_tau)){
  if(is.na(DO_tau$tau[k])){
    DO_use$tau[k] = DO_use$tau[k-1]
  }else{
    DO_use$tau[k] = DO_tau$tau[k]
  }
}

DO_use1 <- DO_use %>% filter(site==1 & as.Date(Datetime)<as.Date('2019-10-01'))
DO_use2 <- DO_use %>% filter(site==2 & as.Date(Datetime)<as.Date('2019-10-01'))

site_DO <- matrix(data=NA,2,nrow(DO_use1))
site_DO[1,] <- DO_use1$dissolvedOxygen
site_DO[2,] <- DO_use2$dissolvedOxygen

site_tau <- matrix(data=NA,2,nrow(DO_use1))
site_tau[1,] <- DO_use1$tau
site_tau[2,] <- DO_use2$tau

site_PAR <- matrix(data=NA,2,nrow(DO_use1))
site_PAR[1,] <- DO_use1$PARMean
site_PAR[2,] <- DO_use2$PARMean

site_wind <- matrix(data=NA,2,nrow(DO_use1))
site_wind[1,] <- DO_use1$windSpeedMean10m
site_wind[2,] <- DO_use2$windSpeedMean10m

site_z <- matrix(data=NA,2,nrow(DO_use1))
site_z[1,] <- DO_use1$thermoDepth
site_z[2,] <- DO_use2$thermoDepth

site_sat <- matrix(data=NA,2,nrow(DO_use1))
site_sat[1,] <- DO_use1$dissolvedOxygenSat
site_sat[2,] <- DO_use2$dissolvedOxygenSat


RandomWalk = "
model{
  
  #### Data Model
  for(t in 1:nt){
    for(s in 1:ns){
    y[s,t] ~ dnorm(x[s,t],tau_obs[s,t])    #observation error
    }
  }
  
  #### Process Model
  for(t in 2:nt){
    for(s in 1:ns){
    x[s,t]~dnorm(x[s,t-1],tau_add[s])    #process error
    }
  }
  
  #### Priors
  for(s in 1:ns){
  x[s,1] ~ dnorm(x_ic[s],tau_ic[s])        #initial state
  tau_add[s] ~ dgamma(a_add,r_add)
  }
}
"



data <- list(y=site_DO,
             tau_obs = site_tau,
             nt=ncol(site_DO),
             ns=2,
             x_ic=c(9,5),
             tau_ic=c(1/0.1^2,1/0.1^2),
             a_add=0.1,
             r_add=0.01)

# nchain = 3
# init <- list()
# for(i in 1:nchain){
#   y.samp = sample(DO_use$dissolvedOxygen,length(DO_use$dissolvedOxygen),replace=TRUE)
#   init[[i]] <- list(tau_add=1/var(diff(log(y.samp)),na.rm=TRUE))
# }

j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         # inits = init,
                         n.chains = 3)

#assess convergence
jags.out   <- coda.samples (model = j.model,
                            variable.names = "tau_add",
                            n.iter = 5000,
                            thin = 5)
plot(jags.out)
gelman.diag(jags.out)
gelman.plot(jags.out)

#after convergence, take more samples including x
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add"),
                            n.iter = 10000)

burnin = 5000                                ## determine convergence
# jags.out <- window(jags.out,start=burnin)  ## remove burn-in

#plot confidence intervals
time.rng = c(1,length(DO_use$Datetime)) ## adjust to zoom in and out
out <- as.matrix(jags.out)
out <- out[burnin:nrow(out),]

out <- data.frame(out)
x.cols1 <- out %>% select(contains('x.1'))     
x.cols2 <- out %>% select(contains('x.2'))
ci1 <- apply(x.cols1,2,quantile,c(0.025,0.5,0.975)) 
ci2 <- apply(x.cols2,2,quantile,c(0.025,0.5,0.975)) 

ci2[which(ci2 <0)] <- 0

rwalk <- data.frame(Datetime = DO_use1$Datetime,
                    DO1 = DO_use1$dissolvedOxygen,
                    DO2 = DO_use2$dissolvedOxygen,
                    ci1_l = ci1[1,],
                    ci1_u = ci1[3,],
                    mu1 = ci1[2,],
                    ci2_l = ci2[1,],
                    ci2_u = ci2[3,],
                    mu2 = ci2[2,])

ggplot(rwalk,aes(Datetime,DO1))+
  geom_point(size=1)+
  geom_line(aes(Datetime,mu1),size=1)+
  geom_point(aes(Datetime,DO2),size=1)+
  geom_line(aes(Datetime,mu2),size=1)+
  geom_ribbon(aes(ymin=ci1_l,ymax=ci1_u),alpha=0.5)+
  geom_ribbon(aes(ymin=ci2_l,ymax=ci2_u),alpha=0.5)+
  scale_y_continuous(limits=c(0,NA))

#Add in other variables now------------------------------

#DO[t] = DO[t-1] + GPP[t-1] - R[t-1] + F[t-1]
#GPP = I*IP
#R = R0
#F = kO2(DO-DOsat)/z
#kO2 = 2.07 + 0.215*wind^1.7

do_data <- DO_test

do_pred <- rep(NA,nrow(do_data) )
do_pred[1] <- do_data$dissolvedOxygen[1]
IP = .001
R0 = 5
for(i in 2:nrow(do_data)){
  do_pred[i] <- do_pred[i-1] + (do_data$PARMean[i] * IP  - R0 - 
                                  (2.07 + 0.215 * do_data$windSpeedMean10m[i] ^ 1.7) * ((do_pred[i-1] - do_data$dissolvedOxygenSat[i-1])/do_data$thermoDepth[i-1])) * (1/(6 * 24))
}
combined <- tibble(startDateTime = do_data$startDateTime,
                   prediction = do_pred,
                   observed = do_data$dissolvedOxygen) %>% 
  pivot_longer(-startDateTime, names_to = "Value", values_to = "DO") 
ggplot(combined) +
  geom_point(aes(x = startDateTime, y = DO, color = Value))


DO_proc <- DO_use %>% 
  filter(as.Date(Datetime) < as.Date('2019-10-01'))

I = DO_proc$PARMean
wind = DO_proc$windSpeedMean10m
y = DO_proc$dissolvedOxygen
y_sat = DO_proc$dissolvedOxygenSat
z = DO_proc$thermoDepth

IP = 0.001
R0 = 10

mu = rep(NA,nrow(DO_proc))
mu[1] = y[1]

for(k in 2:nrow(DO_proc)){
  mu[k] = mu[k-1] + (I[k]*IP - R0 - (2.07 + 0.215*wind[k]^1.7)*(mu[k-1]-y_sat[k-1])/z[k-1])*1/24
}

df <- data.frame(Datetime = DO_proc$Datetime,DO = DO_proc$dissolvedOxygen,mu)

ggplot(df,aes(Datetime,mu))+
  geom_point()+
  geom_point(data=DO_use,aes(Datetime,dissolvedOxygen),color='red')

process = "
model{
  
  #### Data Model
  for(t in 1:nt){
    for(s in 1:ns){
      y[s,t] ~ dnorm(x[s,t],tau_obs[s,t])   #observation error
    }
  }

  #### Process Model
  for(t in 2:nt){
    for(s in 1:ns){
      mu[s,t] = x[s,t-1] + (I[s,t]*IP[s] - R0[s] - (2.07 + 0.215*wind[s,t]^1.7)*(x[s,t-1]-y_sat[s,t])/z[s,t])*1/24   
      x[s,t] ~ dnorm(mu[s,t],tau_add[s])
    }
  }

  #### Priors
  for(s in 1:ns){
    x[s,1] ~ dnorm(x_ic[s],tau_ic[s])         #initial state
    tau_add[s] ~ dgamma(a_add,r_add)     #process error
    IP[s] ~ dunif(0,1)          #IP distribution
    R0[s] ~ dunif(0,10)          #R0 distribution
  }
}
"

data <- list(y=site_DO, 
             y_sat = site_sat,
             I = site_PAR,
             wind = site_wind,
             z = site_z,
             tau_obs = site_tau,
             nt=ncol(site_tau),
             ns=2,
             x_ic=c(9,5),
             tau_ic=c(1/0.1^2,1/0.1^2),
             a_add=0.1,
             r_add=0.01)

# nchain = 3
# init <- list()
# for(i in 1:nchain){
#   y.samp = sample(DO_use$dissolvedOxygen,length(DO_use$dissolvedOxygen),replace=TRUE)
#   init[[i]] <- list(tau_add=1/var(diff(log(y.samp)),na.rm=TRUE))
# }

j.model   <- jags.model (file = textConnection(process),
                         data = data,
                         # inits = init,
                         n.chains = 3)

#assess convergence
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","IP","R0","alpha_site"),
                            n.iter = 5000,
                            thin = 5)
plot(jags.out)
gelman.diag(jags.out)
gelman.plot(jags.out)

#after convergence, take more samples including x
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","IP","R0"),
                            n.iter = 10000)

burnin = 5000                                ## determine convergence
# jags.out <- window(jags.out,start=burnin)  ## remove burn-in

#plot confidence intervals
time.rng = c(1,length(DO_use$Datetime)) ## adjust to zoom in and out
out <- as.matrix(jags.out)
out <- out[burnin:nrow(out),]

out <- data.frame(out)
x.cols1 <- out %>% select(contains('x.1'))     
x.cols2 <- out %>% select(contains('x.2'))
ci1 <- apply(x.cols1,2,quantile,c(0.025,0.5,0.975)) 
ci2 <- apply(x.cols2,2,quantile,c(0.025,0.5,0.975)) 

plot_model <- data.frame(Datetime = DO_use1$Datetime,
                    DO_obs1 = DO_use1$dissolvedOxygen,
                    DO_obs2 = DO_use2$dissolvedOxygen,
                    ci1_l = ci1[1,],
                    ci1_u = ci1[3,],
                    mu1 = ci1[2,],
                    ci2_l = ci2[1,],
                    ci2_u = ci2[3,],
                    mu2 = ci2[2,],
                    ID='Train')

ggplot(plot_model,aes(Datetime,DO_obs1))+
  geom_point(size=1)+
  geom_line(aes(Datetime,mu1),size=1)+
  geom_point(aes(Datetime,DO_obs2),size=1)+
  geom_line(aes(Datetime,mu2),size=1)+
  geom_ribbon(aes(ymin=ci1_l,ymax=ci1_u),alpha=0.5)+
  geom_ribbon(aes(ymin=ci2_l,ymax=ci2_u),alpha=0.5)+
  scale_y_continuous(limits=c(0,NA))


#FORECAST-------------------------------


### settings
n = 1000         ## set number of Monte Carlo draws
DO_cast1 <- DO_use %>% filter(site==1 & as.Date(Datetime)>=as.Date('2019-10-01'))
DO_cast2 <- DO_use %>% filter(site==2 & as.Date(Datetime)>=as.Date('2019-10-01'))


cast_DO <- matrix(data=NA,2,nrow(DO_cast1))
cast_DO[1,] <- DO_cast1$dissolvedOxygen
cast_DO[2,] <- DO_cast2$dissolvedOxygen

cast_tau <- matrix(data=NA,2,nrow(DO_cast1))
cast_tau[1,] <- DO_cast1$tau
cast_tau[2,] <- DO_cast2$tau

cast_PAR <- matrix(data=NA,2,nrow(DO_cast1))
cast_PAR[1,] <- DO_cast1$PARMean
cast_PAR[2,] <- DO_cast2$PARMean
cast_PAR[which(cast_PAR<0)] <- 0

cast_wind <- matrix(data=NA,2,nrow(DO_cast1))
cast_wind[1,] <- DO_cast1$windSpeedMean10m
cast_wind[2,] <- DO_cast2$windSpeedMean10m

cast_z <- matrix(data=NA,2,nrow(DO_cast1))
cast_z[1,] <- DO_cast1$thermoDepth
cast_z[2,] <- DO_cast2$thermoDepth

cast_sat <- matrix(data=NA,2,nrow(DO_cast1))
cast_sat[1,] <- DO_cast1$dissolvedOxygenSat
cast_sat[2,] <- DO_cast2$dissolvedOxygenSat

##` @param IC    Initial Conditions
##` @param r     Intrinsic growth rate
##` @param Kg    Across-site ('global') mean carrying capacity
##` @param alpha Site random effect
##` @param beta  Slope of precipitation effect on K
##` @param ppt   Precipitation forecast
##` @param Q     Process error (default = 0 for deterministic runs)
##` @param n     Size of Monte Carlo ensemble
forecast1 <- function(IC,I,IP,R0,wind,y_sat,z,n,PE){
  x <- list()
  x[[1]] <- matrix(NA,n,nrow(DO_cast1))  ## storage
  x[[2]] <- matrix(NA,n,nrow(DO_cast1))
  mu = matrix(NA,n,1)
  prev <- IC           ## initialize
  for(t in 1:nrow(DO_cast1)){
    for(s in 1:2){
      mu = prev[s] + (I[s,t]*IP[s] - R0[s] - (2.07 + 0.215*wind[s,t]^1.7)*(prev[s]-y_sat[s,t])/z[s,t])*1/24
      x[[s]][,t] = rnorm(n,mu,PE[s])                        ## predict next step                                 ## update IC
    }
    prev = c(x[[1]][,t],x[[2]][,t])
  }
  return(x)
}

#DETERMINISTIC------------------------------

## parameters
params <- data.frame(out) %>% 
  select(-contains('x'))

IP_param <- params %>% 
  select(1:2)
IP_mean <- colMeans(IP_param)

R0_param <- params %>% 
  select(3:4)
R0_mean <- colMeans(R0_param)

tau_param <- params %>% 
  select(5:6)
tau_mean <- colMeans(tau_param)

x_df1 <- data.frame(out) %>% 
  select(contains('x.1'))
x_df2 <- data.frame(out) %>% 
  select(contains('x.2'))
## initial conditions

IC1 <- as.matrix(x_df1)
IC2 <- as.matrix(x_df2)
IC <- list(as.matrix(x_df1),as.matrix(x_df2))

cast_det <- forecast1(IC=c(mean(IC[[1]][,ncol(x_df1)]),mean(IC[[2]][,ncol(x_df2)])),
                   I=cast_PAR,
                   IP=IP_mean,
                   R0=R0_mean,
                   PE=c(0,0),
                   wind = cast_wind,
                   y_sat = cast_sat,
                   z = cast_z,
                   n=1)


# cast_ci <- apply(cast,2,quantile,c(0.025,0.5,0.975)) 

plot_det <- data.frame(Datetime = DO_cast1$Datetime,
                       DO_obs1 = DO_cast1$dissolvedOxygen,
                       DO_obs2 = DO_cast2$dissolvedOxygen,
                       ci1_l = rep(NA,nrow(DO_cast1)),
                       ci1_u = rep(NA,nrow(DO_cast1)),
                       mu1 = colMeans(cast_det[[1]]),
                       ci2_l = rep(NA,nrow(DO_cast1)),
                       ci2_u = rep(NA,nrow(DO_cast1)),
                       mu2 = colMeans(cast_det[[2]]),
                       ID='Deterministic')

ggplot(plot_det,aes(Datetime,mu1))+
  geom_line()+
  geom_line(aes(Datetime,mu2))

plot_all <- rbind(plot_model,plot_det) %>% 
  mutate(group = case_when(as.Date(Datetime)<as.Date('2019-06-26')~'Train',
                           as.Date(Datetime)>=as.Date('2019-06-26')~'Cast'))

ggplot(plot_all,aes(Datetime,DO_obs1))+
  geom_point(size=1)+
  geom_point(aes(Datetime,DO_obs2),size=1)+
  geom_line(aes(Datetime,mu1))+
  geom_line(aes(Datetime,mu2))+
  geom_ribbon(aes(ymin=ci1_l,ymax=ci1_u,fill=group),alpha=0.2)+
  geom_ribbon(aes(ymin=ci2_l,ymax=ci2_u,fill=group),alpha=0.2)


#INITIAL CONDITION--------------------------------
forecast2 <- function(IC,I,IP,R0,wind,y_sat,z,n,PE){
  x <- list()
  x[[1]] <- matrix(NA,n,nrow(DO_cast1))  ## storage
  x[[2]] <- matrix(NA,n,nrow(DO_cast1))
  mu = matrix(NA,n,1)
  prev <- IC           ## initialize
  for(t in 1:nrow(DO_cast1)){
    for(s in 1:2){
      mu = prev[s,] + (I[s,t]*IP[s] - R0[s] - (2.07 + 0.215*wind[s,t]^1.7)*(prev[s,]-y_sat[s,t])/z[s,t])*1/24
      x[[s]][,t] = rnorm(n,mu,PE[s])                        ## predict next step                                 ## update IC
    }
    prev = data.frame(x[[1]][,t],x[[2]][,t]) %>% t()
  }
  return(x)
}


prow = sample.int(nrow(params),1000,replace=TRUE)


cast_IC <- forecast2(IC=data.frame(IC[[1]][prow,ncol(x_df1)],IC[[2]][prow,ncol(x_df1)]) %>% t(),
                    I=cast_PAR,
                    IP=IP_mean,
                    R0=R0_mean,
                    PE=c(0,0),
                    wind = cast_wind,
                    y_sat = cast_sat,
                    z = cast_z,
                    n=1000)

IC_ci1 <- apply(cast_IC[[1]],2,quantile,c(0.025,0.5,0.975)) 
IC_ci2 <- apply(cast_IC[[2]],2,quantile,c(0.025,0.5,0.975)) 

plot_IC <- data.frame(Datetime = DO_cast1$Datetime,
                       DO_obs1 = DO_cast1$dissolvedOxygen,
                       DO_obs2 = DO_cast2$dissolvedOxygen,
                       ci1_l = IC_ci1[1,],
                       ci1_u = IC_ci1[3,],
                       mu1 = IC_ci1[2,],
                       ci2_l = IC_ci2[1,],
                       ci2_u = IC_ci2[3,],
                       mu2 = IC_ci2[2,],
                       ID='Initial Condition')

ggplot(plot_IC,aes(Datetime,mu1))+
  geom_line()+
  geom_line(aes(Datetime,mu2))

plot_all <- rbind(plot_model,plot_det,plot_IC) %>% 
  mutate(group = case_when(as.Date(Datetime)<as.Date('2019-06-26')~'Train',
                           as.Date(Datetime)>=as.Date('2019-06-26')~'Cast'))

ggplot(plot_all,aes(Datetime,DO_obs1))+
  geom_point(size=1)+
  geom_point(aes(Datetime,DO_obs2),size=1)+
  geom_line(aes(Datetime,mu1))+
  geom_line(aes(Datetime,mu2))+
  geom_ribbon(aes(ymin=ci1_l,ymax=ci1_u,fill=ID),alpha=0.2)+
  geom_ribbon(aes(ymin=ci2_l,ymax=ci2_u,fill=ID),alpha=0.2)

#PARAMETER--------------------------------
forecast3 <- function(IC,I,IP,R0,wind,y_sat,z,n,PE){
  x <- list()
  x[[1]] <- matrix(NA,n,nrow(DO_cast1))  ## storage
  x[[2]] <- matrix(NA,n,nrow(DO_cast1))
  mu = matrix(NA,n,1)
  prev <- IC           ## initialize
  for(t in 1:nrow(DO_cast1)){
    for(s in 1:2){
      mu = prev[s,] + (I[s,t]*IP[s,] - R0[s,] - (2.07 + 0.215*wind[s,t]^1.7)*(prev[s,]-y_sat[s,t])/z[s,t])*1/24
      x[[s]][,t] = rnorm(n,mu,PE[s])                        ## predict next step                                 ## update IC
    }
    prev = data.frame(x[[1]][,t],x[[2]][,t]) %>% t()
  }
  return(x)
}

prow = sample.int(nrow(params),1000,replace=TRUE)

cast_param <- forecast3(IC=data.frame(IC[[1]][prow,ncol(x_df1)],IC[[2]][prow,ncol(x_df1)]) %>% t(),
                    I=cast_PAR,
                    IP=IP_param[prow,] %>% t(),
                    R0=R0_param[prow,] %>% t(),
                    PE=c(0,0),
                    wind = cast_wind,
                    y_sat = cast_sat,
                    z = cast_z,
                    n=1000)

param_ci1 <- apply(cast_param[[1]],2,quantile,c(0.025,0.5,0.975)) 
param_ci2 <- apply(cast_param[[2]],2,quantile,c(0.025,0.5,0.975)) 

plot_param <- data.frame(Datetime = DO_cast1$Datetime,
                      DO_obs1 = DO_cast1$dissolvedOxygen,
                      DO_obs2 = DO_cast2$dissolvedOxygen,
                      ci1_l = param_ci1[1,],
                      ci1_u = param_ci1[3,],
                      mu1 = param_ci1[2,],
                      ci2_l = param_ci2[1,],
                      ci2_u = param_ci2[3,],
                      mu2 = param_ci2[2,],
                      ID='Parameter')

ggplot(plot_param,aes(Datetime,mu1))+
  geom_line()+
  geom_line(aes(Datetime,mu2))+
  geom_ribbon(aes(ymin=ci2_l,ymax=ci2_u,fill=ID),alpha=0.2)

plot_all <- rbind(plot_model,plot_det,plot_IC,plot_param) %>% 
  mutate(group = case_when(as.Date(Datetime)<as.Date('2019-06-26')~'Train',
                           as.Date(Datetime)>=as.Date('2019-06-26')~'Cast'))

ggplot(plot_all,aes(Datetime,DO_obs1))+
  geom_point(size=1)+
  geom_point(aes(Datetime,DO_obs2),size=1)+
  geom_line(aes(Datetime,mu1))+
  geom_line(aes(Datetime,mu2))+
  geom_ribbon(aes(ymin=ci1_l,ymax=ci1_u,fill=ID),alpha=0.2)+
  geom_ribbon(aes(ymin=ci2_l,ymax=ci2_u,fill=ID),alpha=0.2)
#PROCESS--------------------------------
prow = sample.int(nrow(params),1000,replace=TRUE)

cast_proc <- forecast3(IC=data.frame(IC[[1]][prow,ncol(x_df1)],IC[[2]][prow,ncol(x_df1)]) %>% t(),
                       I=cast_PAR,
                       IP=IP_param[prow,] %>% t(),
                       R0=R0_param[prow,] %>% t(),
                       PE=1/sqrt(tau_param[prow,] %>% t()),
                       wind = cast_wind,
                       y_sat = cast_sat,
                       z = cast_z,
                       n=1000)

proc_ci1 <- apply(cast_proc[[1]],2,quantile,c(0.025,0.5,0.975)) 
proc_ci2 <- apply(cast_proc[[2]],2,quantile,c(0.025,0.5,0.975)) 

plot_proc <- data.frame(Datetime = DO_cast1$Datetime,
                         DO_obs1 = DO_cast1$dissolvedOxygen,
                         DO_obs2 = DO_cast2$dissolvedOxygen,
                         ci1_l = proc_ci1[1,],
                         ci1_u = proc_ci1[3,],
                         mu1 = proc_ci1[2,],
                         ci2_l = proc_ci2[1,],
                         ci2_u = proc_ci2[3,],
                         mu2 = proc_ci2[2,],
                         ID='Process')

ggplot(plot_proc,aes(Datetime,mu1))+
  geom_line()+
  geom_line(aes(Datetime,mu2))+
  geom_ribbon(aes(ymin=ci2_l,ymax=ci2_u,fill=ID),alpha=0.2)

plot_all <- rbind(plot_model,plot_det,plot_IC,plot_param,plot_proc) %>% 
  mutate(group = case_when(as.Date(Datetime)<as.Date('2019-06-26')~'Train',
                           as.Date(Datetime)>=as.Date('2019-06-26')~'Cast'))

ggplot(plot_all,aes(Datetime,DO_obs1))+
  geom_point(size=1)+
  geom_point(aes(Datetime,DO_obs2),size=1)+
  geom_line(aes(Datetime,mu1))+
  geom_line(aes(Datetime,mu2))+
  geom_ribbon(aes(ymin=ci1_l,ymax=ci1_u,fill=ID),alpha=0.2)+
  geom_ribbon(aes(ymin=ci2_l,ymax=ci2_u,fill=ID),alpha=0.2)
plot_all$ID <- factor(plot_all$ID,
                         levels=c('Deterministic',
                                  'Process',
                                  'Parameter',
                                  'Initial Condition',
                                  'Train'))

pal <- c("Deterministic" = "NA",
         "Process" = "#FF9999",
         "Parameter" = "#9999FF",
         "Initial Condition" = "#99FF99",
         "Train" = "888888")

ggplot(plot_all,aes(Datetime,DO_obs1))+
  geom_ribbon(aes(ymin=ci1_l,ymax=ci1_u,fill=ID),alpha=0.5)+
  geom_point(size=0.5)+
  geom_line(aes(Datetime,mu1),size=1)+
  geom_ribbon(aes(ymin=ci2_l,ymax=ci2_u,fill=ID),alpha=0.5)+
  geom_point(aes(Datetime,DO_obs2),size=0.5)+
  geom_line(aes(Datetime,mu2),size=1)+
  scale_fill_manual(name='Error',
                    breaks=c('Process','Parameter','Initial Condition'),
                    values = pal)+
  theme_classic()+
  theme(text=element_text(size=12),
        legend.position = c(.1,.2))+
  labs(y='DO (mg/L)',
       x='Day')+
  scale_x_datetime(breaks = '1 day',
                   labels = date_format("%d"))
 
