library(tidyverse)
library(rjags)
library(lubridate)
library(coda)
library(ecoforecastR)
library(imputeTS)

DO <- read.csv('CRAM_dissolved_oxygen_1hr_subset.csv')

DO$Datetime <- ymd_hms(DO$startDateTime)

ggplot(DO,aes(Datetime,dissolvedOxygen))+
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


#Confidence interval function
ciEnvelope <- function (x, ylo, yhi, ...) 
{
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi), ylo[1])), 
          border = NA, ...)
}

RandomWalk = "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs[t])    #observation error
  }
  
  #### Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)    #process error
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)        #initial state
  tau_add ~ dgamma(a_add,r_add)
}
"
hist(DO$dissolvedOxygen)
hist(log(DO$dissolvedOxygen)) 

data <- list(y=DO_use$dissolvedOxygen,
             tau_obs = DO_use$tau,
             n=nrow(DO_use),
             x_ic=5,
             tau_ic=1/0.1^2,
             a_add=1,
             r_add=1)

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(DO_use$dissolvedOxygen,length(DO_use$dissolvedOxygen),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp)),na.rm=TRUE))
}

j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         inits = init,
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

burnin = 1000                                ## determine convergence
# jags.out <- window(jags.out,start=burnin)  ## remove burn-in

#plot confidence intervals
time.rng = c(1,length(DO_use$Datetime)) ## adjust to zoom in and out
out <- as.matrix(jags.out)
out <- out[burnin:nrow(out),]
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) 

plot(DO_use$Datetime,ci[2,],type='n',ylim=range(DO_use$dissolvedOxygen,na.rm=TRUE),ylab="DO",xlim=DO_use$Datetime[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ciEnvelope(DO_use$Datetime,ci[1,],ci[3,],col='#9999CC')
points(DO_use$Datetime,DO_use$dissolvedOxygen,pch="+",cex=0.5)


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
  filter(Datetime < as.Date('2019-06-26'))

I = DO_proc$PARMean
wind = DO_proc$windSpeedMean10m
y = DO_proc$dissolvedOxygen
y_sat = DO_proc$dissolvedOxygenSat
z = DO_proc$thermoDepth

IP = 0.001
R0 = 0.25

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
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs[t])   #observation error
  }

  #### Process Model
  for(t in 2:n){
    mu[t] = x[t-1] + (I[t]*IP - R0 - (2.07 + 0.215*wind[t]^1.7)*(x[t-1]-y_sat[t])/z[t])*1/24   
    x[t] ~ dnorm(mu[t],tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)         #initial state
  tau_add ~ dgamma(a_add,r_add)     #process error
  IP ~ dunif(0,1)          #IP distribution
  R0 ~ dunif(0,10)          #R0 distribution
}
"
# IP ~ dunif(0.0,1.0)
# Ro ~ dunif(0.0,50.0)
# 
# # mu_IP = 0.01,
# # tau_IP = 0.167,
# # mu_R0 = 0.25,
# # tau_R0 = 1.67
# 
# 
# tau_add <- 1 / pow(sd_add,2)
# sd_add ~ dunif(lower_add,upper_add)


data <- list(y=DO_proc$dissolvedOxygen, 
             y_sat = DO_proc$dissolvedOxygenSat,
             I = DO_proc$PARMean,
             wind = DO_proc$windSpeedMean10m,
             z = DO_proc$thermoDepth,
             tau_obs = DO_proc$tau,
             n=nrow(DO_proc),
             x_ic=9,
             tau_ic=1/0.1^2,
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
                            variable.names = c("tau_add","IP","R0"),
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
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) 

plot(DO_proc$Datetime,ci[2,],type='n',ylim=range(DO_proc$dissolvedOxygen,na.rm=TRUE),ylab="DO",xlim=DO_use$Datetime[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ciEnvelope(DO_proc$Datetime,ci[1,],ci[3,],col='#9999CC')
points(DO_proc$Datetime,DO_proc$dissolvedOxygen,pch="+",cex=0.5)

plot_model <- data.frame(CI_L = ci[1,],
                         CI_U = ci[3,],
                         mu = ci[2,],
                         Datetime = DO_proc$Datetime,
                         DO_obs = DO_proc$dissolvedOxygen,
                         ID = 'Train')

ggplot(plot_model,aes(Datetime,DO_obs))+
  geom_point(size=2)+
  geom_line(aes(Datetime,mu))+
  geom_ribbon(aes(ymin=CI_L,ymax=CI_U),alpha=0.2,fill='blue')


#FORECAST-------------------------------


### settings
n = 1000         ## set number of Monte Carlo draws
DO_cast <- DO_use %>% 
  filter(as.Date(Datetime) >= as.Date('2019-06-26'))
time = as.integer(as.factor(as.character(DO_cast$Datetime)))

##` @param IC    Initial Conditions
##` @param r     Intrinsic growth rate
##` @param Kg    Across-site ('global') mean carrying capacity
##` @param alpha Site random effect
##` @param beta  Slope of precipitation effect on K
##` @param ppt   Precipitation forecast
##` @param Q     Process error (default = 0 for deterministic runs)
##` @param n     Size of Monte Carlo ensemble
forecast <- function(IC,I,IP,R0,wind,y_sat,z,n,PE){
  x <- matrix(NA,n,length(time))  ## storage
  mu = NA
  prev <- IC           ## initialize
  for(k in 1:length(time)){
    mu = prev + (I[k]*IP - R0 - (2.07 + 0.215*wind[k]^1.7)*(prev-y_sat[k])/z[k])*1/24
    x[,k] = rnorm(n,mu,PE)                        ## predict next step
    prev = x[,k]                                  ## update IC
  }
  return(x)
}

#DETERMINISTIC------------------------------

## parameters
params <- data.frame(out) %>% 
  select(-contains('x'))
param.mean <- apply(params,2,mean)

x_df <- data.frame(out) %>% 
  select(contains('x'))
## initial conditions
IC <- as.matrix(x_df)

cast_det <- forecast(IC=mean(IC[,ncol(IC)]),
                   I=DO_cast$PARMean,
                   IP=param.mean["IP"],
                   R0=param.mean["R0"],
                   PE=0,
                   wind = DO_cast$windSpeedMean10m,
                   y_sat = DO_cast$dissolvedOxygenSat,
                   z = DO_cast$thermoDepth,
                   n=1)

# cast_ci <- apply(cast,2,quantile,c(0.025,0.5,0.975)) 

plot_det <- data.frame(CI_L = NA,
                       CI_U = NA,
                       mu = colMeans(cast_det),
                       Datetime = DO_cast$Datetime,
                       DO_obs = DO_cast$dissolvedOxygen,
                       ID = 'Deterministic')

ggplot(plot_det,aes(Datetime,mu))+
  geom_line()

plot_all <- rbind(plot_model,plot_det) %>% 
  mutate(group = case_when(as.Date(Datetime)<as.Date('2019-06-26')~'Train',
                           as.Date(Datetime)>=as.Date('2019-06-26')~'Cast'))

ggplot(plot_all,aes(Datetime,DO_obs))+
  geom_point(size=2)+
  geom_line(aes(Datetime,mu))+
  geom_ribbon(aes(ymin=CI_L,ymax=CI_U,fill=group),alpha=0.2)


#INITIAL CONDITION--------------------------------
prow = sample.int(nrow(params),1000,replace=TRUE)


cast_IC <- forecast(IC=IC[prow,ncol(IC)],
                     I=DO_cast$PARMean,
                     IP=param.mean["IP"],
                     R0=param.mean["R0"],
                     PE=0,
                     wind = DO_cast$windSpeedMean10m,
                     y_sat = DO_cast$dissolvedOxygenSat,
                     z = DO_cast$thermoDepth,
                     n=1000)

IC_ci <- apply(cast_IC,2,quantile,c(0.025,0.5,0.975)) 

plot_IC <- data.frame(CI_L = IC_ci[1,],
                       CI_U = IC_ci[3,],
                       mu = IC_ci[2,],
                       Datetime = DO_cast$Datetime,
                       DO_obs = DO_cast$dissolvedOxygen,
                       ID = 'Initial Condition')

ggplot(plot_IC,aes(Datetime,mu))+
  geom_line()

plot_all <- rbind(plot_model,plot_det,plot_IC) %>% 
  mutate(group = case_when(as.Date(Datetime)<as.Date('2019-06-26')~'Train',
                           as.Date(Datetime)>=as.Date('2019-06-26')~'Cast'))

ggplot(plot_all,aes(Datetime,DO_obs))+
  geom_point(size=2)+
  geom_line(aes(Datetime,mu))+
  geom_ribbon(aes(ymin=CI_L,ymax=CI_U,fill=ID),alpha=0.2)

#PARAMETER--------------------------------
prow = sample.int(nrow(params),1000,replace=TRUE)

forecastN(IC=IC[prow,"N[6,30]"],  ## sample IC
          r=params[prow,"r_global"],  ## sample parameters
          Kg=params[prow,"K_global"],
          alpha=params[prow,"alpha_site[6]"],
          beta=params[prow,"beta"],
          ppt=ppt.mean,
          Q=0,
          n=Nmc)

cast_param <- forecast(IC=IC[prow,ncol(IC)],
                    I=DO_cast$PARMean,
                    IP=params[prow,"IP"],
                    R0=params[prow,"R0"],
                    PE=0,
                    wind = DO_cast$windSpeedMean10m,
                    y_sat = DO_cast$dissolvedOxygenSat,
                    z = DO_cast$thermoDepth,
                    n=1000)

param_ci <- apply(cast_param,2,quantile,c(0.025,0.5,0.975)) 

plot_param <- data.frame(CI_L = param_ci[1,],
                      CI_U = param_ci[3,],
                      mu = param_ci[2,],
                      Datetime = DO_cast$Datetime,
                      DO_obs = DO_cast$dissolvedOxygen,
                      ID = 'Parameter')


plot_all <- rbind(plot_model,plot_det,plot_IC,plot_param) %>% 
  mutate(group = case_when(as.Date(Datetime)<as.Date('2019-06-26')~'Train',
                           as.Date(Datetime)>=as.Date('2019-06-26')~'Cast'))

ggplot(plot_all,aes(Datetime,DO_obs))+
  geom_point(size=2)+
  geom_line(aes(Datetime,mu))+
  geom_ribbon(aes(ymin=CI_L,ymax=CI_U,fill=ID),alpha=0.2)

#PROCESS--------------------------------
prow = sample.int(nrow(params),1000,replace=TRUE)

cast_proc <- forecast(IC=IC[prow,ncol(IC)],
                       I=DO_cast$PARMean,
                       IP=params[prow,"IP"],
                       R0=params[prow,"R0"],
                       PE=1/sqrt(params[prow,"tau_add"]),
                       wind = DO_cast$windSpeedMean10m,
                       y_sat = DO_cast$dissolvedOxygenSat,
                       z = DO_cast$thermoDepth,
                       n=1000)

proc_ci <- apply(cast_proc,2,quantile,c(0.025,0.5,0.975)) 

plot_proc <- data.frame(CI_L = proc_ci[1,],
                         CI_U = proc_ci[3,],
                         mu = proc_ci[2,],
                         Datetime = DO_cast$Datetime,
                         DO_obs = DO_cast$dissolvedOxygen,
                         ID = 'Process')


plot_all <- rbind(plot_model,plot_det,plot_IC,plot_param,plot_proc) %>% 
  mutate(group = case_when(as.Date(Datetime)<as.Date('2019-06-26')~'Train',
                           as.Date(Datetime)>=as.Date('2019-06-26')~'Cast'))
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

ggplot(plot_all,aes(Datetime,DO_obs))+
  geom_ribbon(aes(ymin=CI_L,ymax=CI_U,fill=ID),alpha=0.5)+
  geom_point(size=1)+
  geom_line(aes(Datetime,mu),size=1)+
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
  
