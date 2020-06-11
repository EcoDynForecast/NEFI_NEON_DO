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


DO_proc <- DO_use %>% 
  filter(Datetime < as.Date('2019-06-26'))

I = DO_proc$PARMean
wind = DO_proc$windSpeedMean10m
y = DO_proc$dissolvedOxygen
y_sat = DO_proc$dissolvedOxygenSat
z = DO_proc$thermoDepth

IP = 0.001
R0 = 0.25
Ib = 4
IR = 0.0001
P = 10

mu = rep(NA,nrow(DO_proc))
mu[2] = DO_proc$dissolvedOxygen[1]

I_history <- NA

I_mat <- list()

for(t in 3:nrow(DO_proc)){
    I_history[t] = I[t-2]*(exp(-Ib*4)+exp(-Ib*3)+exp(-Ib*2))
    mu[t] = mu[t-1] + (P*(1-exp(-IP*I[t]/P)) - R0 + (I[t-1]+exp(-Ib)*I_history[t])*IR - (2.07 + 0.215*wind[t]^1.7)*(mu[t-1]-y_sat[t-1])/z[t-1])*1/24
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
  for(t in 3:n){
    Ih[t] = I[t-2]*(exp(-Ib*4)+exp(-Ib*3)+exp(-Ib*2))
    mu[t] = x[t-1] + (P*(1-exp(-IP*I[t]/P)) - R0 + (I[t-1]+exp(-Ib)*Ih[t])*IR - (2.07 + 0.215*wind[t]^1.7)*(x[t-1]-y_sat[t])/z[t])*1/24   
    x[t] ~ dnorm(mu[t],tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)         #initial state
  x[2] ~ dnorm(x[1],tau_ic)         #second step
  tau_add ~ dgamma(a_add,r_add)     #process error
  IP ~ dunif(0,1)          #IP distribution
  R0 ~ dunif(0,10)                  #R0 distribution
  IR ~ dunif(0,0.1)                 #IP distribution
  Ib ~ dunif(4,10)      #R0 distribution
  P ~ dunif(0.1,50)
}
"

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
                            variable.names = c("tau_add","IP","R0","IR","Ib","P"),
                            n.iter = 5000,
                            thin = 5)
plot(jags.out)
gelman.diag(jags.out)
gelman.plot(jags.out)

#after convergence, take more samples including x
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","IP","R0","IR","Ib","P"),
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

forecast <- function(IC,I,IP,Ib,IR,R0,P,wind,y_sat,z,n,PE){
  x <- matrix(NA,n,length(time))  ## storage
  mu = NA
  Ih = NA
  prev <- IC           ## initialize
  for(t in 1:length(time)){
    if(t == 1){
      Ii = DO_proc$PARMean[(nrow(DO_proc)-1)]
      Ii2 = DO_proc$PARMean[(nrow(DO_proc))]
    }else if (t==2){
      Ii = DO_proc$PARMean[(nrow(DO_proc))]
      Ii2 = I[t-1]
    }else{
      Ii = I[t-2]
      Ii2 = I[t]
    }
    Ih = Ii*(exp(-Ib*4)+exp(-Ib*3)+exp(-Ib*2))
    mu = prev + (P*(1-exp(-IP*I[t]/P)) - R0 + (Ii2+exp(-Ib)*Ih)*IR - (2.07 + 0.215*wind[t]^1.7)*(prev-y_sat[t])/z[t])*1/24
    x[,t] = rnorm(n,mu,PE)                        ## predict next step
    prev = x[,t]                                  ## update IC
  }
  return(x)
}

#DETERMINISTIC------------------------------

## parameters
params <- data.frame(out) %>% 
  select(-contains('x'),-contains('Ih'))
param.mean <- apply(params,2,mean)

x_df <- data.frame(out) %>% 
  select(contains('x'))

## initial conditions
IC <- as.matrix(x_df)


cast_det <- forecast(IC=mean(IC[,ncol(IC)]),
                     I=DO_cast$PARMean,
                     IP=param.mean["IP"],
                     R0=param.mean["R0"],
                     Ib=param.mean["Ib"],
                     IR=param.mean["IR"],
                     P=param.mean["P"],
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
                    Ib=param.mean["Ib"],
                    IR=param.mean["IR"],
                    P=param.mean["P"],
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

ggplot(plot_all %>% filter(as.Date(Datetime)>=as.Date('2019-06-26')),aes(Datetime,DO_obs))+
  geom_point(size=2)+
  geom_line(data=plot_all,aes(Datetime,mu))+
  geom_ribbon(aes(ymin=CI_L,ymax=CI_U,fill=ID),alpha=0.2)

#PARAMETER--------------------------------
prow = sample.int(nrow(params),1000,replace=TRUE)

cast_param <- forecast(IC=IC[prow,ncol(IC)],
                       I=DO_cast$PARMean,
                       IP=params[prow,"IP"],
                       R0=params[prow,"R0"],
                       IR=params[prow,"IR"],
                       Ib=params[prow,"Ib"],
                       P=params[prow,"P"],
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
                      Ib=params[prow,"Ib"],
                      IR=params[prow,"IR"],
                      P=params[prow,"P"],
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
  geom_line(data=plot_all %>% filter(ID=='Deterministic'|ID=='Train'),aes(Datetime,mu),size=1)+
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
