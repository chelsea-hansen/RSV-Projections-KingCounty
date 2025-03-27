rm(list=ls())
library(deSolve)
library(tidyverse)
library(zoo)
library(imputeTS)
library(lhs)
library(tidycensus)
library(cdcfluview)
library(cowplot)
library(pracma)
library(lhs)
library(bbmle)
library(cNORM)

source("R/MSIRS_immunization_dynamics.R")

data = readRDS("DATA/fixed_parameters.rds")
parmset = data[[1]]
yinit=data[[2]]
yinit.vector=data[[3]]


timeseries = readRDS("DATA/time_series_public.rds") #data has been adjusted to account for testing changes during the pandemic 

dates = timeseries %>% filter(date<'2023-10-01') %>% select(date)

time_series_full = timeseries %>% filter(date<'2023-10-01') %>% 
  mutate(hosp = round(hosp_rate/100000*2267000))
time_series_full = c(time_series_full$hosp)
plot(time_series_full)

time_series_pre = timeseries %>% filter(date<'2020-04-01') %>% 
  mutate(hosp = round(hosp_rate/100000*2267000))

time_series_pre = c(time_series_pre$hosp)
plot(time_series_pre)

time_series_post = timeseries %>% filter(date>='2020-04-01',date<'2023-10-01') %>% 
  mutate(hosp = round(hosp_rate/100000*2267000))

time_series_post = c(time_series_post$hosp)
plot(time_series_post)



age_dist = readRDS("DATA/age_distribution_public.rds")
age_dist = c(age_dist$proportion)
age_dist

#no immunizations for calibration 
parmset$monoclonal_birth =  rep(0, nrow(timeseries) + 104)
parmset$monoclonal_catchup_01 =  rep(0, nrow(timeseries) + 104)
parmset$monoclonal_catchup_23 =  rep(0, nrow(timeseries) + 104)
parmset$monoclonal_catchup_45 =  rep(0, nrow(timeseries) + 104)
parmset$monoclonal_catchup_67 =  rep(0, nrow(timeseries) + 104)
parmset$maternal_vax <- rep(0, nrow(timeseries) + 104)
parmset$senior_vax_75 <- rep(0, nrow(timeseries) + 104)
parmset$senior_vax_60_74 <- rep(0, nrow(timeseries) + 104)
parmset$introductions <- rep(parmset$seed, nrow(timeseries) + 104)
parmset$npi <- rep(1, nrow(timeseries) + 104)


# Fit to pre-pandemic years with 100 starting values ----------------------------------------------
start_lower = c(1.94,-4.6,0.92,-4.6,-9.2,-9.6,-9.2,-9.2,-6.9)
start_upper = c(2.99,-0.1,1.95,-2.1,-4.6,-4.9,-4.6,-4.6,4.6)

conf_int = cbind(start_lower, start_upper)
set.seed(123)

h=100
lhs<-maximinLHS(h,9)

grid_parms <- cbind(
  beta = lhs[,1]*(conf_int[1,2]-conf_int[1,1])+conf_int[1,1],
  b1 = lhs[,2]*(conf_int[2,2]-conf_int[2,1])+conf_int[2,1],
  phi = lhs[,3]*(conf_int[3,2]-conf_int[3,1])+conf_int[3,1],
  RI = lhs[,4]*(conf_int[4,2]-conf_int[4,1])+conf_int[4,1],
  RC = lhs[,5]*(conf_int[5,2]-conf_int[5,1])+conf_int[5,1],
  RA = lhs[,6]*(conf_int[6,2]-conf_int[6,1])+conf_int[6,1],
  RS60 = lhs[,7]*(conf_int[7,2]-conf_int[7,1])+conf_int[7,1],
  RS75 = lhs[,8]*(conf_int[8,2]-conf_int[8,1])+conf_int[8,1],
  theta = lhs[,9]*(conf_int[9,2]-conf_int[9,1])+conf_int[9,1])


fit_times <- seq(1, length(time_series_pre) + 104, by = 1)
parameter_sets=data.frame()
fits = list()

for(g in 1:nrow(grid_parms)){
  fitmodel <-  function(parameters,dat) {
    protrans <- parameters[1] # parameter for baseline transmission rate
    baseline.txn.rate = exp(protrans)
    amp <- parameters[2] # parameter for seasonal amplitude
    b1 <-  exp(amp)
    trans <- parameters[3] # parameter for seasonal phase
    phi <-  (2*pi*(exp(trans))) / (1+exp(trans)) # transform to between 0 and 2pi
    #Age-specific reporting fractions 
    report_infants <- 1 / (1 + exp(-parameters[4])) 
    report_children <- 1 / (1 + exp(-parameters[5]))
    report_adults <- 1 / (1 + exp(-parameters[6]))
    report_seniors60 <- 1 / (1 + exp(-parameters[7]))
    report_seniors75 <- 1 / (1 + exp(-parameters[8]))
    theta <- 1 / (1 + exp(-parameters[9]))
    
    results <- ode(y=yinit.vector, method = "ode45", times=fit_times,
                   func=MSIRS_immunization_dynamics,
                   parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi))
    
    
    t0 <- length(time_series_pre)
    al <- nrow(yinit)
    output <- tail(results,t0)
    St <- output[,-1]
    I1 <- St[,grep('I1', colnames(St))]
    I2 <- St[,grep('I2', colnames(St))]
    I3 <- St[,grep('I3', colnames(St))]
    I4 <- St[,grep('I4', colnames(St))]
    R1 <- St[,grep('R1', colnames(St))]
    R2 <- St[,grep('R2', colnames(St))]
    R3 <- St[,grep('R3', colnames(St))]
    R4 <- St[,grep('R4', colnames(St))]
    S1 <- St[,grep('S1', colnames(St))]
    S2 <- St[,grep('S2', colnames(St))]
    S3 <- St[,grep('S3', colnames(St))]
    S0 <- St[,grep('S0', colnames(St))]
    M0 <- St[,grep('M0', colnames(St))]
    Si<- St[,grep('Si', colnames(St))]
    Mn1<- St[,grep('Mn1', colnames(St))]
    Mn2<- St[,grep('Mn2', colnames(St))]
    Mv1<- St[,grep('Mv1', colnames(St))]
    Mv2<- St[,grep('Mv2', colnames(St))]
    N1<- St[,grep('N1', colnames(St))]
    N2<- St[,grep('N2', colnames(St))]
    Vs1<- St[,grep('Vs1', colnames(St))]
    Vs2<- St[,grep('Vs2', colnames(St))]
    
    contact2 = parmset$npi[105:length( parmset$npi)]
    intro2 = parmset$introductions[105:length( parmset$introductions)]  
    
    lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
    for (t in 1:t0){
      beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2*contact2[t]
      lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}
    
    
    hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
    hosp2 <- hosp1 * 0.4
    hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors60, report_seniors75)
    
    H1 <- matrix(0, nrow = t0, ncol = al)
    for (i in 1:al) {
      H1[, i] <-
        hosp1[i] * parmset$RRHm * parmset$sigmaM * M0[, i] * lambda1[, i] +
        hosp1[i] * parmset$RRHm * parmset$RRHn1 * parmset$sigmaM * Mn1[, i] * lambda1[, i] +
        hosp1[i] * parmset$RRHm * parmset$RRHn2 *  Mn2[, i] * lambda1[, i] +
        hosp1[i] * parmset$RRHm * parmset$RRHv1 * parmset$sigmaM * Mv1[, i] * lambda1[, i] +
        hosp1[i] * parmset$RRHm * parmset$RRHv2 * Mv1[, i] * lambda1[, i] +
        hosp1[i] * parmset$RRHn1 * N1[, i] * lambda1[, i] +
        hosp1[i] * parmset$RRHn2 * N2[, i] * lambda1[, i] +
        hosp1[i] * S0[, i] * lambda1[, i] +
        hosp1[i] * Si[, i] * lambda1[, i] +
        hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
        hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
        hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
        hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs1[, i] * lambda1[, i] +
        hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs2[, i] * lambda1[, i]
    }
    
    H <- rowSums(H1)
    
    H2 <- cbind(H1[,1],H1[,2],H1[,3],H1[,4],H1[,5],H1[,6], rowSums(H1[, 7:8]), rowSums(H1[, 9:12]), H1[, 13],H1[,14])
    age_dist2 <- colSums(H2)
    
    LLall <- sum(dnbinom(x = dat, mu =H, size=1/theta, log = TRUE))
    LLmulti <- dmultinom(x = age_dist2, prob = age_dist, log = TRUE)
    
    LL <- LLall + LLmulti
    
    
    return(LL)
  }
  
  fitLL <- optim(par = grid_parms[g,],
                 fn = fitmodel, 
                 dat = time_series_pre,  
                 control = list(fnscale=-1, maxit=5000),
                 hessian = TRUE)
  fits[[g]] = fitLL
  output =c(fitLL$par,"LL"=fitLL$value)
  parameter_sets = rbind(output,parameter_sets)
}

saveRDS(fits,"Submission Documents/Revisions/full optim object.rds")

names(parameter_sets) = c(colnames(grid_parms),"LL")
parameter_sets$LL = parameter_sets$LL*-1
parameter_sets=parameter_sets %>% mutate(parameter_set = row.names(.))
saveRDS(parameter_sets,"DATA/100_Fitted_Parameter_Sets.rds")

parameter_sets = readRDS("DATA/100_Fitted_Parameter_Sets.rds") 

#### Select all parameters sets within 2 -LL of the minimum 
min_ll = parameter_sets %>% filter(LL<=min(LL)+2)%>% 
  mutate(parameter_set2 = row.names(.))


#Plot Trajecories 
select_sets = min_ll
trajectories = data.frame()
age_distributions=data.frame()

dates1 = seq(from=as.Date('2017-07-02'),to=as.Date('2020-04-01'),by='week')

for(s in 1:nrow(select_sets)){
  baseline.txn.rate = exp(select_sets[s,1])
  b1 <- exp(select_sets[s,2])
  phi= (2*pi*(exp(select_sets[s,3]))) / (1+exp(select_sets[s,3]))
  report_infants <- 1 / (1 + exp(-select_sets[s,4])) 
  
  report_children <-1 / (1 + exp(-select_sets[s,5]))
  report_adults <- 1 / (1 + exp(-select_sets[s,6]))
  report_seniors60 <- 1 / (1 + exp(-select_sets[s,7]))
  report_seniors75 <- 1 / (1 + exp(-select_sets[s,8]))
  
  
  #run and plot fitted parameters
  results <- ode(y=yinit.vector, method = "ode45", times=fit_times,
                 func=MSIRS_immunization_dynamics,
                 parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi))
  
  t0 <- length(time_series_pre)
  al <- nrow(yinit)
  output <- tail(results,t0)
  St <- output[,-1]
  I1 <- St[,grep('I1', colnames(St))]
  I2 <- St[,grep('I2', colnames(St))]
  I3 <- St[,grep('I3', colnames(St))]
  I4 <- St[,grep('I4', colnames(St))]
  R1 <- St[,grep('R1', colnames(St))]
  R2 <- St[,grep('R2', colnames(St))]
  R3 <- St[,grep('R3', colnames(St))]
  R4 <- St[,grep('R4', colnames(St))]
  S1 <- St[,grep('S1', colnames(St))]
  S2 <- St[,grep('S2', colnames(St))]
  S3 <- St[,grep('S3', colnames(St))]
  S0 <- St[,grep('S0', colnames(St))]
  M0 <- St[,grep('M0', colnames(St))]
  Si<- St[,grep('Si', colnames(St))]
  Mn1<- St[,grep('Mn1', colnames(St))]
  Mn2<- St[,grep('Mn2', colnames(St))]
  Mv1<- St[,grep('Mv1', colnames(St))]
  Mv2<- St[,grep('Mv2', colnames(St))]
  N1<- St[,grep('N1', colnames(St))]
  N2<- St[,grep('N2', colnames(St))]
  Vs1<- St[,grep('Vs1', colnames(St))]
  Vs2<- St[,grep('Vs2', colnames(St))]
  
  
  contact2 = parmset$npi[105:length( parmset$npi)]
  intro2 = parmset$introductions[105:length( parmset$introductions)]  
  
  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0){
    beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2*contact2[t]
    lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}
  
  
  hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
  hosp2 <- hosp1 * 0.4
  hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors60, report_seniors75)
  
  H1 <- matrix(0, nrow = t0, ncol = al)
  for (i in 1:al) {
    H1[, i] <-
      hosp1[i] * parmset$RRHm * parmset$sigma3 * M0[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn1 * parmset$sigma3 * Mn1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn2 *  Mn2[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv1 * parmset$sigma3 * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv2 * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn1 * N1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn2 * N2[, i] * lambda1[, i] +
      hosp1[i] * S0[, i] * lambda1[, i] +
      hosp1[i] * Si[, i] * lambda1[, i] +
      hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs2[, i] * lambda1[, i]
  }
  
  H <- data.frame(H=rowSums(H1))
  plot(H$H)
  H$date = dates1
  H$parameter_set = s 
  trajectories = rbind(H,trajectories)
  
  H2 <- cbind(H1[,1],H1[,2],H1[,3],H1[,4],H1[,5],H1[,6], rowSums(H1[, 7:8]), rowSums(H1[, 9:12]), H1[, 13],H1[,14])
  age_dist2 <- c(colSums(H2)/sum(H2),s)
  age_distributions = rbind(age_dist2, age_distributions) 
  
  
}

#Plot time series 
plot_traj = ggplot()+
  theme_bw()+
  geom_line(data=trajectories,aes(x=date,y=H,group=parameter_set,color=as.factor(parameter_set)))+
  geom_line(aes(x=dates1,y=time_series_pre),color="black",linewidth=1)
plot_traj

#plot age distribution 
ages = c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1-4yrs","5-59yrs","60-74yrs","75+yrs")
names(age_distributions) = c(ages,"set")

age_dist_obs = data.frame(age=ages, proportion=age_dist)
age_long = age_distributions %>% 
  pivot_longer(cols=c(`<2m`:`75+yrs`),names_to="age",values_to="proportion") %>% 
  mutate(age=factor(age,levels=ages))

plot_age=ggplot()+
  theme_bw()+
  geom_point(data=age_long,aes(x=age,y=proportion),size=3,shape=1,color="grey")+
  geom_point(data=age_dist_obs,aes(x=age,y=proportion),size=3,color="red")
plot_age




# Fit_NPI period ----------------------------------------------------------
start_lower = c(-2,-2,-2)
start_upper = c(2,2,2)
conf_int = cbind(start_lower, start_upper)
set.seed(123)
h=1 #only using one set of starting values for each parameter set 
lhs<-maximinLHS(h,3)

npi_grid <- cbind(
  npi1 = lhs[,1]*(conf_int[1,2]-conf_int[1,1])+conf_int[1,1],
  npi2 = lhs[,2]*(conf_int[2,2]-conf_int[2,1])+conf_int[2,1],
  npi3 = lhs[,3]*(conf_int[3,2]-conf_int[3,1])+conf_int[3,1])


fit_times2 <- seq(1, length(time_series_full) + 104, by = 1)
npi_sets = data.frame()

for(n in 1:1){
  for(m in 1:nrow(min_ll)){
    fitpand <-  function(parameters,dat) {
      
      npi1 <-1 / (1 + exp(-parameters[1])) 
      npi2 <- 1 / (1 + exp(-parameters[2])) 
      npi3 <- 1 / (1 + exp(-parameters[3])) 
      
      
      baseline.txn.rate = exp(min_ll[m,1])
      b1 <- exp(min_ll[m,2])
      phi= (2*pi*(exp(min_ll[m,3]))) / (1+exp(min_ll[m,3]))
      report_infants <- 1 / (1 + exp(-min_ll[m,4])) 
      report_children <-1 / (1 + exp(-min_ll[m,5]))
      report_adults <- 1 / (1 + exp(-min_ll[m,6]))
      report_seniors60 <- 1 / (1 + exp(-min_ll[m,7]))
      report_seniors75 <-1 / (1 + exp(-min_ll[m,8]))
      theta <- 1 / (1 + exp(-min_ll[m,9]))
      
      
      introductions = data.frame(intros=c(rep(parmset$seed,248),rep(0,44),rep(NA,24),rep(parmset$seed,200))) %>% 
        mutate(intros = na_interpolation(intros, method="linear"))
      introductions = introductions$intros
      parmset$introductions = introductions
      
      #reductions in contacts during the pandemic 
      npi = data.frame(npis=c(rep(1,248),rep(npi1,52),rep(NA,22),rep(npi2,13),rep(npi3,26),rep(NA,13),rep(1,200)))%>% 
        mutate(npis= na_interpolation(npis, method="linear"))
      npi=npi$npis
      parmset$npi = npi
      
      # Simulate the model with initial conditions and timesteps defined above, and parameter values from function call
      results <- ode(y=yinit.vector, method = "ode45", times=fit_times2,
                     func=MSIRS_immunization_dynamics,
                     parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi))
      
      t0 <- length(time_series_full)
      al <- nrow(yinit)
      output <- tail(results,t0)
      St <- output[,-1]
      I1 <- St[,grep('I1', colnames(St))]
      I2 <- St[,grep('I2', colnames(St))]
      I3 <- St[,grep('I3', colnames(St))]
      I4 <- St[,grep('I4', colnames(St))]
      R1 <- St[,grep('R1', colnames(St))]
      R2 <- St[,grep('R2', colnames(St))]
      R3 <- St[,grep('R3', colnames(St))]
      R4 <- St[,grep('R4', colnames(St))]
      S1 <- St[,grep('S1', colnames(St))]
      S2 <- St[,grep('S2', colnames(St))]
      S3 <- St[,grep('S3', colnames(St))]
      S0 <- St[,grep('S0', colnames(St))]
      M0 <- St[,grep('M0', colnames(St))]
      Si<- St[,grep('Si', colnames(St))]
      Mn1<- St[,grep('Mn1', colnames(St))]
      Mn2<- St[,grep('Mn2', colnames(St))]
      Mv1<- St[,grep('Mv1', colnames(St))]
      Mv2<- St[,grep('Mv2', colnames(St))]
      N1<- St[,grep('N1', colnames(St))]
      N2<- St[,grep('N2', colnames(St))]
      Vs1<- St[,grep('Vs1', colnames(St))]
      Vs2<- St[,grep('Vs2', colnames(St))]
      
      contact2 = parmset$npi[105:length(parmset$npi)]
      intro2 = parmset$introductions[105:length(parmset$introductions)]  
      
      lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
      for (t in 1:t0){
        beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2*contact2[t]
        lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}
      
      hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
      hosp2 <- hosp1 * 0.4
      hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors60, report_seniors75)
      
      H1 <- matrix(0, nrow = t0, ncol = al)
      for (i in 1:al) {
        H1[, i] <-
          hosp1[i] * parmset$RRHm * parmset$sigmaM * M0[, i] * lambda1[, i] +
          hosp1[i] * parmset$RRHm * parmset$RRHn1 * parmset$sigmaM * Mn1[, i] * lambda1[, i] +
          hosp1[i] * parmset$RRHm * parmset$RRHn2 *  Mn2[, i] * lambda1[, i] +
          hosp1[i] * parmset$RRHm * parmset$RRHv1 * parmset$sigmaM * Mv1[, i] * lambda1[, i] +
          hosp1[i] * parmset$RRHm * parmset$RRHv2 * Mv1[, i] * lambda1[, i] +
          hosp1[i] * parmset$RRHn1 * N1[, i] * lambda1[, i] +
          hosp1[i] * parmset$RRHn2 * N2[, i] * lambda1[, i] +
          hosp1[i] * S0[, i] * lambda1[, i] +
          hosp1[i] * Si[, i] * lambda1[, i] +
          hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
          hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
          hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
          hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs1[, i] * lambda1[, i] +
          hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs2[, i] * lambda1[, i]
      }
      
      H <- rowSums(H1)
      
      LL <- sum(dnbinom(x = dat, mu =H, size=1/theta, log = TRUE))
      return(LL)
      
    }
    
    pandLL <- optim(par = npi_grid[n,],
                    fn = fitpand,  # the distance function to optimize
                    dat = time_series_full,  # the dataset to fit to (dpois function)
                    control = list(fnscale=-1,maxit=5000))#
    
    output =c(pandLL$par,"LL"=pandLL$value, npi_set=n, parameter_set = m)
    npi_sets = rbind(output,npi_sets)
  }}

names(npi_sets) = c("npi1","npi2","npi3","LL","npi_set","parameter_set")
npi_sets = npi_sets %>% mutate(LL = LL*-1)
saveRDS(npi_sets,"DATA/fitted_NPIs.rds")
npi_fits = readRDS("DATA/fitted_NPIs.rds") 

min_ll$parameter_set2 = as.character(min_ll$parameter_set2)
npi_fits$parameter_set = as.character(npi_fits$parameter_set)

combo = min_ll %>% 
  right_join(npi_fits, by=c("parameter_set2"="parameter_set")) %>% 
  select(beta,b1,phi,RI,RC,RA,RS60,RS75,npi1,npi2,npi3, parameter_set2,npi_set)


## Plot full period with the 8 parameter sets 
fit_times3 <- seq(1, nrow(timeseries) + 104, by = 1)
trajectories=data.frame()
save_npi=data.frame()
for(n in 1:nrow(combo)){
  
  baseline.txn.rate = exp(combo[n,1])
  b1 <- exp(combo[n,2])
  phi= (2*pi*(exp(combo[n,3]))) / (1+exp(combo[n,3]))
  report_infants <- 1 / (1 + exp(-combo[n,4])) 
  report_children <-1 / (1 + exp(-combo[n,5]))
  report_adults <- 1 / (1 + exp(-combo[n,6]))
  report_seniors60 <- 1 / (1 + exp(-combo[n,7]))
  report_seniors75 <-1 / (1 + exp(-combo[n,8]))
  npi1 = 1 / (1 + exp(-combo[n,9]))
  npi2 = 1 / (1 + exp(-combo[n,10]))
  npi3 = 1 / (1 + exp(-combo[n,11]))
  
  
  introductions = data.frame(intros=c(rep(parmset$seed,248),rep(0,44),rep(NA,24),rep(parmset$seed,500))) %>% 
    mutate(intros = na_interpolation(intros, method="linear"))
  introductions = introductions$intros
  
  parmset$introductions = introductions
  
  #reductions in contacts during the pandemic 
  npi = data.frame(npis=c(rep(1,248),rep(npi1,52),rep(NA,22),rep(npi2,13),rep(npi3,26),rep(NA,13),rep(1,500)))%>% 
    mutate(npis= na_interpolation(npis, method="linear"))
  npi=npi$npis
  parmset$npi = npi
  
  # Simulate the model with initial conditions and timesteps defined above, and parameter values from function call
  results <- ode(y=yinit.vector, method = "ode45", times=fit_times3,
                 func=MSIRS_immunization_dynamics,
                 parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi))
  
  t0 <- nrow(timeseries)
  al <- nrow(yinit)
  output <- tail(results,t0)
  St <- output[,-1]
  I1 <- St[,grep('I1', colnames(St))]
  I2 <- St[,grep('I2', colnames(St))]
  I3 <- St[,grep('I3', colnames(St))]
  I4 <- St[,grep('I4', colnames(St))]
  R1 <- St[,grep('R1', colnames(St))]
  R2 <- St[,grep('R2', colnames(St))]
  R3 <- St[,grep('R3', colnames(St))]
  R4 <- St[,grep('R4', colnames(St))]
  S1 <- St[,grep('S1', colnames(St))]
  S2 <- St[,grep('S2', colnames(St))]
  S3 <- St[,grep('S3', colnames(St))]
  S0 <- St[,grep('S0', colnames(St))]
  M0 <- St[,grep('M0', colnames(St))]
  Si<- St[,grep('Si', colnames(St))]
  Mn1<- St[,grep('Mn1', colnames(St))]
  Mn2<- St[,grep('Mn2', colnames(St))]
  Mv1<- St[,grep('Mv1', colnames(St))]
  Mv2<- St[,grep('Mv2', colnames(St))]
  N1<- St[,grep('N1', colnames(St))]
  N2<- St[,grep('N2', colnames(St))]
  Vs1<- St[,grep('Vs1', colnames(St))]
  Vs2<- St[,grep('Vs2', colnames(St))]
  
  
  contact2 = parmset$npi[105:length(parmset$npi)]
  intro2 = parmset$introductions[105:length(parmset$introductions)]  
  
  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0){
    beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2*contact2[t]
    lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}
  
  hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
  hosp2 <- hosp1 * 0.4
  hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors60, report_seniors75)
  
  H1 <- matrix(0, nrow = t0, ncol = al)
  for (i in 1:al) {
    H1[, i] <-
      hosp1[i] * parmset$RRHm * parmset$sigmaM * M0[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn1 * parmset$sigmaM * Mn1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn2 *  Mn2[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv1 * parmset$sigmaM * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv2 * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn1 * N1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn2 * N2[, i] * lambda1[, i] +
      hosp1[i] * S0[, i] * lambda1[, i] +
      hosp1[i] * Si[, i] * lambda1[, i] +
      hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs2[, i] * lambda1[, i]
  }
  
  H2 = cbind(rowSums(H1[,1:3]),
             rowSums(H1[,4:6]),
             rowSums(H1[,7:8]),
             rowSums(H1[,9:12]),
             H1[,13], H1[,14])
  
  H = data.frame(H2)
  names(H)=c("<6m","6-11m","1-4yrs","5-59yrs","60-74yrs","75+yrs")
  H$All = rowSums(H2)
  H$date = timeseries$date
  
  H = H %>% 
    mutate(sample=n) %>% 
    pivot_longer(cols=c("<6m":"All"),names_to="Age",values_to="value")
  
  
  trajectories=rbind(H,trajectories) %>% filter(!is.na(date))
  save_cont = data.frame(npi=contact2[1:nrow(timeseries)]) %>% mutate(date=timeseries$date, sample=n)
  save_npi = rbind(save_cont,save_npi)
}


final_plot = ggplot()+
  theme_bw()+
  geom_line(data=trajectories %>% filter(Age=="All"),aes(x=date,y=value/22.67, group=sample,color="Model Trajectory"))+
  geom_point(data=timeseries ,aes(x=date, y=hosp_rate,color="Observed Data"),shape=1)+
  geom_line(data=save_npi, aes(x=date, y=npi*100/22.67,group=sample,color="Contact Intensity"))+
  scale_color_manual(name=NULL,values=c("blue","grey70","black"))+
  scale_y_continuous(
    name = "RSV Hospitalization Rate", 
    sec.axis = sec_axis(~.*22.67, name = "% of Contacts") # Secondary axis
  )+
  geom_text(aes(x=as.Date("2020-12-01"),y=3.35,label="c1"))+
  geom_text(aes(x=as.Date("2021-10-15"),y=4,label="c2"))+
  geom_text(aes(x=as.Date("2022-2-15"),y=3,label="c3"))+
  labs(x=NULL)+
  theme(axis.text=element_text(size=12),
        axis.title = element_text(size=15),
        legend.text = element_text(size=12))
final_plot




# Slices - pre-pandemic parameters-----------------------------------------------------------------
slices=data.frame()
for(m in 1:nrow(min_ll)){
  slice_beta = seq(from=min_ll[m,1]*.8, to=min_ll[m,1]*1.2,by=0.005)
  slice_beta = cbind("beta"=slice_beta, min_ll[m,2:8],slice="beta")
  
  slice_b1 = seq(from=min_ll[m,2]*1.2, to=min_ll[m,2]*0.8, by=0.005)
  slice_b1 = cbind("b1"=slice_b1, min_ll[m,c(1,3:8)],slice="b1")
  
  slice_phi = seq(from=min_ll[m,3]*.8, to=min_ll[m,3]*1.2,by=0.001)
  slice_phi = cbind("phi"=slice_phi, min_ll[m,c(1:2,4:8)],slice="phi")
  
  slice_RI = seq(from=min_ll[m,4]*1.2, to=min_ll[m,4]*0.8,by=0.01)
  slice_RI = cbind("RI"=slice_RI, min_ll[m,c(1:3,5:8)],slice="RI")
  
  slice_RC = seq(from=min_ll[m,5]*1.2, to=min_ll[m,5]*0.8,by=0.01)
  slice_RC = cbind("RC"=slice_RC, min_ll[m,c(1:4,6:8)],slice="RC")
  
  slice_RA = seq(from=min_ll[m,6]*1.2, to=min_ll[m,6]*0.8,by=0.01)
  slice_RA = cbind("RA"=slice_RA, min_ll[m,c(1:5,7:8)],slice="RA")
  
  slice_RS60 = seq(from=min_ll[m,7]*1.2, to=min_ll[m,7]*0.8,by=0.01)
  slice_RS60 = cbind("RS60"=slice_RS60, min_ll[m,c(1:6,8)],slice="RS60")
  
  slice_RS75 = seq(from=min_ll[m,8]*1.2, to=min_ll[m,8]*0.8,by=0.01)
  slice_RS75 = cbind("RS75"=slice_RS75, min_ll[m,c(1:7)],slice="RS75")
  
  slice_new = bind_rows(slice_beta,slice_b1,slice_phi,slice_RI,slice_RC,slice_RA,slice_RS60,slice_RS75) %>% mutate(parameter_set2=m, theta=min_ll[m,9])
  slices=rbind(slices,slice_new)
}

slices_ll = data.frame()

for(s in 1:nrow(slices)){
  
  protrans <- slices[s,1] # parameter for baseline transmission rate
  baseline.txn.rate = exp(protrans) 
  amp <-  slices[s,2]# parameter for seasonal amplitude
  b1 <-  exp(amp)
  trans <- slices[s,3] # parameter for seasonal phase
  phi <-  (2*pi*(exp(trans))) / (1+exp(trans)) # transform to between 0 and 2pi
  #Age-specific reporting fractions 
  report_infants <- 1 / (1 + exp(-slices[s,4])) 
  report_children <- 1 / (1 + exp(-slices[s,5]))
  report_adults <- 1 / (1 + exp(-slices[s,6]))
  report_seniors60 <- 1 / (1 + exp(-slices[s,7]))
  report_seniors75 <- 1 / (1 + exp(-slices[s,8]))
  theta <- 1 / (1 + exp(-slices[s,11]))
  
  results <- ode(y=yinit.vector, method = "ode45", times=fit_times,
                 func=MSIRS_immunization_dynamics,
                 parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi))
  
  
  t0 <- length(time_series_pre)
  al <- nrow(yinit)
  output <- tail(results,t0)
  St <- output[,-1]
  I1 <- St[,grep('I1', colnames(St))]
  I2 <- St[,grep('I2', colnames(St))]
  I3 <- St[,grep('I3', colnames(St))]
  I4 <- St[,grep('I4', colnames(St))]
  R1 <- St[,grep('R1', colnames(St))]
  R2 <- St[,grep('R2', colnames(St))]
  R3 <- St[,grep('R3', colnames(St))]
  R4 <- St[,grep('R4', colnames(St))]
  S1 <- St[,grep('S1', colnames(St))]
  S2 <- St[,grep('S2', colnames(St))]
  S3 <- St[,grep('S3', colnames(St))]
  S0 <- St[,grep('S0', colnames(St))]
  M0 <- St[,grep('M0', colnames(St))]
  Si<- St[,grep('Si', colnames(St))]
  Mn1<- St[,grep('Mn1', colnames(St))]
  Mn2<- St[,grep('Mn2', colnames(St))]
  Mv1<- St[,grep('Mv1', colnames(St))]
  Mv2<- St[,grep('Mv2', colnames(St))]
  N1<- St[,grep('N1', colnames(St))]
  N2<- St[,grep('N2', colnames(St))]
  Vs1<- St[,grep('Vs1', colnames(St))]
  Vs2<- St[,grep('Vs2', colnames(St))]
  
  
  contact2 = parmset$npi[105:length( parmset$npi)]
  intro2 = parmset$introductions[105:length( parmset$introductions)]  
  
  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0){
    beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2*contact2[t]
    lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}
  
  
  
  hosp1 <- c(report_infants, report_infants*0.59, report_infants*0.33, report_infants*0.2, report_infants*0.15, report_infants*0.15, rep(report_children, 2), rep(.001, 6))
  hosp2 <- hosp1 * 0.4
  hosp3 <- c(rep(0.001, 8), rep(report_adults, 4), report_seniors60, report_seniors75)
  
  H1 <- matrix(0, nrow = t0, ncol = al)
  for (i in 1:al) {
    H1[, i] <-
      hosp1[i] * parmset$RRHm * parmset$sigmaM * M0[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn1 * parmset$sigmaM * Mn1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn2 *  Mn2[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv1 * parmset$sigmaM * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv2 * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn1 * N1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn2 * N2[, i] * lambda1[, i] +
      hosp1[i] * S0[, i] * lambda1[, i] +
      hosp1[i] * Si[, i] * lambda1[, i] +
      hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs2[, i] * lambda1[, i]
  }
  
  H <- rowSums(H1)
  H2 <- cbind(H1[,1],H1[,2],H1[,3],H1[,4],H1[,5],H1[,6], rowSums(H1[, 7:8]), rowSums(H1[, 9:12]), H1[, 13],H1[,14])
  age_dist2 <- colSums(H2)
  
  LLall <- sum(dnbinom(x = time_series_pre, mu =H, size=1/theta, log = TRUE))
  LLmulti <- dmultinom(x = age_dist2, prob = age_dist, log = TRUE)
  
  LL <- (LLall + LLmulti)*-1
  
  output = cbind(slices[s,],LL=LL, parameters=m)
  slices_ll = rbind(slices_ll, output)
}


saveRDS(slices_ll,"DATA/slices_core.rds")
slices_ll = readRDS("DATA/slices_core.rds")


slices_long = slices_ll %>% 
  select(-parameters) %>% 
  pivot_longer(cols=c(beta:RS75),names_to="parameter",values_to="value")


set1_plot= ggplot(slices_long %>% filter(parameter_set2==1) %>% filter(slice==parameter, LL<=min(LL)+4))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=4)+
  geom_hline(aes(yintercept = min(LL)+2))+
  geom_hline(aes(yintercept = min(LL)+2.5))
set1_plot

set2_plot= ggplot(slices_long %>% filter(parameter_set2==2) %>% filter(slice==parameter, LL<=min(LL)+4))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=4)+
  geom_hline(aes(yintercept = min(LL)+2))
set2_plot

set3_plot= ggplot(slices_long %>% filter(parameter_set2==3) %>% filter(slice==parameter, LL<=min(LL)+4))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=4)+
  geom_hline(aes(yintercept = min(LL)+2))
set3_plot

set4_plot= ggplot(slices_long %>% filter(parameter_set2==4) %>% filter(slice==parameter, LL<=min(LL)+4))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=4)+
  geom_hline(aes(yintercept = min(LL)+2))
set4_plot

set5_plot= ggplot(slices_long %>% filter(parameter_set2==5) %>% filter(slice==parameter, LL<=min(LL)+4))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=4)+
  geom_hline(aes(yintercept = min(LL)+2))
set5_plot

set6_plot= ggplot(slices_long %>% filter(parameter_set2==6) %>% filter(slice==parameter, LL<=min(LL)+4))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=4)+
  geom_hline(aes(yintercept = min(LL)+2))
set6_plot

set7_plot= ggplot(slices_long %>% filter(parameter_set2==7) %>% filter(slice==parameter, LL<=min(LL)+4))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=4)+
  geom_hline(aes(yintercept = min(LL)+2))
set7_plot
set8_plot= ggplot(slices_long %>% filter(parameter_set2==8) %>% filter(slice==parameter, LL<=min(LL)+4))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=4)+
  geom_hline(aes(yintercept = min(LL)+2))
set8_plot



# Get all values within -4LL of the minimum
get_boundaries = slices_long %>% 
  group_by(parameter_set2) %>% 
  filter(LL<=min(LL)+4) %>% 
  ungroup() %>% 
  group_by(parameter_set2, parameter) %>% 
  summarize(min = min(value),
            max=max(value))


# Slices - NPI periods  ---------------------------------------------------
npi_slices=data.frame()
for(n in 1:nrow(npi_fits)){
  npi_fits = npi_fits %>% arrange(parameter_set)
  slice_npi1 = seq(from=npi_fits[n,1]*.5, to=npi_fits[n,1]*1.5,by=0.01)
  slice_npi1 = cbind("npi1"=slice_npi1, npi_fits[n,2:3],slice="npi1")
  
  slice_npi2 = seq(from=npi_fits[n,2]*.5, to=npi_fits[n,2]*1.5,by=0.01)
  slice_npi2 = cbind("npi2"=slice_npi2, npi_fits[n,c(1,3)],slice="npi2")
  
  slice_npi3 = seq(from=npi_fits[n,3]*.5, to=npi_fits[n,3]*1.5,by=0.01)
  slice_npi3 = cbind("npi3"=slice_npi3, npi_fits[n,1:2],slice="npi3")
  
  npi_slices2 = bind_rows(slice_npi1,slice_npi2,slice_npi3) %>% mutate(parameter_set2=n)
  npi_slices = rbind(npi_slices,npi_slices2)
}

npi_slices = npi_slices %>% 
  mutate(parameter_set2 = as.character(parameter_set2)) %>% 
  left_join(min_ll %>% select(-LL,-parameter_set),by=c("parameter_set2"))



slices_npi_ll = data.frame()

for(n in 1:nrow(npi_slices)){
  
  npi1 <-1 / (1 + exp(-npi_slices[n,1]))
  npi2 <- 1 / (1 + exp(-npi_slices[n,2]))
  npi3 <- 1 / (1 + exp(-npi_slices[n,3]))
  
  baseline.txn.rate = exp(npi_slices[n,6])
  b1 <- exp(npi_slices[n,7])
  phi= (2*pi*(exp(npi_slices[n,8]))) / (1+exp(npi_slices[n,8]))
  report_infants <- 1 / (1 + exp(-npi_slices[n,9])) 
  report_children <-1 / (1 + exp(-npi_slices[n,10]))
  report_adults <- 1 / (1 + exp(-npi_slices[n,11]))
  report_seniors60 <- 1 / (1 + exp(-npi_slices[n,12]))
  report_seniors75 <-1 / (1 + exp(-npi_slices[n,13]))
  theta <-1 / (1 + exp(-npi_slices[n,14]))
  
  
  introductions = data.frame(intros=c(rep(parmset$seed,248),rep(0,44),rep(NA,24),rep(parmset$seed,200))) %>% 
    mutate(intros = na_interpolation(intros, method="linear"))
  introductions = introductions$intros
  parmset$introductions = introductions
  
  #reductions in contacts during the pandemic 
  npi = data.frame(npis=c(rep(1,248),rep(npi1,52),rep(NA,22),rep(npi2,13),rep(npi3,26),rep(NA,13),rep(1,200)))%>% 
    mutate(npis= na_interpolation(npis, method="linear"))
  npi=npi$npis
  parmset$npi = npi
  
  # Simulate the model with initial conditions and timesteps defined above, and parameter values from function call
  results <- ode(y=yinit.vector, method = "ode45", times=fit_times2,
                 func=MSIRS_immunization_dynamics,
                 parms=c(parmset,baseline.txn.rate=baseline.txn.rate,b1=b1,phi=phi))
  
  t0 <- length(time_series_full)
  al <- nrow(yinit)
  output <- tail(results,t0)
  St <- output[,-1]
  I1 <- St[,grep('I1', colnames(St))]
  I2 <- St[,grep('I2', colnames(St))]
  I3 <- St[,grep('I3', colnames(St))]
  I4 <- St[,grep('I4', colnames(St))]
  R1 <- St[,grep('R1', colnames(St))]
  R2 <- St[,grep('R2', colnames(St))]
  R3 <- St[,grep('R3', colnames(St))]
  R4 <- St[,grep('R4', colnames(St))]
  S1 <- St[,grep('S1', colnames(St))]
  S2 <- St[,grep('S2', colnames(St))]
  S3 <- St[,grep('S3', colnames(St))]
  S0 <- St[,grep('S0', colnames(St))]
  M0 <- St[,grep('M0', colnames(St))]
  Si<- St[,grep('Si', colnames(St))]
  Mn1<- St[,grep('Mn1', colnames(St))]
  Mn2<- St[,grep('Mn2', colnames(St))]
  Mv1<- St[,grep('Mv1', colnames(St))]
  Mv2<- St[,grep('Mv2', colnames(St))]
  N1<- St[,grep('N1', colnames(St))]
  N2<- St[,grep('N2', colnames(St))]
  Vs1<- St[,grep('Vs1', colnames(St))]
  Vs2<- St[,grep('Vs2', colnames(St))]
  
  
  
  contact2 = parmset$npi[105:length( parmset$npi)]
  intro2 = parmset$introductions[105:length( parmset$introductions)]  
  
  lambda1=matrix(0,nrow=t0,ncol=al)#Force of infection
  for (t in 1:t0){
    beta <-  baseline.txn.rate/(parmset$dur.days1/7)/(sum(yinit)^(1-parmset$q))*parmset$c2*contact2[t]
    lambda1[t,] <- as.vector((1+b1*cos(2*pi*(t-phi*52.1775)/52.1775))*((I1[t,]+parmset$rho1*I2[t,]+parmset$rho2*I3[t,]+parmset$rho2*I4[t,]+intro2[t])%*%beta)/sum(St[t,]))}
  
  H1 <- matrix(0, nrow = t0, ncol = al)
  for (i in 1:al) {
    H1[, i] <-
      hosp1[i] * parmset$RRHm * parmset$sigmaM * M0[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn1 * parmset$sigmaM * Mn1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHn2 *  Mn2[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv1 * parmset$sigmaM * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHm * parmset$RRHv2 * Mv1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn1 * N1[, i] * lambda1[, i] +
      hosp1[i] * parmset$RRHn2 * N2[, i] * lambda1[, i] +
      hosp1[i] * S0[, i] * lambda1[, i] +
      hosp1[i] * Si[, i] * lambda1[, i] +
      hosp2[i] * parmset$sigma1 * S1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma2 * S2[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * S3[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs1[, i] * lambda1[, i] +
      hosp3[i] * parmset$sigma3 * parmset$RRHs * Vs2[, i] * lambda1[, i]
  }
  
  H <- rowSums(H1)#[145:325]
  LL <- sum(dnbinom(x = time_series_full, mu =H, size=1/theta, log = TRUE))
  
  output = cbind(npi_slices[n,],LL=LL)
  slices_npi_ll = rbind(slices_npi_ll, output)
}

slices_npi_ll$LL = slices_npi_ll$LL*-1
saveRDS(slices_npi_ll,"DATA/slices_npi.rds")
slices_npi_ll = readRDS("DATA/slices_npi.rds")

slices_npi_long = slices_npi_ll %>%
  select(npi1,npi2,npi3,slice,parameter_set2,LL) %>% 
  pivot_longer(cols=c(npi1:npi3),names_to="parameter",values_to="value") %>% 
  filter(slice==parameter)


npi_set1= ggplot(slices_npi_long %>% filter(parameter_set2==1) %>% filter(LL<=min(LL)+3))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=3)+
  geom_hline(aes(yintercept = min(LL)+2))
npi_set1

npi_set2= ggplot(slices_npi_long %>% filter(parameter_set2==2) %>% filter(LL<=min(LL)+3))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=3)+
  geom_hline(aes(yintercept = min(LL)+2))
npi_set2

npi_set3= ggplot(slices_npi_long %>% filter(parameter_set2==3) %>% filter(LL<=min(LL)+3))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=3)+
  geom_hline(aes(yintercept = min(LL)+2))
npi_set3

npi_set4= ggplot(slices_npi_long %>% filter(parameter_set2==4) %>% filter(LL<=min(LL)+3))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=3)+
  geom_hline(aes(yintercept = min(LL)+2))
npi_set4

npi_set5= ggplot(slices_npi_long %>% filter(parameter_set2==5) %>% filter(LL<=min(LL)+3))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=3)+
  geom_hline(aes(yintercept = min(LL)+2))
npi_set5

npi_set6= ggplot(slices_npi_long %>% filter(parameter_set2==6) %>% filter(LL<=min(LL)+3))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=3)+
  geom_hline(aes(yintercept = min(LL)+2))
npi_set6

npi_set7= ggplot(slices_npi_long %>% filter(parameter_set2==7) %>% filter(LL<=min(LL)+3))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=3)+
  geom_hline(aes(yintercept = min(LL)+2))
npi_set7

npi_set8= ggplot(slices_npi_long %>% filter(parameter_set2==8) %>% filter(LL<=min(LL)+3))+
  theme_bw()+
  geom_point(aes(x=value, y=LL))+
  facet_wrap(~slice,scales="free",ncol=3)+
  geom_hline(aes(yintercept = min(LL)+2))
npi_set8

#Get all values within -4LL of the minimum
npi_boundaries = slices_npi_long %>% 
  group_by(parameter_set2) %>% 
  filter(LL<=min(LL)+4) %>% 
  ungroup() %>% 
  group_by(parameter_set2, parameter) %>% 
  summarize(min = min(value),
            max=max(value)) %>% 
  mutate(parameter_set2 = as.integer(parameter_set2))


# Combine Slices to get 8 parameter Sets  ---------------------------------

boundaries = rbind(get_boundaries, npi_boundaries) %>% 
  mutate(parameter = factor(parameter, levels=c("beta","b1","phi","RI","RC","RA","RS60","RS75","npi1","npi2","npi3"))) %>% 
  mutate(flag = ifelse(parameter=="phi" & min<1,1,0)) %>% #drop the phase that's in a different cycle 
  filter(flag==0)

#get the global minimum and maximum from each set 
boundaries2 = boundaries %>% 
  group_by(parameter) %>% 
  summarize(min=min(min),max=max(max)) %>% 
  mutate(parameter = factor(parameter, levels=c("beta","b1","phi","RI","RC","RA","RS60","RS75","npi1","npi2","npi3")))

boundary_plot=ggplot(boundaries)+
  theme_bw()+
  geom_errorbar(aes(xmin=min,xmax=max,y=parameter_set2),width=0,linewidth=2)+
  facet_wrap(~parameter,scales="free")
boundary_plot
# LHS to get 100 combinations for each set  -------------------------------

set.seed(123)

h=100
lhs<-maximinLHS(h,11)

apply_lhs = function(conf_int){
  grid_parms <- cbind(
    beta = lhs[,1]*(conf_int[1,2]-conf_int[1,1])+conf_int[1,1],
    b1 = lhs[,2]*(conf_int[2,2]-conf_int[2,1])+conf_int[2,1],
    phi = lhs[,3]*(conf_int[3,2]-conf_int[3,1])+conf_int[3,1],
    RI = lhs[,4]*(conf_int[4,2]-conf_int[4,1])+conf_int[4,1],
    RC = lhs[,5]*(conf_int[5,2]-conf_int[5,1])+conf_int[5,1],
    RA = lhs[,6]*(conf_int[6,2]-conf_int[6,1])+conf_int[6,1],
    RS60 = lhs[,7]*(conf_int[7,2]-conf_int[7,1])+conf_int[7,1],
    RS75 = lhs[,8]*(conf_int[8,2]-conf_int[8,1])+conf_int[8,1],
    npi1 = lhs[,9]*(conf_int[9,2]-conf_int[9,1])+conf_int[9,1],
    npi2 = lhs[,10]*(conf_int[10,2]-conf_int[10,1])+conf_int[10,1],
    npi3 = lhs[,11]*(conf_int[11,2]-conf_int[11,1])+conf_int[11,1])
  return(grid_parms)
}


conf_int = as.matrix(boundaries2[1:11,2:3])
samples = apply_lhs(conf_int)

grid_parms = data.frame(samples) %>% mutate(sample=row.names(.))
saveRDS(grid_parms, "DATA/fitted_parameters.rds")
