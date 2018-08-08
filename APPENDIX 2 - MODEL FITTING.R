#APPENDIX 2 - MODEL FITTING
###########################

#CLEAR WORKSPACE
rm(list=ls())
setwd("C:/Users/Lausdeus Chiegboka/Documents/ox docs/PLACEMENT PROJECT/MAIN/data")
getwd()

data<-read.csv('South-West Malaria data from 2004 to 2011.csv') #Data attached as annex A to this appendix

par(mfrow=c(1,1))
plot(data$Cases, type="o",lwd=3, col="red", ylim=c(0,2500000))
plot(data$CumCases)


#LOAD REQUIRED LIBRARIES
library(deSolve)
library(gtools)

malaria <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    mub=((1/54)+(2.26/100))               # birth rate (1/life expectancy)
    mud=(1/54)               # death rate (1/life expectancy)
    muc=10.3/100                 # Severe Malaria case fatality rate
    omega=(1/5)                # rate of loss of immunity = 1/(average duration of immunity)
    alpha=(1/(7.5/365))        # rate of movement from liver stage to asexual blood stage = 1/(average prepatent period)
    beta=(1/(2.5/365))         # Rate of movement  from asexual blood stage to gametocyte stage = 1/(latent period - prepatent period)
    deltaone=exp(logdeltaone) #(1/(5/365))       # Rate of development of severe disease = 1/duration of clinical disease before onset of severe malaria
    deltatwo=exp(logdeltatwo) #(1/(3/365))       # Rate of treatment = 1/time lapse in infectious stage before seeking treatment
    deltathree=exp(logdeltathree) #(1/(2.5/52))      # Rate of symptom resolution without treatment = 1/duration of clinical disease without treatment
    vee=exp(logvee) #(8.9/100)             # Proportion of clinical malaria that receive treatment = treatment coverage of ACTs
    chi=exp(logchi) #(2/100)                # Proportion of clinical malaria cases that become severe
    gamma=exp(loggamma) #(1/(2/365))          # Rate of treatment in severe malaria = 1/time lapse before seeking treatment of severe disease
    sigmaone=(1/(32/24/365))   # Rate of recovery following treatment = 1/parasite clearance time
    sigmatwo=exp(logsigmatwo)      # Rate of recrudescence = 1/time to recrudescence
    rho=(4.2/100)                 # Proportion of treated cases with treatment failure (recrudescence)
    kappa=exp(logkappa) #(90/100)                # Proportion of clinically immune that develop latent infection
    eps=(1/(18/12))                    # Rate of natural recovery from infectiousness (among the clinically immune)
    bee=0.5*365  #exp(logbee)   #                 # Mosquito biting rate per man/year
    
    
    amp=inv.logit(logitamp)                         # relative amplitude of seasonal forcing
    phi=inv.logit(logitphi)                         # YEAR of peak in seasonal forcing
    
    llin_eff = 0.5      # efficacy of llins in reducing disease incidence
    
    
    #mosquito parameters
    mubm=(1/(12.24/365))      # Mosquito birth rate = 1/life expectancy of mosquitoes
    mum=(1/(12.24/365))        # Mosquito death rate = 1/life expectancy of mosquitoes
    mgamma=(1/(11/365))        # Rate of movement from exposed to infectiousness in mosquitoes = 1/latent period
    ai=(5/100) #inv.logit(logitai)                    # Per bite probability of transmission from infected mosquito to human
    aim=(47/100) #inv.logit(logitaim)                   # Per bite probability of transmission from infectious human to mosquito 
    
    # DEFINE VARIABLES
    P <- (S+L+B+C+D+Tr+R+Lr+Br+G)
    Pm <- (Sm+Em+Im)
    
    seas<-1+amp*cos(2*pi*(t-phi))
    
    #starting from 1% coverage in 2004
    llin_cov=c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07,0.2,0.2) #to relect published estimates of llin coverage
    llin<-NULL
    for(jj in 1:length(llin_cov)){
      llin <- c(llin,1-llin_eff*llin_cov[jj])
    }
    lam<-NULL
    for(hh in 1:length(llin)){
      lam<-c(lam,llin[hh]*seas*bee*ai*Im/Pm)
    }
    mlam<-seas*bee*aim*(C+D+G)/P
    
    for(dd in 1:length(lam)){
      births<-(mub*P+muc*D)
      # RATE OF CHANGE - HUMAN
      dS <-births+omega*(1-kappa)*R-mud*S-lam[dd]*S
      dL <-lam[dd]*S-mud*L-alpha*L
      dB <-alpha*L-mud*B-beta*B
      dC <-beta*B-mud*C-deltaone*chi*C-deltatwo*vee*C-(1-chi)*(1-vee)*deltathree*C
      dD <-deltaone*chi*C-(mud+muc)*D-gamma*D
      dTr <-deltatwo*vee*C+gamma*D-mud*Tr-sigmaone*(1-rho)*Tr-sigmatwo*rho*Tr
      dR <-sigmaone*(1-rho)*Tr+eps*G+(1-chi)*(1-vee)*deltathree*C-mud*R-omega*(1-kappa)*R-lam[dd]*kappa*R
      dLr <-lam[dd]*kappa*R-mud*Lr-alpha*Lr
      dBr <-alpha*Lr-mud*Br-beta*Br
      dG <-sigmatwo*rho*Tr+beta*Br-mud*G-eps*G
      dCumInc<-beta*B+deltaone*chi*C
      # RATE OF CHANGE - MOSQUITO
      dSm <-mubm*Pm-mum*Sm-mlam*Sm
      dEm <-mlam*Sm-mum*Em-mgamma*Em
      dIm <-mgamma*Em-mum*Im
    }
    
    # RETURN THE RATE OF CHANGE
    output<-(c(dS, dL, dB, dC, dD, dTr, dR, dLr, dBr, dG, dCumInc, dSm, dEm, dIm))
    list(output)
  })
}

# MODEL INITIAL CONDITIONS - HUMAN
initP<-26299606  # population size
initL<-0         # Liver (exo-erythrocytic) stage compartment
initB<-0   # Asexual blood (erythrocytic) stage compartment
initC<-(98/100)*1059718    # Gametocyte (sexual) stage compartment with clinical disease
initD<-(2/100)*1059718          # Severe malaria compartment
initTr<-0      # Treatment compartment
initR<-0       # Recovered
initLr<-0         # Liver stage compartment (with clinical immunity)
initBr<-0         # Asexual blood stage compartment (with clinical immunity)
initG<-0         # Gametocyte stage compartment (with clinical immunity)

#INITIAL SUSCEPTIBLE COMPARTMENT - HUMAN
initS<-initP-initL-initB-initC-initD-initTr-initR-initLr-initBr-initG #Susceptible, fully naive


#MODEL INITIAL CONDITIONS - MOSQUITO
initPm<-26299606       # Mosquito population size
initEm<-1               # Initial exposed mosquito population
initIm<-597001          # Initial infected mosquito population

#INITIAL SUSCEPTIBLE COMPARTMENT - MOSQUITO
initSm<-initPm-initEm-initIm        #Initial susceptible malaria population

#DEFINE THE MODEL VARIABLES
start <- c(S = initS, L = initL, B = initB, C = initC, D = initD, Tr = initTr, R = initR, Lr = initLr, Br = initBr, G = initG, CumInc=1059718
           , Sm = initSm, Em = initEm, Im = initIm)



## The parameters 
parms <- c(logdeltaone=log(1/(5/365)), logdeltatwo=log(1/(3/365)), logdeltathree=log(1/(2.5/52)), logvee=log(8.9/100), logchi=log(2/100),  loggamma=log(1/(2/365)), logsigmatwo=log(1/(42/365)), logkappa=log(90/100), logitamp=logit(0.1), logitphi=logit(0.99))

period<-7 #years

#vector of time steps
times <- seq(0, period, by = 1) # in years, by 1 year intervals.

#RUN THE MODEL
run_out<-ode(times=times, y=start, func=malaria,parms=parms)

#Calculating incidence 
#Comparing model plot Inc and data Cases
inc<-c(run_out[1,12])
for (tt in 2:length(times)){
  inc<-c(inc,(run_out[tt,12] - run_out[tt-1, 12]))
}

#factor in changing report rate over the years
#report<-19%
report<-c(1,0.52,0.62,0.3,0.1,0.23,0.7,0.8)
par(mfrow=c(1,1))
plot(times, data$Cases, type="o",lwd=3, col="red", ylim=c(0,2500000))
lines(times,report*inc, col=1, type="o",lwd=2,lty=2)

malaria.sse<-function(data,parms){
  model<-ode(times=times, y=start, func=malaria,parms=parms)
  inc<-c(run_out[1,12])
  for (tt in 2:length(times)){
    inc<-c(inc,(run_out[tt,12] - run_out[tt-1, 12]))
  }
  for(ui in 1:7){
    error<-report[ui]*inc[ui]-data$Cases
    sse<-sum(error^2)
  }
  return(sse)
}

malaria.sse(data,parms)


fit<-optim(parms,malaria.sse,data=data)

parmtr<-function(parms){
  deltaone=exp(parms[1])
  deltatwo=exp(parms[2])
  deltathree=exp(parms[3])
  vee=exp(parms[4])
  chi=exp(parms[5])
  gamma=exp(parms[6])
  sigmatwo=exp(parms[7])
  kappa=exp(parms[8])
  amp=inv.logit(parms[9])
  phi=inv.logit(parms[10])
  output<-c(deltaone=deltaone, deltatwo=deltatwo, deltathree=deltathree, vee=vee, chi=chi, gamma=gamma, sigmatwo=sigmatwo, kappa=kappa, amp=amp, phi=phi) 
  return(output)
}

bestpar<-parmtr(fit$par)
startpar<-parmtr(parms)
cbind(startpar,bestpar)

run_out<-ode(times=times, y=start, func=malaria,parms=parms)

run_fit<-ode(times=times, y=start, func=malaria,parms=fit$par)
incfit<-c(run_fit[1,12])
for (tt in 2:length(times)){
  incfit<-c(incfit,run_fit[tt,12] - run_fit[tt-1,12])
}

# Plot the output
par(mfrow=c(1,1))
plot(2004+times,data$Cases, type="b",lwd=2,main = "Malaria Transmission Fitting Plot",xlab = "Time in years",ylab="Incident cases",ylim=c(0,3000000))
lines(2004+times,report*inc,col="red",lwd=2,lty=1,main = "Malaria Transmission Fitting Plot",xlab = "Time in years",ylab="Incident cases",ylim=c(0,3000000))
legend("topleft", c("data", "model"), fill = c("black","red"))


