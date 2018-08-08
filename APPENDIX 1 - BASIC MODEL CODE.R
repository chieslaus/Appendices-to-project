#APPENDIX 1 - BASIC MODEL CODE
##############################

#MALARIA MATHEMATICAL MODEL AND COST-EFFECTIVENESS ANALYSIS 2018: R CODE 2018
version #R version 3.4.3 (2017-11-30) nickname       Kite-Eating Tree

#CLEAR WORKSPACE
rm(list=ls())
setwd("C:/Users/Lausdeus Chiegboka/Documents/ox docs/PLACEMENT PROJECT/MAIN/data")
getwd()

data<-read.csv('South-West Malaria data from 2004 to 2011.csv')
plot(data$Cases)
plot(data$CumCases)

#LOAD REQUIRED LIBRARIES
library(deSolve)

year_start <- 13 # Start time for model simulation
year_stop<- 23 #years

times <- seq(year_start, year_stop, by = 1) # in years, by 1 year intervals.

## The model parameters 
parameters <- c(mud=(1/54),                 # Death rate (1/life expectancy)
                muc=(10.3/100),             # Severe Malaria case fatality rate
                mub=((1/54)+(2.6/100)),     # Birth rate (1/life expectancy+population growth rate)
                omega=(1/5),                # Rate of loss of immunity = 1/(average duration of immunity)
                alpha=(1/(7.5/365)),        # Rate of movement from liver stage to asexual blood stage = 1/(average prepatent period)
                beta=(1/(2.5/365)),         # Rate of movement  from asexual blood stage to gametocyte stage = 1/(latent period - prepatent period)
                deltaone=(1/(5/365)),       # Rate of development of severe disease = 1/duration of clinical disease before onset of severe malaria
                deltatwo=(1/(3/365)),       # Rate of treatment = 1/time lapse in infectious stage before seeking treatment
                deltathree=(1/(2.5/52)),    # Rate of symptom resolution without treatment = 1/duration of clinical disease without treatment
                vee=(8.9/100),              # Treatment coverage i.e., Proportion of clinical malaria that receive treatment = treatment coverage of ACTs
                chi=(2/100),                # Proportion of clinical malaria cases that become severe
                gamma=(1/(2/365)),          # Rate of treatment in severe malaria = 1/time lapse before seeking treatment of severe disease
                sigmaone=(1/(32/24/365)),   # Rate of recovery following treatment = 1/parasite clearance time
                sigmatwo=(1/(42/365)),      # Rate of recrudescence = 1/time to recrudescence
                rho=(4.2/100),              # Proportion of treated cases with treatment failure (recrudescence)
                kappa=(90/100),             # Proportion of clinically immune that develop latent infection
                eps=(1/(18/12)),            # Rate of natural recovery from infectiousness (among the clinically immune)
                bee=0.5*365,                # Mosquito biting rate per man/year
                amp=0.1,                    # Relative amplitude of seasonal forcing
                phi=0.99,                   # Year of peak in seasonal forcing
                
                llin_eff = 0.5,              # Efficacy of llins in reducing disease incidence
                irs_eff = 0,                # Efficacy of IRS in reducing disease incidence
                irs_cov_i = 0,              # Intervention coverage of IRS
                year_interv = 0,            # Years after first case when the intervention starts
                
                #mosquito parameters
                mubm=(1/(12.24/365)),       # Mosquito birth rate = 1/life expectancy of mosquitoes
                mum=(1/(12.24/365)),        # Mosquito death rate = 1/life expectancy of mosquitoes
                mgamma=(1/(11/365)),        # Rate of movement from exposed to infectiousness in mosquitoes = 1/latent period
                ai=(5/100),                 # Per bite probability of transmission from infected mosquito to human
                aim=(47/100)                # Per bite probability of transmission from infectious human to mosquito 
                
)


# MODEL INITIAL CONDITIONS - HUMAN
initP<-26299606                             # Population size
initL<-0                                    # Liver (exo-erythrocytic) stage compartment
initB<-0                                    # Asexual blood (erythrocytic) stage compartment
initC<-(98/100)*1059718                     # Gametocyte (sexual) stage compartment with clinical disease
initD<-(2/100)*1059718                      # Severe malaria compartment
initTr<-0                                   # Treatment compartment
initR<-0                                    # Recovered compartment
initLr<-0                                   # Liver stage compartment (with clinical immunity)
initBr<-0                                   # Asexual blood stage compartment (with clinical immunity)
initG<-0                                    # Gametocyte stage compartment (with clinical immunity)

#INITIAL SUSCEPTIBLE COMPARTMENT - HUMAN
initS<-initP-initL-initB-initC-initD-initTr-initR-initLr-initBr-initG #Susceptible, fully naive


#MODEL INITIAL CONDITIONS - MOSQUITO
initPm<-26299606                            # Mosquito population size
initEm<-1                                   # Initial exposed mosquito population
initIm<-597001                              # Initial infected mosquito population

#INITIAL SUSCEPTIBLE COMPARTMENT - MOSQUITO
initSm<-initPm-initEm-initIm                # Initial susceptible malaria population

#DEFINE THE MODEL VARIABLES
start <- c(S = initS, L = initL, B = initB, C = initC, D = initD, Tr = initTr, R = initR, Lr = initLr, Br = initBr, G = initG, CumInc=1059718
           , Sm = initSm, Em = initEm, Im = initIm)


malaria <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    # DEFINE VARIABLES
    P <- (S+L+B+C+D+Tr+R+Lr+Br+G)
    Pm <- (Sm+Em+Im)
    
    seas<-1+amp*cos(2*pi*(t-phi))
    #llin_cov <- (t>=(year_start+year_interv))*llin_cov_i
    llin_cov=c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.2,0.2,0.2,0.387,0.387,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53)
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


#RUN THE MODEL
run_out<-ode(times=times, y=start, func=malaria,parms=parameters)

#Reporting parameter = Proportion of cases reported
report<-0.8

#MODEL OUTPUT
#Dataframe for run_out
RUN_OUT<-data.frame(run_out)

#calculating the human/host population size at each time step
POPh<-RUN_OUT$S+RUN_OUT$L+RUN_OUT$B+RUN_OUT$C+RUN_OUT$D+RUN_OUT$Tr+RUN_OUT$R+RUN_OUT$Lr+RUN_OUT$Br+RUN_OUT$G

#calculating the vector size at each time step
POPm<-RUN_OUT$Sm+RUN_OUT$Em+RUN_OUT$Im

#Cumulative incidence
Cum<-RUN_OUT$CumInc #cumulative report
Cum_inc<-round((Cum/POPh)*1000) # cumulative incidence per 1000 of population
Cum_inc

#Calculating incidence 
inc<-c(run_out[1,12])
for (tt in 2:length(times)){
  inc<-c(inc,(run_out[tt,12] - run_out[tt-1, 12]))
}

#Calculating Deaths from Malaria
Deaths<-(parameters["muc"]*RUN_OUT$D)
TotalDeaths<-round(sum(parameters["muc"]*RUN_OUT$D))


#Prevalence of gametocytes
Gam<-RUN_OUT$C+RUN_OUT$D+RUN_OUT$G

#Percentage of gametocytes in the population
PercGam<-Gam/POPh*100

#Total parasite prevalence (all blood stages, asexual and sexual)
ParPrev<-RUN_OUT$B+RUN_OUT$C+RUN_OUT$D+RUN_OUT$G+RUN_OUT$Br

#Percentage parasite prevalence
ParPrevPerc<-(ParPrev/POPh)*100

#Percentage of infectious cases with immunity
PercInfImm<- RUN_OUT$G/Gam*100

#MODEL PLOTS

#options(scipen=1000)

# a simple plot of the model output
plot(run_out)

par(mfrow=c(1,1))
#Plot of cumulative incidence over the time period
plot(2004+times, report*Cum_inc, col="red", type="b", lwd="2", ylim = c(0,1000), xlab="Time in Years", ylab="Cumulative Incidence (in thousands)", main="Cumulative Incidence")

par(mfrow=c(1,1))
#Plot of incidence over time
plot(2004+times, report*inc, col="red", type="l", lwd="1.5", xlab = "Time in Years", ylab = "Incident Cases", main="Reported Incident Cases over Time")

#Plot of Parasite Prevalence Percentage over time
plot(2004+times, ParPrevPerc, ylim=c(0,90), col="red", type="b", lwd="2", xlab = "Time in Years", ylab = "Parasite prevalence (%)", main="Parasite Prevalence over Time", cex.main=0.8)
