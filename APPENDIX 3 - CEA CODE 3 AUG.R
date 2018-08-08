#APPENDIX 3 - CEA MODEL CODE
############################

#MALARIA MATHEMATICAL MODEL AND COST-EFFECTIVENESS ANALYSIS 2018: R CODE 2018
version #R version 3.4.3 (2017-11-30) nickname       Kite-Eating Tree

#CLEAR WORKSPACE
rm(list=ls())
setwd("C:/Users/Lausdeus Chiegboka/Documents/ox docs/PLACEMENT PROJECT/MAIN/data")
getwd()

data<-read.csv('South-West Malaria data from 2004 to 2011.csv')
plot(data$CumCases)

#LOAD REQUIRED LIBRARIES
library(deSolve)

year_start <- 13  # Start time for model simulation
year_stop<-23     # years

times <- seq(year_start, year_stop, by = 1) # in years, by 1 year intervals.

# MODEL INITIAL CONDITIONS - HUMAN
initP<-26299606         # Population size
initL<-0                # Liver (exo-erythrocytic) stage compartment
initB<-0                # Asexual blood (erythrocytic) stage compartment
initC<-(98/100)*1059718 # Gametocyte (sexual) stage compartment with clinical disease
initD<-(2/100)*1059718  # Severe malaria compartment
initTr<-0               # Treatment compartment
initR<-0                # Recovered compartment
initLr<-0               # Liver stage compartment (with clinical immunity)
initBr<-0               # Asexual blood stage compartment (with clinical immunity)
initG<-0                # Gametocyte stage compartment (with clinical immunity)

#INITIAL SUSCEPTIBLE COMPARTMENT - HUMAN
initS<-initP-initL-initB-initC-initD-initTr-initR-initLr-initBr-initG #Susceptible, fully naive


#MODEL INITIAL CONDITIONS - MOSQUITO
initPm<-26299606        # Mosquito population size
initEm<-1               # Initial exposed mosquito population
initIm<-597001          # Initial infected mosquito population

#INITIAL SUSCEPTIBLE COMPARTMENT - MOSQUITO
initSm<-initPm-initEm-initIm        #Initial susceptible malaria population

#DEFINE THE MODEL VARIABLES
start <- c(S = initS, L = initL, B = initB, C = initC, D = initD, Tr = initTr, R = initR, Lr = initLr, Br = initBr, G = initG, CumInc=1059718
           , Sm = initSm, Em = initEm, Im = initIm)

#Base case (low LLIN intervention coverage, no IRS)
#################################################
malaria <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    # DEFINE VARIABLES
    P <- (S+L+B+C+D+Tr+R+Lr+Br+G)
    Pm <- (Sm+Em+Im)
    
    seas<-1+amp*cos(2*pi*(t-phi))
    llin_cov=c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.2,0.2,0.2,0.387,0.387,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53,0.53)
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


## The base case model parameters 
base_par <- c(mud=(1/54),                 # Death rate (1/life expectancy)
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
              
              llin_eff = 0.5,              # Efficacy of LLINs in reducing disease incidence
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

#INTERVENTION 1: SCALE-UP OF LLIN TO 95%
########################################
llin_malaria <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    # DEFINE VARIABLES
    P <- (S+L+B+C+D+Tr+R+Lr+Br+G)
    Pm <- (Sm+Em+Im)
    
    seas<-1+amp*cos(2*pi*(t-phi))
    llin_cov=c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.2,0.2,0.2,0.387,0.387,0.53,0.53,0.53,0.572,0.614,0.656,0.698,0.74,0.782,0.824,0.866,0.908,0.95)
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


## LLIN SCALE-UP TO 95%
#######################
llin_par <- c(mud=(1/54),                 # Death rate (1/life expectancy)
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
              llin_cov_i = 0.75,          # Intervention coverage of LLIN
              llin_eff = 0.5,             # Efficacy of LLINs in reducing disease incidence
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

llin_scenario_malaria <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    # DEFINE VARIABLES
    P <- (S+L+B+C+D+Tr+R+Lr+Br+G)
    Pm <- (Sm+Em+Im)
    
    seas<-1+amp*cos(2*pi*(t-phi))
    llin_cov=c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.2,0.2,0.2,0.387,0.387,0.53,0.53,0.53,0.552,0.574,0.596,0.618,0.640,0.662,0.684,0.706,0.728,0.75)
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

## INTERVENTION 2: SCALE-UP OF LLIN TO 75%
##########################################
llin_scenario <- c(mud=(1/54),                 # Death rate (1/life expectancy)
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
                   llin_cov_i = 0.75,          # Intervention coverage of LLIN
                   llin_eff = 0.5,             # Efficacy of LLINs in reducing disease incidence
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

#INTERVENTION 3: SCALE-UP TO 95% LLIN+50%IRS
############################################
irs_malaria <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    # DEFINE VARIABLES
    P <- (S+L+B+C+D+Tr+R+Lr+Br+G)
    Pm <- (Sm+Em+Im)
    
    seas<-1+amp*cos(2*pi*(t-phi))
    irs_cov <- (t>=(year_start+year_interv))*irs_cov_i
    irs <- (1-irs_eff*irs_cov)
    llin_cov=c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.2,0.2,0.2,0.387,0.387,0.53,0.53,0.53,0.572,0.614,0.656,0.698,0.74,0.782,0.824,0.866,0.908,0.95)
    llin<-NULL
    for(jj in 1:length(llin_cov)){
      llin <- c(llin,1-llin_eff*llin_cov[jj])
    }
    lam<-NULL
    for(hh in 1:length(llin)){
      lam<-c(lam,llin[hh]*irs*seas*bee*ai*Im/Pm)
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

#LLIN+IRS model parameters
irs_par <- c(mud=(1/54),                 # Death rate (1/life expectancy)
             muc=(10.3/100),             # Severe Malaria case fatality rate
             mub=((1/54)+(2.6/100)),     # Birth rate (1/life expectancy + population growth rate)
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
             
             llin_cov_i = 0.95,           # Intervention coverage of LLIN
             llin_eff = 0.5,              # Efficacy of LLINs in reducing disease incidence
             irs_eff = 0.14,             # Efficacy of IRS in reducing disease incidence
             irs_cov_i = 0.5,            # Intervention coverage of IRS
             year_interv = 0,            # Years after first case when the intervention starts
             
             #mosquito parameters
             mubm=(1/(12.24/365)),       # Mosquito birth rate = 1/life expectancy of mosquitoes
             mum=(1/(12.24/365)),        # Mosquito death rate = 1/life expectancy of mosquitoes
             mgamma=(1/(11/365)),        # Rate of movement from exposed to infectiousness in mosquitoes = 1/latent period
             ai=(5/100),                 # Per bite probability of transmission from infected mosquito to human
             aim=(47/100)                # Per bite probability of transmission from infectious human to mosquito 
             
)

# INTERVENTION 4: SCALE-UP TO 75% LLIN+50% IRS
##############################################
#LLIN+IRS model parameters
new_irs_scenario_malaria <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    # DEFINE VARIABLES
    P <- (S+L+B+C+D+Tr+R+Lr+Br+G)
    Pm <- (Sm+Em+Im)
    
    seas<-1+amp*cos(2*pi*(t-phi))
    irs_cov <- (t>=(year_start+year_interv))*irs_cov_i
    irs <- (1-irs_eff*irs_cov)
    llin_cov=c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.2,0.2,0.2,0.387,0.387,0.53,0.53,0.53,0.552,0.574,0.596,0.618,0.640,0.662,0.684,0.706,0.728,0.75)
    llin<-NULL
    for(jj in 1:length(llin_cov)){
      llin <- c(llin,1-llin_eff*llin_cov[jj])
    }
    lam<-NULL
    for(hh in 1:length(llin)){
      lam<-c(lam,llin[hh]*irs*seas*bee*ai*Im/Pm)
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

new_irs_scenario <- c(mud=(1/54),                 # Death rate (1/life expectancy)
                      muc=(10.3/100),             # Severe Malaria case fatality rate
                      mub=((1/54)+(2.6/100)),     # Birth rate (1/life expectancy + population growth rate)
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
                      
                      llin_cov_i = 0.75,           # Intervention coverage of LLIN
                      llin_eff = 0.5,              # Efficacy of LLINs in reducing disease incidence
                      irs_eff = 0.14,             # Efficacy of IRS in reducing disease incidence
                      irs_cov_i = 0.5,            # Intervention coverage of IRS
                      year_interv = 0,            # Years after first case when the intervention starts
                      
                      #mosquito parameters
                      mubm=(1/(12.24/365)),       # Mosquito birth rate = 1/life expectancy of mosquitoes
                      mum=(1/(12.24/365)),        # Mosquito death rate = 1/life expectancy of mosquitoes
                      mgamma=(1/(11/365)),        # Rate of movement from exposed to infectiousness in mosquitoes = 1/latent period
                      ai=(5/100),                 # Per bite probability of transmission from infected mosquito to human
                      aim=(47/100)                # Per bite probability of transmission from infectious human to mosquito 
                      
)

#RUN THE MODEL
##############
out_base<-ode(y=start,times=times, func=malaria,parms=base_par)
out_llin <- ode(y=start, times = times, func = llin_malaria, parms = llin_par)
out_irs <- ode(y=start, times = times, func = irs_malaria, parms = irs_par)
out_llin_scenario <- ode(y=start, times = times, func = llin_scenario_malaria, parms = llin_scenario)
out_new_irs_scenario <- ode(y=start, times = times, func = new_irs_scenario_malaria, parms = new_irs_scenario)

#Population size
poph<-out_base[,"S"]+out_base[,"L"]+out_base[,"B"]+out_base[,"C"]+out_base[,"D"]+out_base[,"Tr"]+out_base[,"R"]+out_base[,"Lr"]+out_base[,"Br"]+out_base[,"G"]
mean_poph<- mean(poph)
mean_poph

report<-0.8 # 80% report rate

#Outputs of all interventions
##############################
#Total cases from model prediction
cumreports_base<- report*out_base[,12]
cumreports_llin<- report*out_llin[,12]
cumreports_llin_scenario<- report*out_llin_scenario[,12]
cumreports_irs<- report*out_irs[,12]
cumreports_new_irs_scenario<- report*out_new_irs_scenario[,12]

Total_case_base <- round(cumreports_base[11]) #the end year of simulation
Total_case_llin <- round(cumreports_llin[11])
Total_case_llin_scenario <- round(cumreports_llin_scenario[11])
Total_case_irs <- round(cumreports_irs[11])
Total_case_new_irs_scenario <- round(cumreports_new_irs_scenario[11])

Total_case_base
Total_case_llin
Total_case_llin_scenario
Total_case_irs
Total_case_new_irs_scenario

#Discounted total cases
# 3% Discount rate
##################
disc_rate<- 0.03 # 3% Discount rate
n<-10 # Number of years
Disc_Total_case_base<-Total_case_base/((1+disc_rate)^n);Disc_Total_case_base
Disc_Total_case_llin<-Total_case_llin/((1+disc_rate)^n); Disc_Total_case_llin
Disc_Total_case_llin_scenario<-Total_case_llin_scenario/((1+disc_rate)^n);Disc_Total_case_llin_scenario
Disc_Total_case_irs<-Total_case_irs/((1+disc_rate)^n);Disc_Total_case_irs
Disc_Total_case_new_irs_scenario<-Total_case_new_irs_scenario/((1+disc_rate)^n);Disc_Total_case_new_irs_scenario

#plot of total cases
Total_cases<- c(Total_case_base/1000000, Total_case_llin/1000000, Total_case_llin_scenario/1000000, 
                Total_case_irs/1000000, Total_case_new_irs_scenario/1000000)
barplot(Total_cases, col=c("red", "brown", "blue", "black", "darkgrey"), xaxt="n", ylim = c(0,35), ylab = "Total cases (Millions)", main = "Predicted Total Number of Cases under Different Interventions", cex.main=0.8)
legend("topright", legend=c("Baseline", "LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", "LLIN(75%)+IRS(50%)"), 
       fill = c("red", "brown", "blue", "black", "darkgrey"), cex=0.7)

#Cumulative incidence per 1000 of population
############################################
cum_inc_base<-round((cumreports_base/poph)*1000)
cum_inc_llin<-round((cumreports_llin/poph)*1000)
cum_inc_llin_scenario<-round((cumreports_llin_scenario/poph)*1000)
cum_inc_irs<-round((cumreports_irs/poph)*1000)
cum_inc_new_irs_scenario<-round((cumreports_new_irs_scenario/poph)*1000)

cum_inc_base
cum_inc_llin
cum_inc_llin_scenario
cum_inc_irs
cum_inc_new_irs_scenario

#FInal cumulative incidence value in base case
endcuminc<-cum_inc_base[11]
endcuminc

# Plot cumulative incidence
###########################
plot(2004+times, cum_inc_base, type="b", lwd=2.0, col="red", xlab="Time in years", ylab="Cumulative incidence (per thousand)", main="Predicted Cumulative Incidence of Malaria Under Different Intervention States (per 1000)", cex.main=0.7)
lines(2004+times, cum_inc_llin, lwd=1.5, col="brown")
lines(2004+times, cum_inc_llin_scenario, lwd=1.5, col="blue", lty=3)
lines(2004+times, cum_inc_irs, lwd=1.5, col="black", lty=3)
lines(2004+times, cum_inc_new_irs_scenario, lwd=1.5, col="darkgrey", lty=3)
legend("topleft", c("Base Case", "LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", "LLIN(75%)+IRS(50%)"), col = c("red","brown","blue","black","darkgrey"), lty=c(1,1,3,3,3), lwd = c(2,1.5,1.5,1.5,1.5), cex=0.7)

# Zooming into the cumulative incidence
#######################################
plot(2004+times, cum_inc_base, type="b", lwd=2.0, col="red", ylim=c(0,100), xlab="Time in years", ylab="Cumulative incidence (per thousand)", main="Predicted Cumulative Incidence of Malaria Under Different Intervention States (per 1000)", cex.main=0.7)
lines(2004+times, cum_inc_llin, lwd=1.5, col="brown")
lines(2004+times, cum_inc_llin_scenario, lwd=1.5, col="blue", lty=3)
lines(2004+times, cum_inc_irs, lwd=1.5, col="black", lty=3)
lines(2004+times, cum_inc_new_irs_scenario, lwd=1.5, col="darkgrey", lty=3)
legend("topleft", c("Base Case", "LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", "LLIN(75%)+IRS(50%)"), col = c("red","brown","blue","black","darkgrey"), lty=c(1,1,3,3,3), lwd = c(2,1.5,1.5,1.5,1.5), cex=0.7)

#Calculating yearly incidence
############################# 
inc_base<-c(out_base[1,12])
for (tt in 2:length(times)){
  inc_base<-c(inc_base,(out_base[tt,12] - out_base[tt-1, 12]))
}

inc_llin<-c(out_llin[1,12])
for (tt in 2:length(times)){
  inc_llin<-c(inc_llin,(out_llin[tt,12] - out_llin[tt-1, 12]))
}
inc_llin_scenario<-c(out_llin_scenario[1,12])
for (tt in 2:length(times)){
  inc_llin_scenario<-c(inc_llin_scenario,(out_llin_scenario[tt,12] - out_llin_scenario[tt-1, 12]))
}

inc_irs<-c(out_irs[1,12])
for (tt in 2:length(times)){
  inc_irs<-c(inc_irs,(out_irs[tt,12] - out_irs[tt-1, 12]))
}

inc_new_irs_scenario<-c(out_new_irs_scenario[1,12])
for (tt in 2:length(times)){
  inc_new_irs_scenario<-c(inc_new_irs_scenario,(out_new_irs_scenario[tt,12] - out_new_irs_scenario[tt-1, 12]))
}
inc_base
inc_llin
inc_llin_scenario
inc_irs
inc_new_irs_scenario

#Incidence rate
rate_inc_base<-(inc_base/poph)*1000
rate_inc_llin<-(inc_llin/poph)*1000
rate_inc_llin_scenario<-(inc_llin_scenario/poph)*1000
rate_inc_irs<-(inc_irs/poph)*1000
rate_inc_new_irs_scenario<-(inc_new_irs_scenario/poph)*1000

rate_inc_base
rate_inc_llin
rate_inc_llin_scenario
rate_inc_irs
rate_inc_new_irs_scenario

#Plot of incidence
plot(2004+times, rate_inc_base, type="l", lwd=2.0, col="red", xlab="Time in years", ylab="Incidence rate (thousands)", main="Predicted Yearly Incidence of Malaria Under Different Intervention States", cex.main=1.0)
lines(2004+times, rate_inc_llin, lwd=1.8, col="brown")
lines(2004+times, rate_inc_llin_scenario, lwd=2.0, col="blue", lty=3)
lines(2004+times, rate_inc_irs, lwd=1.5, col="black", lty=3)
lines(2004+times, rate_inc_new_irs_scenario, lwd=2.0, col="darkgrey", lty=3)
legend("topleft", c("Baseline", "LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", "LLIN(75%)+IRS(50%)"), col = c("red","brown","blue","black","darkgrey"), lty=c(1,1,3,3,3), lwd = c(2,2,2,2,2), cex=1.0)

#Plot of incidence
plot(2004+times, rate_inc_base, type="l", ylim=c(0,100), lwd=2.0, col="red", xlab="Time in years", ylab="Incidence rate", main="Predicted Yearly Incidence of Malaria Under Different Intervention States (per 1000)", cex.main=0.7)
lines(2004+times, rate_inc_llin, lwd=1.5, col="brown")
lines(2004+times, rate_inc_llin_scenario, lwd=1.5, col="blue", lty=3)
lines(2004+times, rate_inc_irs, lwd=1.5, col="black", lty=3)
lines(2004+times, rate_inc_new_irs_scenario, lwd=1.5, col="darkgrey", lty=3)
legend("topleft", c("Base Case", "LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", "LLIN(75%)+IRS(50%)"), col = c("red","brown","blue","black","darkgrey"), lty=c(1,1,3,3,3), lwd = c(2,1.5,1.5,1.5,1.5), cex=0.8)

#Parasite Prevalence
####################
#Total parasite prevalence (all blood stages, asexual and sexual)
ParPrev_base<-out_base[,"B"]+out_base[,"C"]+out_base[,"D"]+out_base[,"G"]+out_base[,"Br"]
ParPrev_llin<-out_llin[,"B"]+out_llin[,"C"]+out_llin[,"D"]+out_llin[,"G"]+out_llin[,"Br"]
ParPrev_llin_scenario<-out_llin_scenario[,"B"]+out_llin_scenario[,"C"]+out_llin_scenario[,"D"]+out_llin_scenario[,"G"]+out_llin_scenario[,"Br"]
ParPrev_irs<-out_irs[,"B"]+out_irs[,"C"]+out_irs[,"D"]+out_irs[,"G"]+out_irs[,"Br"]
ParPrev_new_irs_scenario<-out_new_irs_scenario[,"B"]+out_new_irs_scenario[,"C"]+out_new_irs_scenario[,"D"]+out_new_irs_scenario[,"G"]+out_new_irs_scenario[,"Br"]

#Percentage parasite prevalence
ParPrevPerc_base<-(ParPrev_base/poph)*100
ParPrevPerc_llin<-(ParPrev_llin/poph)*100
ParPrevPerc_llin_scenario<-(ParPrev_llin_scenario/poph)*100
ParPrevPerc_irs<-(ParPrev_irs/poph)*100
ParPrevPerc_new_irs_scenario<-(ParPrev_new_irs_scenario/poph)*100

ParPrevPerc_base
ParPrevPerc_llin
ParPrevPerc_llin_scenario
ParPrevPerc_irs
ParPrevPerc_new_irs_scenario

#Plot percentage parasite prevalence
plot(2004+times, ParPrevPerc_base, type="b", lwd=2.0, col="red", ylim=c(0,100), xlab="Time in years", ylab="Parasite prevalence (%)", main="Predicted Malaria Parasite Prevalence under Different Intervention States", cex.main=0.7)
lines(2004+times, ParPrevPerc_llin, lwd=1.5, col="brown")
lines(2004+times, ParPrevPerc_llin_scenario, lwd=1.5, col="blue")
lines(2004+times, ParPrevPerc_irs, lwd=1.5, col="black", lty=3)
lines(2004+times, ParPrevPerc_new_irs_scenario, lwd=1.5, col="darkgrey", lty=3)
legend("topleft", c("Base Case", "LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", 
                    "LLIN(75%)+IRS(50%)"), col = c("red","brown","blue","black","darkgrey"), lty=c(1,1,3,3,3), lwd = c(2,1.5,1.5,1.5,1.5), cex=0.7)

#Percent reduction in Parasite prevalence
red_llin<-((ParPrevPerc_base-ParPrevPerc_llin)/ParPrevPerc_base)*100
red_llin_scenario<-((ParPrevPerc_base-ParPrevPerc_llin_scenario)/ParPrevPerc_base)*100
red_irs<-((ParPrevPerc_base-ParPrevPerc_irs)/ParPrevPerc_base)*100
red_new_irs_scenario<-((ParPrevPerc_base-ParPrevPerc_new_irs_scenario)/ParPrevPerc_base)*100

plot(2004+times, red_llin, type="l", lwd=2.0, col="brown", ylim=c(0,100), xlab="Time in years", 
     ylab="Reduction of parasite prevalence (%)", main="Predicted Malaria Parasite Prevalence Reduction Compared to Baseline", cex.main=0.8)
lines(2004+times, red_llin_scenario, lwd=2.0, col="blue")
lines(2004+times, red_irs, lwd=2.0, col="black", lty=3)
lines(2004+times, red_new_irs_scenario, lwd=2.0, col="darkgrey", lty=3)
legend("bottomright", c("LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", 
                        "LLIN(75%)+IRS(50%)"), col = c("brown","blue","black","darkgrey"), lty=c(1,1,3,3), lwd = c(2,2,2,2), cex=0.8)



#Mean Percent reduction in Parasite prevalence
mean_red_llin<-mean((ParPrevPerc_base-ParPrevPerc_llin)/ParPrevPerc_base)*100
mean_red_llin_scenario<-mean((ParPrevPerc_base-ParPrevPerc_llin_scenario)/ParPrevPerc_base)*100
mean_red_irs<-mean((ParPrevPerc_base-ParPrevPerc_irs)/ParPrevPerc_base)*100
mean_red_new_irs_scenario<-mean((ParPrevPerc_base-ParPrevPerc_new_irs_scenario)/ParPrevPerc_base)*100
Prevbarplot<-c(mean_red_llin, mean_red_llin_scenario, mean_red_irs, mean_red_new_irs_scenario)
barplot(Prevbarplot, col=c("brown", "blue", "black","darkgrey"), xaxt="n", ylab = "Reduction in parasite prevalence (%)", main = "Mean Cumulative Percent Reduction in Parasite Prevalence over ten years", 
        ylim=c(0, 100), cex.main=0.8)
legend(2,0, legend=c("LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", 
                     "LLIN(75%)+IRS(50%)"), col=c("brown", "blue", "black", "darkgrey"), 
       lty=1, lwd = 5, cex=0.7, bty = "n", xpd=TRUE, inset=c(-0.3,0))


#Death Outcomes (number of deaths)
case_fat<-0.076 # 7.6% case fatality rate

Total_Death_base<-round(case_fat*Total_case_base)
Total_Death_llin<-round(case_fat*Total_case_llin)
Total_Death_llin_scenario<-round(case_fat*Total_case_llin_scenario)
Total_Death_irs<-round(case_fat*Total_case_irs)
Total_Death_new_irs_scenario<-round(case_fat*Total_case_new_irs_scenario)

Total_Death_base
Total_Death_llin
Total_Death_llin_scenario
Total_Death_irs
Total_Death_new_irs_scenario


#Plot deaths
deathplot<- c(Total_Death_base, Total_Death_llin,Total_Death_llin_scenario, Total_Death_irs, Total_Death_new_irs_scenario)/1000
barplot(deathplot, col=c("red", "brown", "blue", "black", "darkgrey"), xaxt="n", ylab = "Predicted deaths (in thousands)", main = "Predicted Deaths Under Each Strategy", ylim=c(0, 2500), cex.main=0.8)
legend("topright", legend=c("Baseline", "LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", 
                            "LLIN(75%)+IRS(50%)"), fill = c("red", "brown", "blue", "black", "darkgrey"), cex=0.7)

#Calculating costs and DALYs
############################
#Costs

C_llin <- 2.61                     # Cost of LLIN per unit ($US)
Num_household <- mean_poph/4.6     # Number of household in the population (average 4.6 people per household) 
Num_household_eff<-Num_household/2 # Number of households for effective coverage
C_irs<- 2.38                       # Cost of IRS per household ($US)
C_diagnosis <- 4.03                # Cost of diagnosis of malaria
C_ACT <- 4.13                      # Unit cost of ACT
C_medical <- C_diagnosis+ C_ACT    # Average medical cost of malaria per case ($US)

C_medical_base<-C_medical*Total_case_base
c_int_base<-C_llin*Num_household_eff*0.53   ## Cost of interventions at Baseline LLIN coverage = 53%

C_medical_llin<-C_medical*Total_case_llin
c_int_llin<-C_llin*Num_household_eff*llin_par[['llin_cov_i']] # Cost of interventions

C_medical_llin_scenario<-C_medical*Total_case_llin_scenario
c_int_llin_scenario<-C_llin*Num_household_eff*llin_scenario[['llin_cov_i']]

C_medical_irs<-C_medical*Total_case_irs
c_int_irs<-C_irs*Num_household*irs_par[['irs_cov_i']]+C_llin*Num_household_eff*irs_par[['llin_cov_i']]

C_medical_new_irs_scenario<-C_medical*Total_case_new_irs_scenario
c_int_new_irs_scenario<-C_irs*Num_household*new_irs_scenario[['irs_cov_i']]+C_llin*Num_household_eff*new_irs_scenario[['llin_cov_i']]


Total_cost_base <- C_medical_base+c_int_base #Base coverage of 53% in 2017
Total_cost_llin <- c_int_llin+C_medical_llin
Total_cost_llin_scenario <- c_int_llin_scenario+C_medical_llin_scenario
Total_cost_irs <- c_int_irs+C_medical_irs
Total_cost_new_irs_scenario <- c_int_new_irs_scenario+C_medical_new_irs_scenario

#Component bar chart for the cost breakdown
costtable<- cbind(c((C_medical_base)/1000000,(c_int_base)/1000000), c((C_medical_llin)/1000000, (c_int_llin)/1000000), 
                  c((C_medical_llin_scenario)/1000000, (c_int_llin_scenario)/1000000), 
                  c((C_medical_irs)/1000000,(c_int_irs)/1000000),c((C_medical_new_irs_scenario)/1000000,(c_int_new_irs_scenario)/1000000)); costtable
colnames(costtable)<-c("Base", "LLIN(95%)", "LLIN(75%)", "LLIN(95%)+IRS", 
                       "LLIN(75%)+IRS")

rownames(costtable)<-c("Medical cost", "Intervention cost")
costtable

#Barplot of overall costs
barplot(costtable,
        col = c("darkgrey","white"), legend =rownames(costtable), ylim = c(0,300), ylab = "Total costs (US$ Million)", 
        main = "Total Cost Breakdown", cex.main=1.2)

#Plot of total costs
costplot<-c((Total_cost_base)/1000000, (Total_cost_llin)/1000000, 
            (Total_cost_llin_scenario)/1000000, (Total_cost_irs)/1000000, (Total_cost_new_irs_scenario)/1000000)
barplot(costplot, col=c("red", "brown", "blue", "black", "darkgrey"), xaxt="n", ylab = "Total costs of interventions (Million US$)", main = "Total Cost of Each Strategy", ylim= c(0,300), cex.main=0.8)
legend("topright", legend=c("Base Case", "LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", 
                            "LLIN(75%)+IRS(50%)"), col=c("red", "brown", "blue", "black", "darkgrey"), lty=1, lwd = 5, cex=0.7)

#Health outcomes
DW<- 0.2             # Average disability weight for malaria
dur_dz<- 2.5/52      # Duration of clinical disease = 2.5 weeks
YLD.case<- DW*dur_dz # Years lived with disability per case 
YLL.death<- 50       # Life years lost per death

# Calculate DALY
################
DALYsloss_base <-Total_Death_base*YLL.death+Total_case_base*YLD.case
DALYsloss_llin <- Total_Death_llin*YLL.death+Total_case_llin*YLD.case
DALYsloss_llin_scenario <- Total_Death_llin_scenario*YLL.death+Total_case_llin_scenario*YLD.case
DALYsloss_irs <- Total_Death_irs*YLL.death+Total_case_irs*YLD.case
DALYsloss_new_irs_scenario <- Total_Death_new_irs_scenario*YLL.death+Total_case_new_irs_scenario*YLD.case

DALYsloss_base
DALYsloss_llin
DALYsloss_llin_scenario
DALYsloss_irs
DALYsloss_new_irs_scenario

#Cost-Effectiveness analysis in terms of DALYs
##############################################
# Scale-up of interventions vs current practice
#Calculations of incremental costs and incremental outomes done on Excel Spreadsheets attached as an annex to this.
Inc_cost_llin<--13292280
Inc_cost_llin_scenario<- -214256222
Inc_cost_irs<- 7766742
Inc_cost_new_irs_scenario<- -9075292


DALY_averted_llin<- 2576152
DALY_averted_llin_scenario<- 100904330
DALY_averted_irs<- 799618
DALY_averted_new_irs_scenario<- 7866455

ICER_llin<- Inc_cost_llin/DALY_averted_llin
ICER_llin_scenario<- Inc_cost_llin_scenario/DALY_averted_llin_scenario
ICER_irs<- Inc_cost_irs/DALY_averted_irs
ICER_new_irs_scenario<- Inc_cost_new_irs_scenario/DALY_averted_new_irs_scenario

ICER_llin
ICER_llin_scenario
ICER_irs
ICER_new_irs_scenario

# Create a matrix for the ICER
Malaria_res<-matrix(NA, ncol=3, nrow=4)
Malaria_res[1,]<-c(Inc_cost_llin, DALY_averted_llin, ICER_llin)
Malaria_res[2,]<-c(Inc_cost_llin_scenario, DALY_averted_llin_scenario, ICER_llin_scenario)
Malaria_res[3,]<-c(Inc_cost_irs, DALY_averted_irs, ICER_irs)
Malaria_res[4,]<-c(Inc_cost_new_irs_scenario, DALY_averted_new_irs_scenario, ICER_new_irs_scenario)
colnames(Malaria_res)<-c("Incremental Cost", "DALYs_averted", "ICER")
rownames(Malaria_res)<-c("LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", 
                         "LLIN(75%)+IRS(50%)")
Malaria_res

# Plot to see results on the ICER plane
plot((Malaria_res[,2])/1000000, (Malaria_res[,1])/1000000, xlim=c(-20,100), ylim=c(-220,100), ylab="Incremental costs (Million US$)", 
     xlab="DALY Averted (in millions)", pch=18, col=c("red", "brown", "blue", "black", "darkgrey"), main="CEA of Interventions for Malaria", cex.main=0.8)
#text((Malaria_res[,2])/1000000, (Malaria_res[,1])/1000000,label=c("LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", "LLIN(75%)+IRS(50%)"), cex=0.7, pos=2)
legend("topright", legend=c("LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", 
                         "LLIN(75%)+IRS(50%)"), fill = c("brown", "blue", "black", "darkgrey"), cex=0.7)

abline(a=0, b=1969, lty=3, col='black', lwd=2) #b = GDP per capita of Nigeria in 2017
abline(v = 0, h = 0, lty=1, col="blue",lwd=1)

#Discounted ICER

Disc_Inc_cost_llin<--9890704
Disc_Inc_cost_llin_scenario<- -159426753
Disc_Inc_cost_irs<- 5779185
Disc_Inc_cost_new_irs_scenario<- -6752869


Disc_DALY_averted_llin<- 1974406
Disc_DALY_averted_llin_scenario<- 77334766
Disc_DALY_averted_irs<- 612841
Disc_DALY_averted_new_irs_scenario<- 6028983


Disc_ICER_llin<- Disc_Inc_cost_llin/Disc_DALY_averted_llin
Disc_ICER_llin_scenario<- Disc_Inc_cost_llin_scenario/Disc_DALY_averted_llin_scenario
Disc_ICER_irs<- Disc_Inc_cost_irs/Disc_DALY_averted_irs
Disc_ICER_new_irs_scenario<- Disc_Inc_cost_new_irs_scenario/Disc_DALY_averted_new_irs_scenario

Disc_ICER_llin
Disc_ICER_llin_scenario
Disc_ICER_irs
Disc_ICER_new_irs_scenario

# Create a matrix for the ICER
Disc_Malaria_res<-matrix(NA, ncol=3, nrow=4)
Disc_Malaria_res[1,]<-c(Disc_Inc_cost_llin, Disc_DALY_averted_llin, Disc_ICER_llin)
Disc_Malaria_res[2,]<-c(Disc_Inc_cost_llin_scenario, Disc_DALY_averted_llin_scenario, Disc_ICER_llin_scenario)
Disc_Malaria_res[3,]<-c(Disc_Inc_cost_irs, Disc_DALY_averted_irs, Disc_ICER_irs)
Disc_Malaria_res[4,]<-c(Disc_Inc_cost_new_irs_scenario, Disc_DALY_averted_new_irs_scenario, Disc_ICER_new_irs_scenario)
colnames(Disc_Malaria_res)<-c("Incremental Cost", "DALYs_averted", "ICER")
rownames(Disc_Malaria_res)<-c("LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", 
                         "LLIN(75%)+IRS(50%)")
Disc_Malaria_res

# Plot to see results on the ICER plane
plot((Disc_Malaria_res[,2])/1000000, (Disc_Malaria_res[,1])/1000000, xlim=c(-50,120), ylim=c(-240,200), ylab="Incremental costs (Million US$)", 
     xlab="DALY Averted (in millions)", pch=18, col=c("red", "brown", "blue", "black", "darkgrey"), main="CEA of Interventions for Malaria", cex.main=0.8)
#text((Malaria_res[,2])/1000000, (Malaria_res[,1])/1000000,label=c("LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", "LLIN(75%)+IRS(50%)"), cex=0.7, pos=2)
legend("topright", legend=c("LLIN only (95%)", "LLIN only (75%)", "LLIN(95%)+IRS(50%)", 
                            "LLIN(75%)+IRS(50%)"), fill = c("brown", "blue", "black", "darkgrey"), cex=0.7)

#abline(a=0, b=1969, lty=3, col='black', lwd=2) #b = GDP per capita of Nigeria in 2017
abline(v = 0, h = 0, lty=1, col="blue",lwd=1)
