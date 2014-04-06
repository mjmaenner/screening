# Function to calculate prevalence and variance from two-phase screening designs
#
# Formula originally presented in:
# Shrout PE and Newman SC, "Design of Two-Phase Prevalence Surveys of Rare Disorders"
# Biometrics ,1989; 45; 549-555
#
#
# Explanation of input variables:
#   Nsamp: the total number of persons screened
#   ppv: probability of disease given a positive screen
#   fnpv: probability of having disease given negative screen
#   PrS: Proportion of sample screening positive
#   fuS: Proportion of second stage follow-up for positive screens
#   fuSnot: Proportion of second-stage follow-up for negative screens

twoPhase<-function(Nsamp, ppv, fnpv, PrS, fuS, fuSNot)
  {
  PrSNot <- 1- PrS
  twoPhase_prev<-ppv*PrS+fnpv*PrSNot
  twoPhase_var<-(1/Nsamp)*((PrS*ppv*(1-ppv))/fuS  +  (PrSNot*fnpv*(1-fnpv))/fuSNot  +  PrS*PrSNot*(ppv-fnpv)^2)
  twoPhase_l95<-twoPhase_prev - 1.96 * sqrt(twoPhase_var)
  twoPhase_u95<-twoPhase_prev + 1.96 * sqrt(twoPhase_var)
return(cbind(twoPhase_prev, twoPhase_var, twoPhase_l95, twoPhase_u95)) 
}


# example usage :
# twoPhase(Nsamp=200,ppv=1,fnpv=0, PrS=.20,fuS=1,fuSNot=1)

# The above scenario would be equivalent to a simple (1-stage) sample:
# require(Hmisc)
# binconf(40, 200, method="asymptotic")
