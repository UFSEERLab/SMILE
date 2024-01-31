# Testing the SMILE functions
# Load SMILE source code before running this script.

#VACCINE DATA IN LAO CAI PROVINCE (1991-2013) for 54 high-risk communes only

years <- 23
vacc.pcts<- c(0.34, 0.24, 0.37, 0.44, 0.44, 0.64, 0.71, 0.69, 0.73, 
              0.72, 0.66, 0.60, 0.57, 0.27, 0.27, 0.16, 0.12, 0.12, 
              0.10, 0.08, 0.05, 0.05, 0.02)

#VACCINE DATA IN LAO CAI PROVINCE (1991-2008) for 54 high-risk communes only
#assume vaccine coverage in 1990 is 0
#years <- 19
#vacc.pcts<- c(0.34, 0.24, 0.37, 0.44, 0.44, 0.64, 0.71, 0.69, 0.73, 
#              0.72, 0.66, 0.60, 0.57, 0.27, 0.27, 0.16, 0.12, 0.12, 0.10)

# Create a vector of vaccine-modified Zeta for the final graph
zeta0.graph <- 1/(1+exp(-(beta0+beta1*vacc.pcts)))
zeta.graph <- rep(zeta0.graph,each=52)
zeta.graph <- c(zeta.graph[1],zeta.graph)

#target.zeta <- 0.88
#target.vacc <- (log(target.zeta)-log(1-target.zeta) -beta0)/beta1
#vacc.vec <- rep(target.vacc,10)#c(.2,.5,.3,.29,.34,.33,.19,.68,.70,.15)

####################################################################################
# Case 1: No age structuring, infection probability is not seasonal, vaccine coverage
N1 <- 20000
tau=1
theta=100
b=0.001
beta0=-3.5
beta1=7.4
years <- 23

sim1.vax	<-	smile.vacc1(b=b, theta=theta,tau=tau,years=years,beta0=beta0,
                      beta1=beta1,vacc.vec=vacc.pcts,N1=N1)

n.weeks <- years*52+1
weeks <- 1:n.weeks
weekm1 <- 2:n.weeks
par(mfrow=c(2,3), mar=c(4,5,2,2), oma=c(2,3,1,1))
plot(x=weekm1, y=sim1.vax[,1], type="l", xlab="Week", ylab="Susceptible", cex.lab=1.5)
plot(x=weekm1, y=sim1.vax[,3], type="l", xlab="Week", ylab="Infected", cex.lab=1.5)
plot(x=weekm1, y=sim1.vax[,2], type="l", xlab="Week", ylab="Immune", cex.lab=1.5)
plot(x=weekm1, y=sim1.vax[,4], type="l", xlab="Week", ylab="LIZ", cex.lab=1.5)
plot(x=weekm1, y=sim1.vax[,5], type="l", xlab="Week", ylab="Environment", cex.lab=1.5)
plot(x=weeks, y=zeta.graph, type="l", xlab="Week", ylab="Zeta by vaccine", cex.lab=1.5)


#Total estimated dead buffalo from 1991-2009
cases.integ1 <- as.integer(sim1.vax[,4])
sum(cases.integ1)
#Total estimated dead buffalo from 1991-2000
sum(cases.integ1[2:521])
#Total estimated dead buffalo from 2001-2009
sum(cases.integ1[522:988])

############################################################################
# Case 2: No age structuring, infection probability is seasonal, vaccine coverage
N1 <- 20000
tau=1
theta=100
b0=-30
b1=0.85
period=3*52
beta0=-3.5
beta1=7.4
years <- 23

sim2.vax	<-	smile.vacc2(b0=b0,b1=b1,period=period,theta=theta,tau=tau,years=years,
                        beta0=beta0,beta1=beta1,vacc.vec=vacc.pcts,N1=N1)

n.weeks <- years*52+1
weeks <- 1:n.weeks
weekm1 <- 2:n.weeks
par(mfrow=c(2,3), mar=c(4,5,2,2), oma=c(2,3,1,1))
plot(x=weekm1, y=sim2.vax[,1], type="l", xlab="Week", ylab="Susceptible", cex.lab=1.5)
plot(x=weekm1, y=sim2.vax[,3], type="l", xlab="Week", ylab="Infected", cex.lab=1.5)
plot(x=weekm1, y=sim2.vax[,2], type="l", xlab="Week", ylab="Immune", cex.lab=1.5)
plot(x=weekm1, y=sim2.vax[,4], type="l", xlab="Week", ylab="LIZ", cex.lab=1.5)
plot(x=weekm1, y=sim2.vax[,5], type="l", xlab="Week", ylab="Environment", cex.lab=1.5)
plot(x=weeks, y=zeta.graph, type="l", xlab="Week", ylab="Zeta by vaccine", cex.lab=1.5)


#Total estimated dead buffalo from 1991-2009
cases.integ2 <- as.integer(sim2.vax[,4])
sum(cases.integ2)
#Total estimated dead buffalo from 1991-2000
sum(cases.integ2[2:521])
#Total estimated dead buffalo from 2001-2009
sum(cases.integ2[522:988])

#####################################################################################
# Case 3: Age structuring, vaccine coverage, infection probability is NOT seasonal

N1 <- 20000
tau=1
theta=100
b=0.001
beta0=-3.5
beta1=7.4
years <- 23

sim3.vax	<-	smile.vacc3(b=b,theta=theta,tau=tau,years=years,
                        beta0=beta0,beta1=beta1,vacc.vec=vacc.pcts,N1=N1)

n.weeks <- years*52+1
weeks <- 1:n.weeks
weekm1 <- 2:n.weeks
par(mfrow=c(2,3), mar=c(4,5,2,2), oma=c(2,3,1,1))
plot(x=weekm1, y=sim3.vax[,1], type="l", xlab="Week", ylab="Susceptible", cex.lab=1.5)
plot(x=weekm1, y=sim3.vax[,3], type="l", xlab="Week", ylab="Infected", cex.lab=1.5)
plot(x=weekm1, y=sim3.vax[,2], type="l", xlab="Week", ylab="Immune", cex.lab=1.5)
plot(x=weekm1, y=sim3.vax[,4], type="l", xlab="Week", ylab="LIZ", cex.lab=1.5)
plot(x=weekm1, y=sim3.vax[,5], type="l", xlab="Week", ylab="Environment", cex.lab=1.5)
plot(x=weeks, y=zeta.graph, type="l", xlab="Week", ylab="Zeta by vaccine", cex.lab=1.5)


#Total estimated dead buffalo from 1991-2009
cases.integ3 <- as.integer(sim3.vax[,4])
sum(cases.integ3)
#Total estimated dead buffalo from 1991-2000
sum(cases.integ3[2:521])
#Total estimated dead buffalo from 2001-2009
sum(cases.integ3[522:988])

#####################################################################################
# Case 4: Age structuring, infection probability has seasonal forcing, vaccine coverage
N1 <- 20000
tau=1
theta=100
b0=-30
b1=0.85
period=3*52
beta0=-3.5
beta1=7.4
years <- 23

sim4.vax	<-	smile.vacc4(b0=b0,b1=b1,period=period,theta=theta,tau=tau,
                        years=years,beta0=beta0,beta1=beta1,vacc.vec=vacc.pcts,N1=N1)

n.weeks <- years*52+1
weeks <- 1:n.weeks
weekm1 <- 2:n.weeks
par(mfrow=c(2,3), mar=c(4,5,2,2), oma=c(2,3,1,1))
plot(x=weekm1, y=sim4.vax[,1], type="l", xlab="Week", ylab="Susceptible", cex.lab=1.5)
plot(x=weekm1, y=sim4.vax[,3], type="l", xlab="Week", ylab="Infected", cex.lab=1.5)
plot(x=weekm1, y=sim4.vax[,2], type="l", xlab="Week", ylab="Immune", cex.lab=1.5)
plot(x=weekm1, y=sim4.vax[,4], type="l", xlab="Week", ylab="LIZ", cex.lab=1.5)
plot(x=weekm1, y=sim4.vax[,5], type="l", xlab="Week", ylab="Environment", cex.lab=1.5)
plot(x=weeks, y=zeta.graph, type="l", xlab="Week", ylab="Zeta by vaccine", cex.lab=1.5)


#Total estimated dead buffalo from 1991-2009
cases.integ4 <- as.integer(sim4.vax[,4])
sum(cases.integ4)
#Total estimated dead buffalo from 1991-2000
sum(cases.integ4[2:521])
#Total estimated dead buffalo from 2001-2009
sum(cases.integ4[522:988])


#####################################################################################
# Case 5: Age structuring, infection probability is seasonal, 
#         survival from other causes, vaccine coverage
N1 <- 20000
tau=1
theta=100
b0=-30
b1=0.85
period=3*52
beta0=-3.5
beta1=7.4
years <- 23

sim5.vax	<-	smile.vacc5(b0=b0,b1=b1,period=period,theta=theta,tau=tau,
                        years=years,beta0=beta0,beta1=beta1,vacc.vec=vacc.pcts,N1=N1)

n.weeks <- years*52+1
weeks <- 1:n.weeks
weekm1 <- 2:n.weeks
par(mfrow=c(2,3), mar=c(4,5,2,2), oma=c(2,3,1,1))
plot(x=weekm1, y=sim5.vax[,1], type="l", xlab="Week", ylab="Susceptible", cex.lab=1.5)
plot(x=weekm1, y=sim5.vax[,3], type="l", xlab="Week", ylab="Infected", cex.lab=1.5)
plot(x=weekm1, y=sim5.vax[,2], type="l", xlab="Week", ylab="Immune", cex.lab=1.5)
plot(x=weekm1, y=sim5.vax[,4], type="l", xlab="Week", ylab="LIZ", cex.lab=1.5)
plot(x=weekm1, y=sim5.vax[,5], type="l", xlab="Week", ylab="Environment", cex.lab=1.5)
plot(x=weeks, y=zeta.graph, type="l", xlab="Week", ylab="Zeta by vaccine", cex.lab=1.5)


#Total estimated dead buffalo from 1991-2009
cases.integ5 <- as.integer(sim5.vax[,4])
sum(cases.integ5)
#Total estimated dead buffalo from 1991-2000
sum(cases.integ5[2:521])
#Total estimated dead buffalo from 2001-2009
sum(cases.integ5[522:988])

#####################################################################################
# Case 6: Age structuring, survival from other causes, vaccine coverage. 
#                                     Infection probability is seasonal.
N1 <- 20000
tau=1
theta=100
b=0.001
beta0=-3.5
beta1=7.4
years <- 23

sim6.vax	<-	smile.vacc6(b=b,theta=theta,tau=tau,years=years,beta0=beta0,
                        beta1=beta1,vacc.vec=vacc.pcts,N1=N1)

n.weeks <- years*52+1
weeks <- 1:n.weeks
weekm1 <- 2:n.weeks
par(mfrow=c(2,3), mar=c(4,5,2,2), oma=c(2,3,1,1))
plot(x=weekm1, y=sim6.vax[,1], type="l", xlab="Week", ylab="Susceptible", cex.lab=1.5)
plot(x=weekm1, y=sim6.vax[,3], type="l", xlab="Week", ylab="Infected", cex.lab=1.5)
plot(x=weekm1, y=sim6.vax[,2], type="l", xlab="Week", ylab="Immune", cex.lab=1.5)
plot(x=weekm1, y=sim6.vax[,4], type="l", xlab="Week", ylab="LIZ", cex.lab=1.5)
plot(x=weekm1, y=sim6.vax[,5], type="l", xlab="Week", ylab="Environment", cex.lab=1.5)
plot(x=weeks, y=zeta.graph, type="l", xlab="Week", ylab="Zeta by vaccine", cex.lab=1.5)


#Total estimated dead buffalo from 1991-2009
cases.integ6 <- as.integer(sim6.vax[,4])
sum(cases.integ6)
#Total estimated dead buffalo from 1991-2000
sum(cases.integ6[2:521])
#Total estimated dead buffalo from 2001-2009
sum(cases.integ6[522:988])
