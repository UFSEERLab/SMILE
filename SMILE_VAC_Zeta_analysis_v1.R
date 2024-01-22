# Testing the SMILE functions
# Load SMILE source code before running this script.

#VACCINE DATA IN LAO CAI PROVINCE (1991-2013) for 54 high-risk communes only

years <- 23

vacc.pcts<- c(0.34, 0.24, 0.37, 0.44, 0.44, 0.64, 0.71, 0.69, 0.73, 
              0.72, 0.66, 0.60, 0.57, 0.27, 0.27, 0.16, 0.12, 0.12, 
              0.10, 0.08, 0.05, 0.05, 0.02)

#VACCINE DATA IN LAO CAI PROVINCE (1991-2008) for 54 high-risk communes only
#assume vaccine coverage in 1990 is 0
years <- 19
vacc.pcts<- c(0.34, 0.24, 0.37, 0.44, 0.44, 0.64, 0.71, 0.69, 0.73, 
              0.72, 0.66, 0.60, 0.57, 0.27, 0.27, 0.16, 0.12, 0.12, 0.10)

zeta0.graph <- 1/(1+exp(-(beta0+beta1*vacc.pcts)))
zeta.graph <- rep(zeta0.graph,each=52)
zeta.graph <- c(zeta.graph[1],zeta.graph)

#target.zeta <- 0.88
#target.vacc <- (log(target.zeta)-log(1-target.zeta) -beta0)/beta1
#vacc.vec <- rep(target.vacc,10)#c(.2,.5,.3,.29,.34,.33,.19,.68,.70,.15)

####################################################################################
# Case 1: No age structuring, infection probability is not seasonal
N1 <- 31000
tau=0.88
b=0.001
theta=1000

beta0=-3.5
beta1=7.4

sim1.0	<-	smile.vacc1(b=b, theta=theta,tau=tau,years=years,beta0=beta0,
                      beta1=beta1,vacc.vec=vacc.vec,N1=N1)

n.weeks <- years*52+1
weeks <- 1:n.weeks
weekm1 <- 2:n.weeks
par(mfrow=c(2,3), mar=c(4,5,2,2), oma=c(2,3,1,1))
plot(x=weekm1, y=sim1.0[,1], type="l", xlab="Week", ylab="Susceptible", cex.lab=1.5)
plot(x=weekm1, y=sim1.0[,3], type="l", xlab="Week", ylab="Infected", cex.lab=1.5)
plot(x=weekm1, y=sim1.0[,2], type="l", xlab="Week", ylab="Immune", cex.lab=1.5)
plot(x=weekm1, y=sim1.0[,4], type="l", xlab="Week", ylab="LIZ", cex.lab=1.5)
plot(x=weekm1, y=sim1.0[,5], type="l", xlab="Week", ylab="Environment", cex.lab=1.5)
plot(x=weeks, y=zeta.graph, type="l", xlab="Week", ylab="Zeta by vaccine", cex.lab=1.5)


#Total estimated dead buffalo from 1991-2009
cases.integ1 <- as.integer(sim1.0[,4])
sum(cases.integ1)
#Total estimated dead buffalo from 1991-2000
sum(cases.integ1[2:521])
#Total estimated dead buffalo from 2001-2009
sum(cases.integ1[522:988])

par(mfrow=c(1,1), mar=c(4,5,2,2), oma=c(2,3,1,1))
plot(x=weekm1, y=sim1.0[,4], type="l", xlab="Week", ylim=c(1, 2), ylab="LIZ", cex.lab=1.5)
plot(x=weekm1, y=cases.integ1, type="l")

