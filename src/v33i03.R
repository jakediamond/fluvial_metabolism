###################################################
### chunk number 1: preliminaries
###################################################
library("FME")


###################################################
### The HIV model, as an R-function
###################################################
HIV_R <- function (pars, V_0 = 50000, dV_0 = -200750, T_0 = 100) {

  derivs <- function( time, y, pars) {
    with (as.list(c(pars,y)), {

      dT = lam - rho*T - bet*T*V
      dI = bet*T*V - delt*I
      dV = n*delt*I - c*V - bet*T*V

      return(list(c(dT, dI, dV), logV = log(V)))
    })
  }

  # initial conditions
  I_0   <- with(as.list(pars), (dV_0 + c*V_0) / (n*delt))
  y     <- c(T = T_0,  I = I_0,  V = V_0)

  times <- c(seq(0, 0.8, 0.1), seq(2, 60, 2))
  out   <- ode(y=y,parms=pars,times=times,func=derivs)

  as.data.frame(out)
}


###################################################
### The HIV model, as a DLL 
###################################################
HIV <- function (pars, V_0 = 50000, dV_0 = -200750, T_0 = 100) {

  I_0 <- with(as.list(pars),(dV_0 + c*V_0) / (n*delt))
  y <- c(T = T_0,  I = I_0,  V = V_0)

  times <- c(0,0.1,0.2,0.4,0.6,0.8,seq(2,60,by=2))
  out <- ode(y = y, parms = pars, times = times, func = "derivshiv",
    initfunc = "inithiv", nout = 1, outnames = c("logV"), dllname = "FME")

  as.data.frame(out)
}

pars = c(bet=0.00002, rho=0.15, delt=0.55, c=5.5, lam=80, n=900)
out <- HIV(pars=pars)


par(mfrow=c(1, 2))
plot(out$time, out$logV, main="Viral load", ylab="log(V)",
  xlab="time", type="b")
plot(out$time, out$T, main="CD4+ T", ylab="-", xlab="time", type="b")
par(mfrow=c(1,1))

###################################################
### The pseudodata for log(viral count), T-cells
###################################################

DataLogV  <- cbind(time=out$time,
             logV = out$logV + rnorm(sd=0.45,n=length(out$logV)),
             sd = 0.45)

ii    <- which (out$time %in% seq(0, 56, by=4))
DataT <- cbind(time=out$time[ii],
               T = out$T[ii] + rnorm(sd=4.5,n=length(ii)),
               sd =4.5)
head(DataT)


###################################################
### The "cost" function 
###################################################

HIVcost <- function (pars) {
  out <- HIV(pars)
  cost <- modCost(model=out, obs=DataLogV, err="sd")
  return( modCost(model=out, obs=DataT   , err="sd", cost=cost))
}

HIVcost(pars)$model
plot(HIVcost(pars))


###################################################
### Sensitivity functions 
###################################################

Sfun <- sensFun(HIVcost,pars)
summary(Sfun)

plot(Sfun, which = c("logV","T"), lwd=2)
pairs(Sfun,which=c("logV","T"),col=c("blue","green"))

###################################################
### Model identifiability (collinearity)
###################################################

ident <- collin(Sfun)
head(ident,n=20)
plot(ident,log="y")

collin(Sfun,parset=c("bet","rho","delt","c","lam","n"))
collin(Sfun,N=5)
collin(Sfun,parset=c("bet","rho","delt","c","lam"),which="logV")
collin(Sfun,parset=c("bet","rho","delt","c","lam"),which="T")


###################################################
### Cost function with 5 parameters; log-transformed
###################################################
HIVcost2 <- function(lpars)
    HIVcost(c(exp(lpars),n=900))


###################################################
### Fit model to data
###################################################

Pars <- pars[1:5]*2  # perturb original pars
Fit  <- modFit(f=HIVcost2, p=log(Pars))
exp(coef(Fit))
deviance(Fit)


###################################################
### Show initial and best fit
###################################################
ini   <- HIV(pars=c(Pars,n=900))
final <- HIV(pars=c(exp(coef(Fit)),n=900))

par(mfrow=c(1,2))
plot(DataLogV, xlab="time", ylab="logV", ylim = c(7,11))
lines(ini$time,ini$logV, lty=2)
lines(final$time,final$logV)
legend("topright", c("data","initial","fitted"),
   lty=c(NA,2,1), pch=c(1,NA,NA))
plot(DataT, xlab="time", ylab="T")
lines(ini$time,ini$T, lty=2)
lines(final$time,final$T)
par(mfrow=c(1,1))

###################################################
### Run MCMC
###################################################
var0 <- Fit$var_ms_unweighted
cov0 <- summary(Fit)$cov.scaled*2.4^2/5

# this takes a while (a minute) ...
MCMC <- modMCMC(f=HIVcost2, p=Fit$par, niter=5000, jump=cov0,
                var0=var0, wvar0=0.1, updatecov=50)

MCMC$pars <- exp(MCMC$pars)
summary(MCMC)

plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000)

###################################################
### Sensitivity ranges
###################################################
sR <- sensRange(func = HIV, parms = pars, parInput = MCMC$par)
plot(summary(sR), xlab = "time")


###################################################
### Monte carlo
###################################################
parRange <- cbind(min=0.75*pars,max=1.25*pars)
crlfun <- function (pars)  return(c(meanVirus = mean(HIV(pars)$V)))

CRL <- modCRL(fun=crlfun, parRange = parRange, num=500)

cor(CRL)[7,]
plot(CRL, ylab="number of virions", trace=TRUE)

