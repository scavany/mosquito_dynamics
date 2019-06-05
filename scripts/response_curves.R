
load('algam_85re.Rdata')


#### 1st biting rate  (Focks et al., 1993)
biterate.1st <- function (t) { # biting rate assuming is 1/fgc
  rho <- 0.216
  HA <-  15725.23
  HH <- 1756481.07
  T05 <- 447.17
  Tk<- t + 273.15
  gc<- rho*(Tk/298)* exp((HA/1.987)*((1/298)-(1/Tk))) /
    (1+exp((HH/1.987)*((1/T05)-(1/Tk))))
  return (gc)
}

#### 2nd biting rate  (Otero et al., 2006)
biterate.2nd <- function (t) { 
  rho <- 0.372
  HA <-  15725.23
  HH <- 1756481.07
  T05 <- 447.17
  Tk<- t + 273.15
  gc<- rho*(Tk/298)* exp((HA/1.987)*((1/298)-(1/Tk))) /
    (1+exp((HH/1.987)*((1/T05)-(1/Tk))))
  return (gc)
}
########### Mean of any given distrbution
meandist<- function(distx,xvec) {
  p=distx/sum(distx)
  sum(xvec*p)
}


iploc <- function(temp,b0,bt) 
  exp(exp(b0 + bt*temp))

### EIP as a function of temperature (Chan & Johanesson, 2012)
### extrinstic incubation period distribution
eip<- function(temp,tmax=xwidth,tb=tby){
  x = seq(0,tmax,tb)
  tau = 4.9
  b0 = 2.9
  bt = -0.08
  mu = iploc(temp,b0,bt)
  eipdf = dlnorm(x, meanlog=log(mu), sdlog=(1/sqrt(tau)))
  epidf_n = eipdf / sum(eipdf)
  mean.eip = meandist(eipdf,x)
  return(list(eipdf_n = epidf_n, mu = mu, tau = tau, mean.eip = mean.eip))
}  

## distribution of mortality rate (Brady etal,2013)
mortalityRT<- function(temp,fieldcorxn) {
  dd<-seq(0,120,length.out=(120*24+2))
  nwdd<-  data.frame(Days=dd,Temperature=rep(temp, (120*24+2)), Study_number=5, Feed_B=2, Feed_S=1)
  nwdd<-  cbind(nwdd,logDay=log(nwdd$Days+1), logTemp=log(nwdd$Temperature+1))  ## +1 avoids log(0) 
  prediction <-as.vector(unlist(predict(algam,newdata = (nwdd), se.fit = TRUE, type = "response")$fit))
  prediction<- prediction[-1]
  prediction[1:24]<- prediction[1:24]/prediction[1]
  prediction[which(prediction>1)]<-1
  prediction[which(prediction<=0.001)]<-0
  diffDeath<- -diff(prediction)
  diffDeath<- diffDeath/sum(diffDeath)
  return(1/sum(dd[2:length(prediction)]*diffDeath) + fieldcorxn) 
}  ##### end function 




