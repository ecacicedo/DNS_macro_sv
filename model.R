# Macro-Finance yield Curve Modeling with Stochastic Volatility
# Author: ecacicedo (expanded from Prof. Laurini's)
# 21/05/2012

# Set up
install.packages('forecast')
install.packages('R2WinBUGS')
library(forecast)
library(R2WinBUGS)

setwd("./path/to/")
jurosmacro=read.table("data.txt",header=T)


## Dynamic Nelson-Siegel (DNS yields only) Model by Diebold & Li
datas=1:length(jurosmacro)

juros=jurosmacro[1:dim(jurosmacro)[1],4:17]/100

n=dim(juros)[1]

vert=c(1, 2, 3, 4,6, 9, 12,15,18,24,27,29,31,33)
dlyields=as.matrix(juros)

b1=double(length(1:dim(juros)[1]))
b2=double(length(1:dim(juros)[1]))
b3=double(length(1:dim(juros)[1]))

for (j in 1:dim(juros)[1]){
  
  n=dim(juros)[1]
  k=dim(juros)[2]
  m <- n + 12
  nd=n*k
  md=m*k
  
  # Neslson-Siegel equation
  
  lambda=.0609
  
  f1=seq(1,1,length=length(vert))
  f2=(1-exp(-lambda*vert))/(lambda*vert)
  f3=((1-exp(-lambda*vert))/(lambda*vert)-exp(-lambda*vert))
  
  # First step: OLS
  dlreg=lm(dlyields[j,]~f1+f2+f3-1)
  b1[j]=dlreg$coef[1]
  b2[j]=dlreg$coef[2]
  b3[j]=dlreg$coef[3]
}

jpeg("betashatdl.jpeg")
par(mfrow=c(3,1))
ts.plot(ts(b1[1:dim(juros)[1]],frequency=12,start = c(2003, 1)),type="l",ylab=expression(hat(beta[1*t])))
ts.plot(ts(b2[1:dim(juros)[1]],frequency=12,start = c(2003, 1)),type="l",ylab=expression(hat(beta[2*t])))
ts.plot(ts(b3[1:dim(juros)[1]],frequency=12,start = c(2003, 1)),type="l",ylab=expression(hat(beta[3*t])))
dev.off()

# Second Step: AR(1)
arb1=ar.ols(b1[1:104], aic = FALSE, order.max = 1)
arb2=ar.ols(b2[1:104], aic = FALSE, order.max = 1)
arb3=ar.ols(b3[1:104], aic = FALSE, order.max = 1)


f11=as.matrix(f1)
f12=as.matrix(f2)
f13=as.matrix(f3)

# Model assessment
m=array(0,dim=c(dim(juros)[1],dim(juros)[2]))

for (kk in 1:n){
  for (i in 1:dim(juros)[2])
    m[kk,i]<- b1[kk]*f11[i,1] + b2[kk]*f12[i,1] + b3[kk]*f13[i,1]
  
}

model1v1=accuracy(m[1:n,1]*10000,dlyields[1:n,1]*10000)
model1v2=accuracy(m[1:n,2]*10000,dlyields[1:n,2]*10000)
model1v3=accuracy(m[1:n,3]*10000,dlyields[1:n,3]*10000)
model1v4=accuracy(m[1:n,4]*10000,dlyields[1:n,4]*10000)
model1v5=accuracy(m[1:n,5]*10000,dlyields[1:n,5]*10000)
model1v6=accuracy(m[1:n,6]*10000,dlyields[1:n,6]*10000)
model1v7=accuracy(m[1:n,7]*10000,dlyields[1:n,7]*10000)
model1v8=accuracy(m[1:n,8]*10000,dlyields[1:n,8]*10000)
model1v9=accuracy(m[1:n,9]*10000,dlyields[1:n,9]*10000)
model1v10=accuracy(m[1:n,10]*10000,dlyields[1:n,10]*10000)
model1v11=accuracy(m[1:n,11]*10000,dlyields[1:n,11]*10000)
model1v12=accuracy(m[1:n,12]*10000,dlyields[1:n,12]*10000)
model1v13=accuracy(m[1:n,13]*10000,dlyields[1:n,13]*10000)
model1v14=accuracy(m[1:n,14]*10000,dlyields[1:n,14]*10000)

comparacao=matrix(0,14,7)
compmatrix=as.data.frame(comparacao)   
colnames(compmatrix)=c("ME","SD","RMSE","MAE","MPE","MAPE","ACF1")
rownames(compmatrix)=c("1", "2", "3", "4", "6", "9", "12", "15", "18", "24", "27", "29", "31", "33")

compmatrix[1,c(1,3,4,5,6)]=    model1v1  
compmatrix[2,c(1,3,4,5,6)]=    model1v2    
compmatrix[3,c(1,3,4,5,6)]=    model1v3    
compmatrix[4,c(1,3,4,5,6)]=    model1v4    
compmatrix[5,c(1,3,4,5,6)]=    model1v5    
compmatrix[6,c(1,3,4,5,6)]=    model1v6    
compmatrix[7,c(1,3,4,5,6)]=    model1v7    
compmatrix[8,c(1,3,4,5,6)]=    model1v8    
compmatrix[9,c(1,3,4,5,6)]=    model1v9    
compmatrix[10,c(1,3,4,5,6)]=   model1v10    
compmatrix[11,c(1,3,4,5,6)]=   model1v11    
compmatrix[12,c(1,3,4,5,6)]=   model1v12    
compmatrix[13,c(1,3,4,5,6)]=   model1v13    
compmatrix[14,c(1,3,4,5,6)]=   model1v14

compmatrix[1,2]=  sd(m[1:n,1]*10000-dlyields[1:n,1]*10000)     
compmatrix[2,2]=  sd(m[1:n,2]*10000-dlyields[1:n,2]*10000)     
compmatrix[3,2]=  sd(m[1:n,3]*10000-dlyields[1:n,3]*10000)     
compmatrix[4,2]=  sd(m[1:n,4]*10000-dlyields[1:n,4]*10000)     
compmatrix[5,2]=  sd(m[1:n,5]*10000-dlyields[1:n,5]*10000)     
compmatrix[6,2]=  sd(m[1:n,6]*10000-dlyields[1:n,6]*10000)     
compmatrix[7,2]=  sd(m[1:n,7]*10000-dlyields[1:n,7]*10000)     
compmatrix[8,2]=  sd(m[1:n,8]*10000-dlyields[1:n,8]*10000)     
compmatrix[9,2]=  sd(m[1:n,9]*10000-dlyields[1:n,9]*10000)     
compmatrix[10,2]= sd(m[1:n,10]*10000-dlyields[1:n,10]*10000)
compmatrix[11,2]= sd(m[1:n,11]*10000-dlyields[1:n,11]*10000)
compmatrix[12,2]= sd(m[1:n,12]*10000-dlyields[1:n,12]*10000)
compmatrix[13,2]= sd(m[1:n,13]*10000-dlyields[1:n,13]*10000)
compmatrix[14,2]= sd(m[1:n,14]*10000-dlyields[1:n,14]*10000)

compmatrix[1,7]=  acf(m[1:n,1]*10000-dlyields[1:n,1]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2]      
compmatrix[2,7]=  acf(m[1:n,2]*10000-dlyields[1:n,2]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2]      
compmatrix[3,7]=  acf(m[1:n,3]*10000-dlyields[1:n,3]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
compmatrix[4,7]=  acf(m[1:n,4]*10000-dlyields[1:n,4]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2]      
compmatrix[5,7]=  acf(m[1:n,5]*10000-dlyields[1:n,5]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
compmatrix[6,7]=  acf(m[1:n,6]*10000-dlyields[1:n,6]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
compmatrix[7,7]=  acf(m[1:n,7]*10000-dlyields[1:n,7]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
compmatrix[8,7]=  acf(m[1:n,8]*10000-dlyields[1:n,8]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
compmatrix[9,7]=  acf(m[1:n,9]*10000-dlyields[1:n,9]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
compmatrix[10,7]= acf(m[1:n,10]*10000-dlyields[1:n,10]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
compmatrix[11,7]= acf(m[1:n,11]*10000-dlyields[1:n,11]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
compmatrix[12,7]= acf(m[1:n,12]*10000-dlyields[1:n,12]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
compmatrix[13,7]= acf(m[1:n,13]*10000-dlyields[1:n,13]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
compmatrix[14,7]= acf(m[1:n,14]*10000-dlyields[1:n,14]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 

compmatrix


### DNS macro-sv model

## Markov-Chain Monte Carlo

outcome.type=1

n.thin = 1

n.iter = 60000
n.burnin =50000


# Maturity matrix
verticesm=matrix(0,dim(juros)[1],14)
for (j in 1:dim(juros)[1]){
  verticesm[j,]=vert
}

# Naming macroeconomic variables
CU=jurosmacro[,2]/100
PI=jurosmacro[,1]/100
FFR=jurosmacro[,3]/100

# Total observations
dias=dim(juros)[1]

counters=double(dias)


for (j in 1:dim(verticesm)[1]){
  counters[j]=dim(verticesm)[2]
}

TM=dim(verticesm)[1]


Maturity=(verticesm)
yields=as.matrix(juros)


inits.beta=rep(0.7,5*TM)

inits<-function(){list(beta=structure(.Data=inits.beta,.Dim = c(dias
                                                                ,6)),taueps=80,itau2=80)}
n=20

# Initial values 42 x 42 diagonal amtrix
pfdata=matrix(data=0,42,42)
for (p in 1:42){
  pfdata[p,p]=10.10
}
R=pfdata


# Initial vales 1 X 104 vector
betap=rep(.01,dias)


# Initial values 42 x 42 diagonal amtrix
Omega=matrix(data=0,42,42)
for (p in 1:42){
  Omega[p,p]=.001
}

# Initial values 1 X 104 vector
betap1=b1
betap2=b2
betap3=b3
betap1[1]=NA
betap2[1]=NA
betap3[1]=NA

svpp1=rep(-3.89,dias)
svpp2=rep(-3.99,dias)
svpp3=rep(-3.09,dias)
svpp4=rep(-3.89,dias)
svpp5=rep(-3.99,dias)
svpp6=rep(-3.09,dias)
svpp=rep(-3.19,dias)

# Initial values 1 X 104 vector
lambdap=rep(-3.58522,dias)
lambdap[1]=NA


imphi=double(42)
imphi[1]<-.3 
imphi[2]<--.01 
imphi[3]<-.04  
imphi[4]<-.6  
imphi[5]<-.065  
imphi[6]<-0.041 
imphi[7]<-0.9  
imphi[8]<-0  
imphi[9]<-0  
imphi[10]<-0 
imphi[11]<-0 
imphi[12]<-0 
imphi[13]<-0 
imphi[14]<-0.9 
imphi[15]<-0 
imphi[16]<-0 
imphi[17]<-0 
imphi[18]<-0 
imphi[19]<-0 
imphi[20]<-0 
imphi[21]<-0.9 
imphi[22]<-0 
imphi[23]<-0 
imphi[24]<-0 
imphi[25]<-0 
imphi[26]<-0 
imphi[27]<-0 
imphi[28]<-0.9
imphi[29]<-0 
imphi[30]<-0 
imphi[31]<-0 
imphi[32]<-0 
imphi[33]<-0 
imphi[34]<-0 
imphi[35]<-0.9
imphi[36]<-0 
imphi[37]<-0 
imphi[38]<-0 
imphi[39]<-0 
imphi[40]<-0 
imphi[41]<-0 
imphi[42]<-0.9

# Run model from WinBUGS
inits<-function(){list(itau2=1,sv0=2,mphi=imphi,
                       beta1=betap1,beta2=betap2,beta3=betap3,theta=svpp,theta1=svpp1,theta2=svpp2,theta3=svpp3,plambda1=lambdap,theta4=svpp4,theta5=svpp5,theta6=svpp6,varphi=Omega)}


Tfinal=120
t=vert
Maturity=Maturity


# Create object with results from WinBUGS estimation
data<-list("TM","yields","Maturity","counters","CU","PI","FFR","R")


# model to be estimated on WinBUGS
program.file.name="bugsmodel.txt"


# Parameter list
parameters<-list("beta1","beta2","beta3","svsigma2","svsigma21","svsigma22","svsigma23","phi1","sv0","sv1","sv01","sv11", "sv02","sv12","sv03","sv13","lambda1","phil11","phil12","sv04","sv14", "sv05","sv15","sv06","sv16","svsigmaCU","svsigmaPI","svsigmaFFR","varphi")


# Run model from WinBUGS
Bayes.fit<- bugs(data, inits, parameters, model.file = program.file.name,n.chains = 1, n.iter = n.iter, n.burnin = n.burnin,n.thin = n.thin,debug = F, DIC = TRUE, digits = 4, codaPkg =FALSE, bugs.directory = "./path/to")

attach.bugs(Bayes.fit)

# Estimation plots	   
jpeg("yields.jpeg")
ts.plot(yields,col=seq(1,17,by=1))
legend(150,.12,vert,col=seq(1,17,by=1),lty=1)
dev.off()


beta1m=matrix(NA,TM)
beta1up=matrix(NA,TM)
beta1ul=matrix(NA,TM)
for (j in 1:TM){
  beta1m[j]=quantile(beta1[,j],.5)
  beta1up[j]=quantile(beta1[,j],.025)
  beta1ul[j]=quantile(beta1[,j],.975)
}


beta2m=matrix(NA,TM)                          
beta2up=matrix(NA,TM)                         
beta2ul=matrix(NA,TM)                         
for (j in 1:TM){                              
  beta2m[j]=quantile(beta2[,j],.5)              
  beta2up[j]=quantile(beta2[,j],.025)           
  beta2ul[j]=quantile(beta2[,j],.975)           
}                                             


beta3m=matrix(NA,TM)                          
beta3up=matrix(NA,TM)                         
beta3ul=matrix(NA,TM)                         
for (j in 1:TM){                              
  beta3m[j]=quantile(beta3[,j],.5)              
  beta3up[j]=quantile(beta3[,j],.025)           
  beta3ul[j]=quantile(beta3[,j],.975)           
}                                             

betam=matrix(NA,TM,3)
betaup=matrix(NA,TM,3)
betaul=matrix(NA,TM,3)


betam[,1]=beta1m
betaup[,1]=beta1up
betaul[,1]=beta1ul

betam[,2]=beta2m
betaup[,2]=beta2up
betaul[,2]=beta2ul

betam[,3]=beta3m
betaup[,3]=beta3up
betaul[,3]=beta3ul


jpeg("beta1.jpeg")
plot(betam[,1],type="l",xlab="Tempo",ylab=expression(beta[1]))
lines(betaup[,1],lty=3)
lines(betaul[,1],lty=3)
title(expression(beta[1]))
dev.off()

jpeg("beta2.jpeg")
plot(betam[,2],type="l",xlab="Tempo",ylab=expression(beta[2]))
lines(betaup[,2],lty=3)
lines(betaul[,2],lty=3)
title(expression(beta[2]))
dev.off()

jpeg("beta3.jpeg")
plot(betam[,3],type="l",xlab="Tempo",ylab=expression(beta[3]))
lines(betaup[,3],lty=3)
lines(betaul[,3],lty=3)
title(expression(beta[3]))
dev.off()


jpeg("beta1c.jpeg")
plot(betam[,1],type="l",xlab="Tempo",ylab=expression(beta[1]))
lines(b1,lty=3)
dev.off()

jpeg("beta2c.jpeg")
plot(betam[,2],type="l",xlab="Tempo",ylab=expression(beta[2]))
lines(b2,lty=3)
dev.off()

jpeg("beta3c.jpeg")
plot(betam[,3],type="l",xlab="Tempo",ylab=expression(beta[3]))
lines(b3,lty=3)
dev.off()



svsigma2m=matrix(NA,TM)
svsigma2up=matrix(NA,TM)
svsigma2ul=matrix(NA,TM)
for (j in 1:TM){
  svsigma2m[j]=quantile(1/svsigma2[2,j],.5)
  svsigma2up[j]=quantile(1/svsigma2[,j],.025)
  svsigma2ul[j]=quantile(1/svsigma2[,j],.975)
}


svsigma21m=matrix(NA,TM)                          
svsigma21up=matrix(NA,TM)                         
svsigma21ul=matrix(NA,TM)                         
for (j in 1:TM){                              
  svsigma21m[j]=quantile(1/svsigma21[,j],.5)              
  svsigma21up[j]=quantile(1/svsigma21[,j],.025)           
  svsigma21ul[j]=quantile(1/svsigma21[,j],.975)           
}                                             


svsigma22m=matrix(NA,TM)                          
svsigma22up=matrix(NA,TM)                         
svsigma22ul=matrix(NA,TM)                         
for (j in 1:TM){                              
  svsigma22m[j]=quantile(1/svsigma22[,j],.5)              
  svsigma22up[j]=quantile(1/svsigma22[,j],.025)           
  svsigma22ul[j]=quantile(1/svsigma22[,j],.975)           
}      


svsigma23m=matrix(NA,TM)                          
svsigma23up=matrix(NA,TM)                         
svsigma23ul=matrix(NA,TM)                         
for (j in 1:TM){                              
  svsigma23m[j]=quantile(1/svsigma23[,j],.5)              
  svsigma23up[j]=quantile(1/svsigma23[,j],.025)           
  svsigma23ul[j]=quantile(1/svsigma23[,j],.975)           
}                                             


svsigmaCUm=matrix(NA,TM)                          
svsigmaCUup=matrix(NA,TM)                         
svsigmaCUul=matrix(NA,TM)                         
for (j in 1:TM){                              
  svsigmaCUm[j]=quantile(1/svsigmaCU[,j],.5)              
  svsigmaCUup[j]=quantile(1/svsigmaCU[,j],.025)           
  svsigmaCUul[j]=quantile(1/svsigmaCU[,j],.975)           
}                                             


svsigmaPIm=matrix(NA,TM)                          
svsigmaPIup=matrix(NA,TM)                         
svsigmaPIul=matrix(NA,TM)                         
for (j in 1:TM){                              
  svsigmaPIm[j]=quantile(1/svsigmaPI[,j],.5)              
  svsigmaPIup[j]=quantile(1/svsigmaPI[,j],.025)           
  svsigmaPIul[j]=quantile(1/svsigmaPI[,j],.975)           
}         

svsigmaFFRm=matrix(NA,TM)                          
svsigmaFFRup=matrix(NA,TM)                         
svsigmaFFRul=matrix(NA,TM)                         
for (j in 1:TM){                              
  svsigmaFFRm[j]=quantile(1/svsigmaFFR[,j],.5)              
  svsigmaFFRup[j]=quantile(1/svsigmaFFR[,j],.025)           
  svsigmaFFRul[j]=quantile(1/svsigmaFFR[,j],.975)           
}                                             

svm=matrix(NA,TM,7)
svup=matrix(NA,TM,7)
svul=matrix(NA,TM,7)


svm[,1]=svsigma2m
svup[,1]=svsigma2up
svul[,1]=svsigma2ul

svm[,2]=svsigma21m
svup[,2]=svsigma21up
svul[,2]=svsigma21ul

svm[,3]=svsigma22m
svup[,3]=svsigma22up
svul[,3]=svsigma22ul

svm[,4]=svsigma23m
svup[,4]=svsigma23up
svul[,4]=svsigma23ul

svm[,5]=svsigmaCUm
svup[,5]=svsigmaCUup
svul[,5]=svsigmaCUul

svm[,6]=svsigmaPIm
svup[,6]=svsigmaPIup
svul[,6]=svsigmaPIul

svm[,7]=svsigmaFFRm
svup[,7]=svsigmaFFRup
svul[,7]=svsigmaFFRul


jpeg("svsigma2.jpeg")
plot(sqrt(svm[,1]),type="l",xlab="Tempo",ylab=expression(svsigma2[1]))
lines(sqrt(svup[,1]),lty=3)
lines(sqrt(svul[,1]),lty=3)
title(expression(svsigma2[1]))
dev.off()

jpeg("svsigma21.jpeg")
plot(sqrt(svm[,2]),type="l",xlab="Tempo",ylab=expression(svsigma2[2]))
lines(sqrt(svup[,2]),lty=3)
lines(sqrt(svul[,2]),lty=3)
title(expression(svsigma2[2]))
dev.off()

jpeg("svsigma22.jpeg")
plot(sqrt(svm[,3]),type="l",xlab="Tempo",ylab=expression(svsigma2[3]))
lines(sqrt(svup[,3]),lty=3)
lines(sqrt(svul[,3]),lty=3)
title(expression(svsigma2[3]))
dev.off()

jpeg("svsigma23.jpeg")
plot(sqrt(svm[,4]),type="l",xlab="Tempo",ylab=expression(svsigma2[4]))
lines(sqrt(svup[,4]),lty=3)
lines(sqrt(svul[,4]),lty=3)
title(expression(svsigma2[4]))
dev.off()

jpeg("svsigmaCU.jpeg")
plot(sqrt(svm[,5]),type="l",xlab="Tempo",ylab=expression(svsigma2[4]))
lines(sqrt(svup[,5]),lty=3)
lines(sqrt(svul[,5]),lty=3)
title(expression(svsigma2[4]))
dev.off()

jpeg("svsigmaPI.jpeg")
plot(sqrt(svm[,6]),type="l",xlab="Tempo",ylab=expression(svsigma2[4]))
lines(sqrt(svup[,6]),lty=3)
lines(sqrt(svul[,6]),lty=3)
title(expression(svsigma2[4]))
dev.off()

jpeg("svsigmaFFR.jpeg")
plot(sqrt(svm[,7]),type="l",xlab="Tempo",ylab=expression(svsigma2[4]))
lines(sqrt(svup[,7]),lty=3)
lines(sqrt(svul[,7]),lty=3)
title(expression(svsigma2[4]))
dev.off()

lambda1m=matrix(NA,TM)                          
lambda1up=matrix(NA,TM)                         
lambda1ul=matrix(NA,TM)                         
for (j in 1:TM){                              
  lambda1m[j]=quantile(lambda1[,j],.5)              
  lambda1up[j]=quantile(lambda1[,j],.025)           
  lambda1ul[j]=quantile(lambda1[,j],.975)           
}


jpeg("lambda1.jpeg")
plot(lambda1m,type="l",xlab="Tempo",ylab=expression(lambda1))
lines(lambda1up,lty=3)
lines(lambda1ul,lty=3)
title(expression(lambda1))
dev.off()

# Creates array with estimated values of the yield curve for every simulation on latent factors and lambda

jurosm2t=yields

nsims=floor((n.iter-n.burnin)/n.thin)
m=array(0,dim=c(nsims,dim(jurosm2t)[1],dim(jurosm2t)[2]))


for (sim in 1:nsims){
  for (kk in 1:TM){
    for (i in 1:counters[kk])
      m[sim,kk,i]<- beta1[sim,kk]+
        beta2[sim,kk]*(1-exp(-lambda1[sim,kk]*Maturity[kk,i]))/(lambda1[sim,kk]*Maturity[kk,i])+
        beta3[sim,kk]*((1-exp(-lambda1[sim,kk]*Maturity[kk,i]))/(lambda1[sim,kk]*Maturity[kk,i])-exp(-lambda1[sim,kk]*Maturity[kk,i]))
    
  }
}


# Matrix with mean and percentiles of estimated values for the yield curve

ypm=matrix(NA,TM,dim(jurosm2t)[2])
ypup=matrix(NA,TM,dim(jurosm2t)[2])
ypul=matrix(NA,TM,dim(jurosm2t)[2])
for (j in 1:TM){
  for (k in 1:counters[TM]){
    ypm[j,k]=quantile((m[,j,k]),.5,na.rm=T)
    ypup[j,k]=quantile((m[,j,k]),.025,na.rm=T)
    ypul[j,k]=quantile((m[,j,k]),.975,na.rm=T)
  }
}

for (kk in 1:TM){
  dayname=paste("anprev",kk,".jpeg",sep="")
  jpeg(dayname)
  plot(Maturity[kk,1:counters[kk]],yields[kk,1:counters[kk]],pch="*",ylab="Juros",xlab="Maturidade")
  lines(Maturity[kk,1:counters[kk]],ypm[kk,1:counters[kk]])
  lines(Maturity[kk,1:counters[kk]],ypup[kk,1:counters[kk]],col=2)
  lines(Maturity[kk,1:counters[kk]],ypul[kk,1:counters[kk]],col=2)
  dev.off()
}

apply(phi1,2,mean)
apply(phi1,2,quantile,.025)
apply(phi1,2,quantile,.975)

1/c(quantile(sv0,.025),mean(sv0),quantile(sv0,.975))
c(quantile(sv1,.025),mean(sv1),quantile(sv01,.975))
1/c(quantile(sv01,.025),mean(sv01),quantile(sv01,.975))
c(quantile(sv11,.025),mean(sv11),quantile(sv11,.975))
1/c(quantile(sv02,.025),mean(sv02),quantile(sv02,.975))
c(quantile(sv12,.025),mean(sv12),quantile(sv12,.975))
1/c(quantile(sv03,.025),mean(sv03),quantile(sv03,.975))
c(quantile(sv13,.025),mean(sv13),quantile(sv13,.975))
c(quantile(phil11,.025),mean(phil11),quantile(phil11,.975))
c(quantile(phil12,.025),mean(phil12),quantile(phil12,.975))


ts.plot(yields,col=seq(1,14,by=1))
legend(200,.18,c(vert),col=seq(1,14,by=1),lty=1)

library(forecast)

# Accuracy metrics for fitted yield curve values compared to original values

tf=length(ypm[,1])
bmodel1v1=accuracy(ypm[1:tf,1]*10000,yields[1:tf,1]*10000)
bmodel1v2=accuracy(ypm[1:tf,2]*10000,yields[1:tf,2]*10000)
bmodel1v3=accuracy(ypm[1:tf,3]*10000,yields[1:tf,3]*10000)
bmodel1v4=accuracy(ypm[1:tf,4]*10000,yields[1:tf,4]*10000)
bmodel1v5=accuracy(ypm[1:tf,5]*10000,yields[1:tf,5]*10000)
bmodel1v6=accuracy(ypm[1:tf,6]*10000,yields[1:tf,6]*10000)
bmodel1v7=accuracy(ypm[1:tf,7]*10000,yields[1:tf,7]*10000)
bmodel1v8=accuracy(ypm[1:tf,8]*10000,yields[1:tf,8]*10000)
bmodel1v9=accuracy(ypm[1:tf,9]*10000,yields[1:tf,9]*10000)
bmodel1v10=accuracy(ypm[1:tf,10]*10000,yields[1:tf,10]*10000)
bmodel1v11=accuracy(ypm[1:tf,11]*10000,yields[1:tf,11]*10000)
bmodel1v12=accuracy(ypm[1:tf,12]*10000,yields[1:tf,12]*10000)
bmodel1v13=accuracy(ypm[1:tf,13]*10000,yields[1:tf,13]*10000)
bmodel1v14=accuracy(ypm[1:tf,14]*10000,yields[1:tf,14]*10000)


bmodel1v1
bmodel1v2
bmodel1v3
bmodel1v4
bmodel1v5
bmodel1v6
bmodel1v7
bmodel1v8
bmodel1v9
bmodel1v10
bmodel1v11
bmodel1v12
bmodel1v13
bmodel1v14


comparacao2=matrix(0,14,7)
bcompmatrix2=as.data.frame(comparacao2)   
colnames(bcompmatrix2)=c("ME","SD","RMSE","MAE","MPE","MAPE","ACF1")
rownames(bcompmatrix2)=c("1", "2", "3", "4", "6", "9", "12", "15", "18", "24", "27", "29", "31", "33")

bcompmatrix2[1,c(1,3,4,5,6)]=	   bmodel1v1  
bcompmatrix2[2,c(1,3,4,5,6)]=	   bmodel1v2    
bcompmatrix2[3,c(1,3,4,5,6)]=	   bmodel1v3    
bcompmatrix2[4,c(1,3,4,5,6)]=	   bmodel1v4    
bcompmatrix2[5,c(1,3,4,5,6)]=	   bmodel1v5    
bcompmatrix2[6,c(1,3,4,5,6)]=	   bmodel1v6    
bcompmatrix2[7,c(1,3,4,5,6)]=	   bmodel1v7    
bcompmatrix2[8,c(1,3,4,5,6)]=	   bmodel1v8    
bcompmatrix2[9,c(1,3,4,5,6)]=	   bmodel1v9    
bcompmatrix2[10,c(1,3,4,5,6)]=   bmodel1v10    
bcompmatrix2[11,c(1,3,4,5,6)]=   bmodel1v11    
bcompmatrix2[12,c(1,3,4,5,6)]=   bmodel1v12    
bcompmatrix2[13,c(1,3,4,5,6)]=   bmodel1v13    
bcompmatrix2[14,c(1,3,4,5,6)]=   bmodel1v14


bcompmatrix2[1,2]=  sd(ypm[1:tf,1]*10000-yields[1:tf,1]*10000)   	 
bcompmatrix2[2,2]=  sd(ypm[1:tf,2]*10000-yields[1:tf,2]*10000)     
bcompmatrix2[3,2]=  sd(ypm[1:tf,3]*10000-yields[1:tf,3]*10000)     
bcompmatrix2[4,2]=  sd(ypm[1:tf,4]*10000-yields[1:tf,4]*10000)     
bcompmatrix2[5,2]=  sd(ypm[1:tf,5]*10000-yields[1:tf,5]*10000)     
bcompmatrix2[6,2]=  sd(ypm[1:tf,6]*10000-yields[1:tf,6]*10000)     
bcompmatrix2[7,2]=  sd(ypm[1:tf,7]*10000-yields[1:tf,7]*10000)     
bcompmatrix2[8,2]=  sd(ypm[1:tf,8]*10000-yields[1:tf,8]*10000)     
bcompmatrix2[9,2]=  sd(ypm[1:tf,9]*10000-yields[1:tf,9]*10000)     
bcompmatrix2[10,2]= sd(ypm[1:tf,10]*10000-yields[1:tf,10]*10000)
bcompmatrix2[11,2]= sd(ypm[1:tf,11]*10000-yields[1:tf,11]*10000)
bcompmatrix2[12,2]= sd(ypm[1:tf,12]*10000-yields[1:tf,12]*10000)
bcompmatrix2[13,2]= sd(ypm[1:tf,13]*10000-yields[1:tf,13]*10000)
bcompmatrix2[14,2]= sd(ypm[1:tf,14]*10000-yields[1:tf,14]*10000)


bcompmatrix2[1,7]=  acf(ypm[1:tf,1]*10000-yields[1:tf,1]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2]    	 
bcompmatrix2[2,7]=  acf(ypm[1:tf,2]*10000-yields[1:tf,2]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2]      
bcompmatrix2[3,7]=  acf(ypm[1:tf,3]*10000-yields[1:tf,3]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
bcompmatrix2[4,7]=  acf(ypm[1:tf,4]*10000-yields[1:tf,4]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2]      
bcompmatrix2[5,7]=  acf(ypm[1:tf,5]*10000-yields[1:tf,5]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
bcompmatrix2[6,7]=  acf(ypm[1:tf,6]*10000-yields[1:tf,6]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
bcompmatrix2[7,7]=  acf(ypm[1:tf,7]*10000-yields[1:tf,7]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
bcompmatrix2[8,7]=  acf(ypm[1:tf,8]*10000-yields[1:tf,8]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
bcompmatrix2[9,7]=  acf(ypm[1:tf,9]*10000-yields[1:tf,9]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
bcompmatrix2[10,7]= acf(ypm[1:tf,10]*10000-yields[1:tf,10]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
bcompmatrix2[11,7]= acf(ypm[1:tf,11]*10000-yields[1:tf,11]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
bcompmatrix2[12,7]= acf(ypm[1:tf,12]*10000-yields[1:tf,12]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
bcompmatrix2[13,7]= acf(ypm[1:tf,13]*10000-yields[1:tf,13]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
bcompmatrix2[14,7]= acf(ypm[1:tf,14]*10000-yields[1:tf,14]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 



tf=length(ypm[,1])
rbmodel1v1=accuracy(ypm[5:tf,1]*10000,yields[5:tf,1]*10000)
rbmodel1v2=accuracy(ypm[5:tf,2]*10000,yields[5:tf,2]*10000)
rbmodel1v3=accuracy(ypm[5:tf,3]*10000,yields[5:tf,3]*10000)
rbmodel1v4=accuracy(ypm[5:tf,4]*10000,yields[5:tf,4]*10000)
rbmodel1v5=accuracy(ypm[5:tf,5]*10000,yields[5:tf,5]*10000)
rbmodel1v6=accuracy(ypm[5:tf,6]*10000,yields[5:tf,6]*10000)
rbmodel1v7=accuracy(ypm[5:tf,7]*10000,yields[5:tf,7]*10000)
rbmodel1v8=accuracy(ypm[5:tf,8]*10000,yields[5:tf,8]*10000)
rbmodel1v9=accuracy(ypm[5:tf,9]*10000,yields[5:tf,9]*10000)
rbmodel1v10=accuracy(ypm[5:tf,10]*10000,yields[5:tf,10]*10000)
rbmodel1v11=accuracy(ypm[5:tf,11]*10000,yields[5:tf,11]*10000)
rbmodel1v12=accuracy(ypm[5:tf,12]*10000,yields[5:tf,12]*10000)
rbmodel1v13=accuracy(ypm[5:tf,13]*10000,yields[5:tf,13]*10000)
rbmodel1v14=accuracy(ypm[5:tf,14]*10000,yields[5:tf,14]*10000)


rbmodel1v1
rbmodel1v2
rbmodel1v3
rbmodel1v4
rbmodel1v5
rbmodel1v6
rbmodel1v7
rbmodel1v8
rbmodel1v9
rbmodel1v10
rbmodel1v11
rbmodel1v12
rbmodel1v13
rbmodel1v14

comparacao3=matrix(0,14,7)
rbcompmatrix2=as.data.frame(comparacao2)   
colnames(rbcompmatrix2)=c("ME","SD","RMSE","MAE","MPE","MAPE","ACF1")
rownames(rbcompmatrix2)=c("1", "2", "3", "4", "6", "9", "12", "15", "18", "24", "27", "29", "31", "33")

rbcompmatrix2[1,c(1,3,4,5,6)]=	   rbmodel1v1  
rbcompmatrix2[2,c(1,3,4,5,6)]=	   rbmodel1v2    
rbcompmatrix2[3,c(1,3,4,5,6)]=	   rbmodel1v3    
rbcompmatrix2[4,c(1,3,4,5,6)]=	   rbmodel1v4    
rbcompmatrix2[5,c(1,3,4,5,6)]=	   rbmodel1v5    
rbcompmatrix2[6,c(1,3,4,5,6)]=	   rbmodel1v6    
rbcompmatrix2[7,c(1,3,4,5,6)]=	   rbmodel1v7    
rbcompmatrix2[8,c(1,3,4,5,6)]=	   rbmodel1v8    
rbcompmatrix2[9,c(1,3,4,5,6)]=	   rbmodel1v9    
rbcompmatrix2[10,c(1,3,4,5,6)]=   rbmodel1v10    
rbcompmatrix2[11,c(1,3,4,5,6)]=   rbmodel1v11    
rbcompmatrix2[12,c(1,3,4,5,6)]=   rbmodel1v12    
rbcompmatrix2[13,c(1,3,4,5,6)]=   rbmodel1v13    
rbcompmatrix2[14,c(1,3,4,5,6)]=   rbmodel1v14


rbcompmatrix2[1,2]=  sd(ypm[5:tf,1]*10000-yields[5:tf,1]*10000)   	 
rbcompmatrix2[2,2]=  sd(ypm[5:tf,2]*10000-yields[5:tf,2]*10000)     
rbcompmatrix2[3,2]=  sd(ypm[5:tf,3]*10000-yields[5:tf,3]*10000)     
rbcompmatrix2[4,2]=  sd(ypm[5:tf,4]*10000-yields[5:tf,4]*10000)     
rbcompmatrix2[5,2]=  sd(ypm[5:tf,5]*10000-yields[5:tf,5]*10000)     
rbcompmatrix2[6,2]=  sd(ypm[5:tf,6]*10000-yields[5:tf,6]*10000)     
rbcompmatrix2[7,2]=  sd(ypm[5:tf,7]*10000-yields[5:tf,7]*10000)     
rbcompmatrix2[8,2]=  sd(ypm[5:tf,8]*10000-yields[5:tf,8]*10000)     
rbcompmatrix2[9,2]=  sd(ypm[5:tf,9]*10000-yields[5:tf,9]*10000)     
rbcompmatrix2[10,2]= sd(ypm[5:tf,10]*10000-yields[5:tf,10]*10000)
rbcompmatrix2[11,2]= sd(ypm[5:tf,11]*10000-yields[5:tf,11]*10000)
rbcompmatrix2[12,2]= sd(ypm[5:tf,12]*10000-yields[5:tf,12]*10000)
rbcompmatrix2[13,2]= sd(ypm[5:tf,13]*10000-yields[5:tf,13]*10000)
rbcompmatrix2[14,2]= sd(ypm[5:tf,14]*10000-yields[5:tf,14]*10000)



rbcompmatrix2[1,7]=  acf(ypm[5:tf,1]*10000-yields[5:tf,1]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2]    	 
rbcompmatrix2[2,7]=  acf(ypm[5:tf,2]*10000-yields[5:tf,2]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2]      
rbcompmatrix2[3,7]=  acf(ypm[5:tf,3]*10000-yields[5:tf,3]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
rbcompmatrix2[4,7]=  acf(ypm[5:tf,4]*10000-yields[5:tf,4]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2]      
rbcompmatrix2[5,7]=  acf(ypm[5:tf,5]*10000-yields[5:tf,5]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
rbcompmatrix2[6,7]=  acf(ypm[5:tf,6]*10000-yields[5:tf,6]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
rbcompmatrix2[7,7]=  acf(ypm[5:tf,7]*10000-yields[5:tf,7]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
rbcompmatrix2[8,7]=  acf(ypm[5:tf,8]*10000-yields[5:tf,8]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
rbcompmatrix2[9,7]=  acf(ypm[5:tf,9]*10000-yields[5:tf,9]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
rbcompmatrix2[10,7]= acf(ypm[5:tf,10]*10000-yields[5:tf,10]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
rbcompmatrix2[11,7]= acf(ypm[5:tf,11]*10000-yields[5:tf,11]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
rbcompmatrix2[12,7]= acf(ypm[5:tf,12]*10000-yields[5:tf,12]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
rbcompmatrix2[13,7]= acf(ypm[5:tf,13]*10000-yields[5:tf,13]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
rbcompmatrix2[14,7]= acf(ypm[5:tf,14]*10000-yields[5:tf,14]*10000,type="correlation",lag.max = 12,plot=FALSE)$acf[2] 
