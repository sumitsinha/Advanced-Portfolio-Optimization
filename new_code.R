#EGARCH MODELLING
egarchsnp.spec = ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model=list(armaOrder=c(0,0)),distribution.model = "sstd")
apple.egarch = ugarchfit(egarchsnp.spec,apple.ts)
ibm.egarch=ugarchfit(egarchsnp.spec,ibm.ts)
lenovo.egarch=ugarchfit(egarchsnp.spec,invgy.ts)
yandex.egarch=ugarchfit(egarchsnp.spec,yandex.ts)
microsoft.egarch=ugarchfit(egarchsnp.spec,microsoft.ts)
bidu.egarch=ugarchfit(egarchsnp.spec,bidu.ts)

#MODEL DESCRIPTION
apple.egarch
ibm.egarch
lenovo.egarch
yandex.egarch
microsoft.egarch
bidu.egarch

#ECDF FUNCTIONS
apple.egarchecdf<-empDist(apple.egarch@fit$residuals/apple.egarch@fit$sigma)
ibm.egarchecdf<-empDist(ibm.egarch@fit$residuals/ibm.egarch@fit$sigma)
lenovo.egarchecdf<-empDist(lenovo.egarch@fit$residuals/lenovo.egarch@fit$sigma)
yandex.egarchecdf<-empDist(yandex.egarch@fit$residuals/yandex.egarch@fit$sigma)
microsoft.egarchecdf<-empDist(microsoft.egarch@fit$residuals/microsoft.egarch@fit$sigma)
bidu.egarchecdf<-empDist(bidu.egarch@fit$residuals/bidu.egarch@fit$sigma)

#COEFFecients
coef(apple.egarch)
coef(apple.garch)

#Tau Matrix
unifVaRegarch<-data.frame(cbind(apple.egarchecdf,ibm.egarchecdf,lenovo.egarchecdf,yandex.egarchecdf,microsoft.egarchecdf,bidu.egarchecdf))
interdependencyMatrixEgarcg<-TauMatrix(unifVaRegarch)
View(interdependencyMatrixEgarcg)

#CVM Structure
cvmegarcg <- RVineStructureSelect(unifVaRegarch, type=0, c(1:4))
set.seed(2) 
uSimEGARCH<-RVineSim(5000,cvmegarcg)
summary(cvmegarcg)
cvmegarcg

#ACF PLots
par(mfrow = c(3, 2))
acf(abs(abs(apple.egarch@fit$residuals/apple.egarch@fit$sigma)),main="APPLE")
acf(abs(abs(ibm.egarch@fit$residuals/ibm.egarch@fit$sigma)),main="IBM")
acf(abs(abs(lenovo.egarch@fit$residuals/lenovo.egarch@fit$sigma)),main="LENOVO")
acf(abs(abs(microsoft.egarch@fit$residuals/microsoft.egarch@fit$sigma)),main="MICROSOFT")
acf(abs(abs(yandex.egarch@fit$residuals/yandex.egarch@fit$sigma)),main="YANDEX")
acf(abs(abs(bidu.egarch@fit$residuals/bidu.egarch@fit$sigma)),main="BIDU")
par(mfrow = c(1, 1))

#Simulation from copula and fitting error dsitribution
zSim.appleEGARCH<-qstd(uSimEGARCH[,1],sd=matrix(coef(apple.egarch))[6,], nu=matrix(coef(apple.egarch ))[7,])/sd(qstd(uSimEGARCH[,1],sd=matrix(coef(apple.egarch))[6,], nu=matrix( coef(apple.egarch))[7,]))
zSim.ibmEGARCH<-qstd(uSimEGARCH[,2],sd=matrix(coef(ibm.egarch))[6,], nu=matrix(coef(ibm.egarch ))[7,])/sd(qstd(uSimEGARCH[,2],sd=matrix(coef(ibm.egarch))[6,], nu=matrix(coef(ibm.egarch))[7,]))
zSim.lenovoEGARCH<-qstd(uSimEGARCH[,3],sd=matrix(coef(lenovo.egarch))[6,], nu=matrix(coef(lenovo.egarch ))[7,])/sd(qstd(uSimEGARCH[,3],sd=matrix(coef(lenovo.egarch))[6,], nu=matrix(coef(lenovo.egarch))[7,]))
zSim.yandexEGARCH<-qstd(uSimEGARCH[,4],sd=matrix(coef(yandex.egarch))[6,], nu=matrix(coef(yandex.egarch ))[7,])/sd(qstd(uSimEGARCH[,4],sd=matrix(coef(yandex.egarch))[6,], nu=matrix(coef(yandex.egarch))[7,]))
zSim.microsoftEGARCH<-qstd(uSimEGARCH[,5],sd=matrix(coef(microsoft.egarch))[6,], nu=matrix(coef(microsoft.egarch ))[7,])/sd(qstd(uSimEGARCH[,5],sd=matrix(coef(microsoft.egarch))[6,], nu=matrix(coef(microsoft.egarch))[7,]))
zSim.biduEGARCH<-qstd(uSimEGARCH[,6],sd=matrix(coef(bidu.egarch))[6,], nu=matrix(coef(bidu.egarch ))[7,])/sd(qstd(uSimEGARCH[,6],sd=matrix(coef(bidu.egarch))[6,], nu=matrix(coef(bidu.egarch))[7,]))

#APPLE ESTIMATES
xEGARCH <-cbind(coef(apple.egarch))
alpha0.appleEGARCH<-xEGARCH[2,1]
alpha1.appleEGARH<-xEGARCH[3,1]
beta1.appleEGARH<-xEGARCH[4,1]
gama1.appleEGARCH<-xEGARCH[5,1]
z.appleEGARCH<-zSim.appleEGARCH
a.appleEGARCH<-array(0,5000)
s.appleEGARCH<-array(0,5000)
s.appleEGARCH.start<-matrix(apple.egarch@fit$sigma)[755,]
a.appleEGARCH.start<- matrix(apple.egarch@fit$residuals/apple.egarch@fit$sigma)[755,]
for(i in 1:5000)
{
  s.appleEGARCH[i]<-sqrt(exp(alpha0.appleEGARCH+alpha1.appleEGARH*abs(a.appleEGARCH.start)+gama1.appleEGARCH*a.appleEGARCH.start + beta1.appleEGARH*log(s.appleEGARCH.start^2)))
  a.appleEGARCH[i]<-z.appleEGARCH[i]*s.appleEGARCH[i]
}
rSim.appleEGARCH<-mean(apple.ts)+a.appleEGARCH
rSim.appleEGARCH1<-apple.ts[755]+a.appleEGARCH


#IBM ESTIMATES
yEGARCH <-cbind(coef(ibm.egarch))
alpha0.ibmEGARCH<-yEGARCH[2,1]
alpha1.ibmEGARH<-yEGARCH[3,1]
beta1.ibmEGARH<-yEGARCH[4,1]
gama1.ibmEGARCH<-yEGARCH[5,1]
z.ibmEGARCH<-zSim.ibmEGARCH
a.ibmEGARCH<-array(0,5000)
s.ibmEGARCH<-array(0,5000)
s.ibmEGARCH.start<-matrix(ibm.egarch@fit$sigma)[755,]
a.ibmEGARCH.start<- matrix(ibm.egarch@fit$residuals/ibm.egarch@fit$sigma)[755,]
for(i in 1:5000)
{
  s.ibmEGARCH[i]<-sqrt(exp(alpha0.ibmEGARCH+alpha1.ibmEGARH*abs(a.ibmEGARCH.start)+gama1.ibmEGARCH*a.ibmEGARCH.start + beta1.ibmEGARH*log(s.ibmEGARCH.start^2)))
  a.ibmEGARCH[i]<-z.ibmEGARCH[i]*s.ibmEGARCH[i]
}
rSim.ibmEGARCH<-mean(ibm.ts)+a.ibmEGARCH
rSim.ibmEGARCH1<-ibm.ts[755]+a.ibmEGARCH

#LENEVO ESTIMATES
zEGARCH <-cbind(coef(lenovo.egarch))
alpha0.lenovoEGARCH<-zEGARCH[2,1]
alpha1.lenovoEGARH<-zEGARCH[3,1]
beta1.lenovoEGARH<-zEGARCH[4,1]
gama1.lenovoEGARCH<-zEGARCH[5,1]
z.lenovoEGARCH<-zSim.lenovoEGARCH
a.lenovoEGARCH<-array(0,5000)
s.lenovoEGARCH<-array(0,5000)
s.lenovoEGARCH.start<-matrix(lenovo.egarch@fit$sigma)[755,]
a.lenovoEGARCH.start<- matrix(lenovo.egarch@fit$residuals/lenovo.egarch@fit$sigma)[755,]
for(i in 1:5000)
{
  s.lenovoEGARCH[i]<-sqrt(exp(alpha0.lenovoEGARCH+alpha1.lenovoEGARH*abs(a.lenovoEGARCH.start)+gama1.lenovoEGARCH*a.lenovoEGARCH.start + beta1.lenovoEGARH*log(s.lenovoEGARCH.start^2)))
  a.lenovoEGARCH[i]<-z.lenovoEGARCH[i]*s.lenovoEGARCH[i]
}
rSim.lenovoEGARCH<-mean(invgy.ts)+a.lenovoEGARCH
rSim.lenovoEGARCH1<-invgy.ts[755]+a.lenovoEGARCH

#YANDEX ESTIMATES
mEGARCH <-cbind(coef(yandex.egarch))
alpha0.yandexEGARCH<-mEGARCH[2,1]
alpha1.yandexEGARH<-mEGARCH[3,1]
beta1.yandexEGARH<-mEGARCH[4,1]
gama1.yandexEGARCH<-mEGARCH[5,1]
z.yandexEGARCH<-zSim.yandexEGARCH
a.yandexEGARCH<-array(0,5000)
s.yandexEGARCH<-array(0,5000)
s.yandexEGARCH.start<-matrix(yandex.egarch@fit$sigma)[755,]
a.yandexEGARCH.start<- matrix(yandex.egarch@fit$residuals/yandex.egarch@fit$sigma)[755,]
for(i in 1:5000)
{
  s.yandexEGARCH[i]<-sqrt(exp(alpha0.yandexEGARCH+alpha1.yandexEGARH*abs(a.yandexEGARCH.start)+gama1.yandexEGARCH*a.yandexEGARCH.start + beta1.yandexEGARH*log(s.yandexEGARCH.start^2)))
  a.yandexEGARCH[i]<-z.yandexEGARCH[i]*s.yandexEGARCH[i]
}
rSim.yandexEGARCH<-mean(yandex.ts)+a.yandexEGARCH
rSim.yandexEGARCH1<-yandex.ts[755]+a.yandexEGARCH

# MICROSOFT ESTIMATES
qEGARCH <-cbind(coef(microsoft.egarch))
alpha0.microsoftEGARCH<-qEGARCH[2,1]
alpha1.microsoftEGARH<-qEGARCH[3,1]
beta1.microsoftEGARH<-qEGARCH[4,1]
gama1.microsoftEGARCH<-qEGARCH[5,1]
z.microsoftEGARCH<-zSim.microsoftEGARCH
a.microsoftEGARCH<-array(0,5000)
s.microsoftEGARCH<-array(0,5000)
s.microsoftEGARCH.start<-matrix(microsoft.egarch@fit$sigma)[755,]
a.microsoftEGARCH.start<- matrix(microsoft.egarch@fit$residuals/microsoft.egarch@fit$sigma)[755,]
for(i in 1:5000)
{
  s.microsoftEGARCH[i]<-sqrt(exp(alpha0.microsoftEGARCH+alpha1.microsoftEGARH*abs(a.microsoftEGARCH.start)+gama1.microsoftEGARCH*a.microsoftEGARCH.start + beta1.microsoftEGARH*log(s.microsoftEGARCH.start^2)))
  a.microsoftEGARCH[i]<-z.microsoftEGARCH[i]*s.microsoftEGARCH[i]
}
rSim.microsoftEGARCH<-mean(microsoft.ts)+a.microsoftEGARCH
rSim.microsoftEGARCH1<-microsoft.ts[755]+a.microsoftEGARCH

#BIDU ESTIMATES
bEGARCH <-cbind(coef(bidu.egarch))
alpha0.biduEGARCH<-bEGARCH[2,1]
alpha1.biduEGARH<-bEGARCH[3,1]
beta1.biduEGARH<-bEGARCH[4,1]
gama1.biduEGARCH<-bEGARCH[5,1]
z.biduEGARCH<-zSim.biduEGARCH
a.biduEGARCH<-array(0,5000)
s.biduEGARCH<-array(0,5000)
s.biduEGARCH.start<-matrix(bidu.egarch@fit$sigma)[755,]
a.biduEGARCH.start<- matrix(bidu.egarch@fit$residuals/bidu.egarch@fit$sigma)[755,]
for(i in 1:5000)
{
  s.biduEGARCH[i]<-sqrt(exp(alpha0.biduEGARCH+alpha1.biduEGARH*abs(a.biduEGARCH.start)+gama1.biduEGARCH*a.biduEGARCH.start + beta1.biduEGARH*log(s.biduEGARCH.start^2)))
  a.biduEGARCH[i]<-z.biduEGARCH[i]*s.biduEGARCH[i]
}
rSim.biduEGARCH<-mean(bidu.ts)+a.biduEGARCH
rSim.biduEGARCH1<-bidu.ts[755]+a.biduEGARCH

#Short term mean calulation
apple.total100=0
ibm.total100=0
lenovo.total100=0
yandex.total100=0
microsoft.total100=0
bidu.total100=0
for(j in 1:100)
{
  apple.total100=apple.total100+apple.ts[655+j]
  ibm.total100=ibm.total100+ibm.ts[655+j]
  lenovo.total100=lenovo.total100+invgy.ts[655+j]
  yandex.total100=yandex.total100+yandex.ts[655+j]
  microsoft.total100=microsoft.total100+microsoft.ts[655+j]
  bidu.total100=bidu.total100+bidu.ts[655+j]
}
apple.mean100=apple.total100/100
ibm.mean100=ibm.total100/100
lenovo.mean100=lenovo.total100/100
yandex.mean100=yandex.total100/100
microsoft.mean100=microsoft.total100/100
bidu.mean100=bidu.total100/100

rSim.appleEGARCH2<-apple.mean100+a.appleEGARCH
rSim.ibmEGARCH2<-ibm.mean100+a.ibmEGARCH
rSim.lenovoEGARCH2<-lenovo.mean100+a.lenovoEGARCH
rSim.yandexEGARCH2<-yandex.mean100+a.yandexEGARCH
rSim.microsoftEGARCH2<-microsoft.mean100+a.microsoftEGARCH
rSim.biduEGARCH2<-bidu.mean100+a.biduEGARCH

# Comuputing final ARRAY
dataEGARCH<-cbind(rSim.appleEGARCH,rSim.ibmEGARCH,rSim.lenovoEGARCH,rSim.yandexEGARCH,rSim.microsoftEGARCH,rSim.biduEGARCH)
tsDataEGARCH <- timeSeries(dataEGARCH)
dataEGARCH1<-cbind(rSim.appleEGARCH1,rSim.ibmEGARCH1,rSim.lenovoEGARCH1,rSim.yandexEGARCH1,rSim.microsoftEGARCH1,rSim.biduEGARCH1)
tsDataEGARCH1 <- timeSeries(dataEGARCH1)
dataEGARCH2<-cbind(rSim.appleEGARCH2,rSim.ibmEGARCH2,rSim.lenovoEGARCH2,rSim.yandexEGARCH2,rSim.microsoftEGARCH2,rSim.biduEGARCH2)
tsDataEGARCH2 <- timeSeries(dataEGARCH2)

#CVAR OPTIMIZAITON USING EGRACH-PCC CALCULATIONS LONG TERM
frontierSpec<-portfolioSpec()
setType(frontierSpec)<-"CVaR" 
setSolver(frontierSpec)<-"solveRglpk.CVAR"

setAlpha(frontierSpec)<-0.01 
setNFrontierPoints(frontierSpec)<-50

frontier1gEGARCH<-portfolioFrontier(data=tsDataEGARCH, spec=frontierSpec,constraints="LongOnly")
summary(frontier1gEGARCH)
VaRtab1gEGARCH <- getTargetRisk(frontier1gEGARCH@portfolio) 
muTab1gEGARCH<-getTargetReturn(frontier1gEGARCH@portfolio) 
weightTab1gEGARCH<- getWeights(frontier1gEGARCH@portfolio)

indE<-1
while(muTab1gEGARCH[indE,1]<0.06/250) 
  indE<-indE+1
optWeights1gEGARCH <- weightTab1gEGARCH[indE,] 
optWeights1gEGARCH
VaRtab1gEGARCH[indE,]
muTab1gEGARCH[indE,]

View(weightTab1gEGARCH)
View(VaRtab1gEGARCH)
View(muTab1gEGARCH)

tailoredFrontierPlot(object=frontier1gEGARCH,mText="EGARCH PCC 99%-CVaR Portfolio2(Long only constraints)",risk="CVaR")

#CVAR OPTIMIZAITON USING EGRACH-PCC CALCULATION SHORT TERM
frontierSpec<-portfolioSpec()
setType(frontierSpec)<-"CVaR" 
setSolver(frontierSpec)<-"solveRglpk.CVAR"

setAlpha(frontierSpec)<-0.01 
setNFrontierPoints(frontierSpec)<-50

frontier1gEGARCH1<-portfolioFrontier(data=tsDataEGARCH1, spec=frontierSpec,constraints="LongOnly")
summary(frontier1gEGARCH1)
VaRtab1gEGARCH1 <- getTargetRisk(frontier1gEGARCH1@portfolio) 
muTab1gEGARCH1<-getTargetReturn(frontier1gEGARCH1@portfolio) 
weightTab1gEGARCH1<- getWeights(frontier1gEGARCH1@portfolio)

indE1<-1
while(muTab1gEGARCH1[indE1,1]<0.06/250) 
  indE1<-indE1+1
optWeights1gEGARCH1 <- weightTab1gEGARCH1[indE1,] 
optWeights1gEGARCH1
VaRtab1gEGARCH1[indE1,]
muTab1gEGARCH1[indE1,]

tailoredFrontierPlot(object=frontier1gEGARCH1,mText="EGARCH SHORT TERM PCC 99%-CVaR Portfolio2(Long only constraints)",risk="CVaR")

#CVAR OPTIMIZAITON USING EGRACH-PCC CALCULATION 100 DAYS RECENT SHORT TERM
frontierSpec<-portfolioSpec()
setType(frontierSpec)<-"CVaR" 
setSolver(frontierSpec)<-"solveRglpk.CVAR"

setAlpha(frontierSpec)<-0.01 
setNFrontierPoints(frontierSpec)<-50

frontier1gEGARCH2<-portfolioFrontier(data=tsDataEGARCH2, spec=frontierSpec,constraints="LongOnly")
summary(frontier1gEGARCH2)
VaRtab1gEGARCH2 <- getTargetRisk(frontier1gEGARCH2@portfolio) 
muTab1gEGARCH2<-getTargetReturn(frontier1gEGARCH2@portfolio) 
weightTab1gEGARCH2<- getWeights(frontier1gEGARCH2@portfolio)

indE2<-1
while(muTab1gEGARCH2[indE2,1]<0.06/250) 
  indE2<-indE2+1
optWeights1gEGARCH2 <- weightTab1gEGARCH2[indE2,] 
optWeights1gEGARCH2
VaRtab1gEGARCH2[indE2,]
muTab1gEGARCH2[indE2,]
View(muTab1gEGARCH2)
tailoredFrontierPlot(object=frontier1gEGARCH2,mText="EGARCH 100 DAYS SHORT TERM PCC 99%-CVaR Portfolio2(Long only constraints)",risk="CVaR")
