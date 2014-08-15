rm(list=ls(all=TRUE))

devtools::install_github("dapr12/fda.classification")
library('fda.classification')

## Growth Data Curves

par(mfrow=c(1,2)) 
matplot(Gwd$Age, Gwd$Male, type="b", pch = 1, main="Male", xlab="Age (years)", ylab="Heigth (cm)")
grid(NA, 20)
matplot(Gwd$Age, Gwd$Female, type="b", pch = 1,main="Female", xlab="Age (years)", ylab="Height (cm)")
grid(NA, 20)

## Protein Data Curves

par(mfrow=c(1,2), font.lab=3, font.axis=1) 
matplot(Gwdprotein$Time, Gwdprotein$Group1, type="l", main="Wheat Data Set", sub="Group 1", xlab="Wave length(nm)", ylab=" d(NIR spectra ( log(1/R) ),1)")
abline(h= 0)

matplot(Gwdprotein$Time, Gwdprotein$Group2, type="l", main="Wheat Data Set", sub="Group 2", xlab="Wave length(nm)", ylab=" d(NIR spectra ( log(1/R) ),1) ")
abline(h= 0)


## Nox Data Curves

par(mfrow=c(1,2), font.lab=3, font.axis=1) 
matplot(Nox$TimeHours, Nox$Working, type="l", main="NOx Levels", sub="Group Working Days", xlab="Time(Hours)", ylab="mgm^3")

matplot(Nox$TimeHours, Nox$NonWorking, type="l",main="NOx Levels", sub="Group Non Working Days", xlab="Time(Hours)", ylab=" mgm^3")

## Plasma Citrate Data Curves

matplot(Plasma$TimeDay, Plasma$Plasma, type="l", main="Plasma Citrate Data", xlab="Time(Hours)", ylab="Concentration")

## Phoneme Data Curves

par(mfrow=c(3,2), font.lab=3, font.axis=1) 
matplot(Phoneme$Time, Phoneme$aa, type="l", main="Phoneme", sub="aa", xlab="Time", ylab="log-periodogram")
matplot(Phoneme$Time, Phoneme$ao, type="l", main="Phoneme", sub="ao", xlab="Time", ylab="log-periodogram")
matplot(Phoneme$Time, Phoneme$dcl, type="l", main="Phoneme", sub="dcl", xlab="Time", ylab="log-periodogram")
matplot(Phoneme$Time, Phoneme$iy, type="l", main="Phoneme", sub="iy", xlab="Time", ylab="log-periodogram")
matplot(Phoneme$Time, Phoneme$sh, type="l", main="Phoneme", sub="sh", xlab="Time", ylab="log-periodogram")


##############################
# Creating FDA Class - Gwt Males 
##############################

fdaobjMale<-fdaclass(Gwd$Male,Gwd$Age,c(1,18))
plot(fdaobjMale)

##############################
# Outliergram 
##############################

OG<-outliergram(fdaobjMale, remove=TRUE)
OG$outliers

#############################################
# FDA Class and Outlies Plasma Citrate Data
#############################################

fdaNox<-fdaclass(Nox$NonWorking,Nox$TimeHours )
plot(fdaNox)
OutliersNox<-outliergram(fdaNox, remove=TRUE)
OutliersNox$outliers

NoxNoOutliers<- fdaclass(OutliersNox$data, OutliersNox$argvals)
plot(NoxNoOutliers)

##############################
# Smooth and CV - Male Growth Data 
##############################

### Smooth and CV
Smoothfda<- smoothfda(fdaobjMale, bandwidth=0.500019, degree=1)
matplot(Smoothfda$YSmooth, type="l", xlab="Time",ylab="x(t)",main="Male Growth Data - Mean and Median")
lines(Smoothfda$Mean,col="black",lwd = 4) 
Mfd<-medianfd(Smoothfda)
lines(Mfd$median,col="red",lwd = 4) 


##############################
# h Select and Smooth and CV 
##############################

Int<-c(1,2)
h<- hselect(NoxNoOutliers, 1, Int)
print(h)

SmoothNox<- smoothfda(fdaNox, bandwidth= h, degree=1)
matplot(SmoothNox$YSmooth, type="l",ylab="x(t)",main="Observations")
lines(SmoothNox$Mean,col="black",lwd = 4) 


##############################
# First Derivative 
##############################

FderivNox<-fderiv(fdaNox,1)
matplot(FderivNox$data, type="l")


#######################################
# Median and Mean - Males Growth Data
#######################################

#Computing the Median
Mfd<-medianfd(Smoothfda)

# Plotting the Median
par(mfrow=c(1,2))
matplot(Smoothfda$YSmooth, type="l",main="Mean - Male Growth Data")
lines(Smoothfda$Mean,col="black",lwd = 4) 
matplot(Smoothfda$YSmooth, type="l", main="Median- Male Growth Data")
lines(Mfd$median, col="red",lwd = 4) 


#######################################
# Median and Mean - Males Growth Data
#######################################

SmoothNox<- smoothfda(fdaNox, bandwidth= h, degree=1)
matplot(SmoothNox$YSmooth, type="l",ylab="x(t)",main="Observations")
lines(SmoothNox$Mean,col="black",lwd = 4) 

#Computing the Median
MfdNox<-medianfd(SmoothNox)

# Plotting the Median
par(mfrow=c(1,2))
matplot(SmoothNox$YSmooth, type="l",main="Mean")
lines(SmoothNox$Mean,col="black",lwd = 4) 
matplot(SmoothNox$YSmooth, type="l", main="Median")
lines(MfdNox$median, col="red",lwd = 4) 


#########################################
# Smooth with BSplines - Males Gwt Data
#########################################

norder<-6
lambda <- 0.01
fdaobj<-fdaobjMale
Lf<-4
Intv<-c(-6,0)
BSplinesGrowth <-smoothbsplines( fdaobj, norder, lambda, Lf=4 )

par(mfrow=c(1,2))
matplot(BSplinesGrowth$fdSmooth$coefs, type="l",main="Smooth Data using B-Splines",xlab="Age (years)", ylab="Heigth (cm)")
GCV <- gcsvc( fdaobj, norder, lambda, Lf, Intv  )



##############################
# PCA
##############################

## D-Plot for Males GWt Data 
Var<-varfd(Smoothfda$YSmooth)
Eg<-eigen(Var)
Dplt<- dpout(Eg$val, plotting=TRUE)


## D-Plot for Nox Data 
Var<-varfd(SmoothNox$YSmooth)
Eg<-eigen(Var)
Dplt<- dpout(Eg$val, plotting=TRUE)


## PCA - PCA Functional Data Analysis
pcaobj2 <- pcafd( fdaobjMale, nharm=2 )
pcaobj3 <- pcafd( fdaobjMale, nharm=3, Plot=TRUE)


## PCA - PCA Functional Data Analysis
pcaobj2 <- pcafd( fdaNox, nharm=4, Plot=TRUE  )


###################################################
# Contour Plot - Density Two Principal Components
##################################################

CPlot<- densityScores(pcaobj2, bivariate=FALSE)
CPlot<- densityScores(pcaobj2, bivariate=TRUE)

#############################
# Kernell Density Estimation 
#############################

fdvector<-fdensity( fdaobjMale, pcaobj2, bandwith=900)
plot(fdvector, type="p")

###################
# Discriminating 
####################

par(mfrow=c(1,2), font.lab=3, font.axis=1) 
matplot(Gwdprotein$Time, Gwdprotein$Group1, type="l", main="Wheat Data Set", sub="Group 1", xlab="Wave length(nm)", ylab=" d(NIR spectra ( log(1/R) ),1)")
abline(h= 0)

matplot(Gwdprotein$Time, Gwdprotein$Group2, type="l", main="Wheat Data Set", sub="Group 2", xlab="Wave length(nm)", ylab=" d(NIR spectra ( log(1/R) ),1) ")
abline(h= 0)

tr <- sample(1:100, 50)
AllGwdprotein<-cbind(Gwdprotein$Group1,Gwdprotein$Group2)
train<- AllGwdprotein[,tr]
test<- AllGwdprotein[,-tr]
Classlearn <- sort(rep(1:2,25)) 
Classifmplsr <- classfd(Classlearn, train, test)

##############################
# Perfect FDA  Discriminating 
##############################


AllGwdprotein<-cbind(Gwdprotein$Group1,Gwdprotein$Group2)
fdaobjProteine<-fdaclass(AllGwdprotein,Gwdprotein$Time)
XXdata<-t(fdaobjProteine$data)
nall=dim(XXdata)[1]

Classes=c(rep(0,50),rep(1,50))

indClass1=which(Classes==0)
indClass0=which(Classes==1)

indClass1All=indClass1
indClass0All=indClass0

indClass1=1:length(indClass1)
indClass0=1:length(indClass0)+length(indClass1)

ntest=30
s=sample(c(0,1), 30, replace = TRUE)
indtest=(1:nall)[s==1]
indClass1=setdiff(indClass1All,indtest)
indClass0=setdiff(indClass0All,indtest)

indtest=sort(indtest)
ntrain=nall-length(indtest)

as.array(indtest)
XXtest=XXdata[indtest,]
XX=XXdata[-indtest,]
Classes<-c(indClass1,indClass0)

fd3<-pdfclasf(XXdata, XXtest, indClass0, indClass1, indtest )
plot(-fd3$Hatpsir, type="l")
fd3$Missclassification

##############################
# Simulate FDA Data 
##############################

simulation1<- simulatefda( nsamples=10, ndrawn=100, rangeval= c(-10,10), mean=0, sigma=1 )
simulation2<- simulatefda( nsamples=50, ndrawn=100, rangeval= c(-5,5), mean=3, sigma=1 )

n<-150
hues = seq(15, 375, length=n+1)
color<-hcl(h=hues, l=65, c=100)[1:n]

par(mfrow = c(2,1))

matplot(simulation1$rangetime,simulation1$simulation,type="l",col =color,xlab="",ylab="x(t)",main="Observations",lty=1)
matplot(simulation2$rangetime,simulation2$simulation,type="l",col =color,xlab="",ylab="x(t)",main="Observations",lty=1)

