#######
# If needed
library(devtools)
install_github("ccrismancox/games2")
######
library(games2)
library(brglm) 
library(detectseparation)
library(numDeriv)
library(readstata13)
library(ggplot2)
library(gridExtra)
library(stringr)
library(knitr)

rm(list=ls())
source("extraFunctions.R")
huth <- read.dta13("../data/huth.dta")

SBI.Sig <- function(dat,names, i, br=FALSE){ #SBI function for bootstrap
  fn <- ifelse(br,
               brglm.fit,
               glm.fit)
  if(br){
    start <- matrix(0,nrow=9)
  }else{
    start <- matrix(s1$coef[13:21], nrow=9)
  }
  
  BS <- dat[i,]
  Y1 <- BS[,1]
  Y2 <- BS[,2][Y1==0]
  BS11 <- cbind(1, BS[,names$U11])
  BS14 <- BS[,names$U14]
  BS2 <- cbind(1,BS[,names$U2])/sqrt(2)
  BS13 <- 1
  mB <- convTest(Y=Y2,X=BS2[Y1==0,], fn=fn)
  if(anyNA(mB) | length(mB) < ncol(BS2)){
    mA <- NA
  }else{
    p4hat <- as.numeric(pnorm(BS2 %*% mB)) #sqrt 2 in the data
    Z11hat <- p4hat *BS14
    Z10hat <- (1-p4hat)*BS13
    mA <- convTest(Y=Y1,X=cbind(BS11, -Z10hat, -Z11hat)/sqrt(1+p4hat^2 + (1-p4hat)^2),
                   fn=fn)
  }
  return(c(mA, mB))
}

#### FIML #### 
s1 <- egame12(outcome~tft+firmflex+dem_att+syear|1|nuclear+ibf+sbf+lbf+milallia+milarm-1|nuclear+ibf+sbf+milallia+milarm+fortrade+stalemat+dem_def,
              data = huth, link = "probit", type = "private",startvals="zero")
#save the data
dat <- s1$model

# create Y
Y1 <- ifelse(dat$outcome==1,
             1,
             0)

Y2 <- ifelse(Y1==0 & dat$outcome==4,
             1,
             0)

dat <- cbind(Y1, Y2, dat)
dat$outcome <- NULL
dat <-  as.matrix(dat)

# sort names
names <-list(U11=c("tft",
                   "firmflex",
                   "dem_att",
                   "syear"),
             U14=c("nuclear",
                   "ibf",
                   "sbf",
                   "lbf",
                   "milallia",
                   "milarm"),
             U2=c("nuclear",
                  "ibf",
                  "sbf",
                  "milallia",
                  "milarm",
                  "fortrade",
                  "stalemat",
                  "dem_def"))







#### SBI- No Correction####

Y1 <- dat[,1]
Y2 <- dat[,2][Y1==0]
dat11 <- cbind(1, dat[,names$U11])
dat14 <- dat[,names$U14]
dat2 <- cbind(1,dat[,names$U2])/sqrt(2)
dat13 <- 1
colnames(dat2)[1] <- "intercept"
dat2part <-  dat2[Y1==0,]
mf.sbi <- as.data.frame(dat2part)
mf.sbi$Y2 <- Y2
mB <- glm(Y2~intercept+nuclear+ibf+sbf+milallia+
            milarm+fortrade+stalemat+dem_def-1, 
          data=mf.sbi,
          x=TRUE,
          family=binomial(link="probit"))
mf.profile <- mf.sbi
car::vif(mB)

p4hat <- as.numeric(pnorm(as.matrix(dat2) %*% mB$coef)) #root 2 is already in the data
Z11hat <- p4hat *dat14
Z10hat <- (1-p4hat)*dat13
dat1 <- cbind(dat11, -Z10hat, -Z11hat)/sqrt(1+p4hat^2 + (1-p4hat)^2)
colnames(dat1) <- c("intercept_sq", names$U11, "itnerceptBD", names$U14)
mf.sbi <- as.data.frame(dat1)
mf.sbi$Y1 <- Y1
mA <- glm(Y1~intercept_sq+tft+firmflex+dem_att+syear+itnerceptBD+nuclear+ibf+sbf+lbf+milallia+milarm-1,
          x=TRUE, data=mf.sbi,
          family=binomial(link="probit"))

fullDesign <- cbind.data.frame(dat1, dat2); colnames(fullDesign) <- paste("V", 1:ncol(fullDesign), sep="")


Bsep.lp <- detect_separation(y=Y2, x=mB$x,
                             family=binomial("probit"), intercept = F, 
                             control=list(linear_program = "primal", purpose="find"))
Asep.lp <- detect_separation(y=Y1, x=mA$x,
                             family=binomial("probit"), intercept = F, 
                             control=list(linear_program = "primal", purpose="find"))

# Try and compute standard errors
OmegaB.inv <- vcov(mA)

U11 <- dat11 %*% mA$coefficients[1:5] #untransformed X
U13 <- dat13 * mA$coefficients[6] #untransformed X
U14 <- dat14 %*% mA$coefficients[7:12] #untransformed X
normalize <- 1+p4hat^2 + (1-p4hat)^2
interior <- (U11-U13*(1-p4hat)-U14*p4hat)/sqrt(normalize)
chain <- (U13-U14)/sqrt(normalize) + ((1-2*p4hat)*((U11-U13*(1-p4hat)-U14*p4hat)))/(normalize^(3/2))
dLdP <- diag(drop(Y1*((dnorm(interior)/pnorm(interior))*( chain)) -
                    (1-Y1)*((dnorm(interior)/(pnorm(interior,lower=F)))*(chain))))
dLdB <-  drop(Y1*dnorm(interior)/pnorm(interior) - (1-Y1)*dnorm(interior)/pnorm(interior,lower=F)) * cbind(dat11,-dat13*(1-p4hat), -dat14*(p4hat))/sqrt(1+p4hat^2 + (1-p4hat)^2)


Omegap <- crossprod(dLdB, dLdP)
SIGMA <- drop(dnorm(as.matrix(dat2) %*% mB$coef)^2) *(dat2 %*% vcov(mB) %*% t(dat2)) # dat2 was already transofrmed
SSE2<- c(sqrt(diag(OmegaB.inv + OmegaB.inv %*% Omegap %*% SIGMA %*% t(Omegap) %*% OmegaB.inv)), sqrt(diag(vcov(mB))))
cat("NAs detectedin SBI Standard Errors, attempting a bootstrap:\t", anyNA(SSE2))

## Moving to bootstrap
s1SBI <- matrix(0, nrow=50, ncol=length(s1$coef))
set.seed(1)
for(i in 1:50){
  c1 <- NA
  while(anyNA(c1)||any(abs(c1)>50)){
    j <- sample(1:nrow(dat), replace=TRUE)
    br=FALSE
    c1 <- SBI.Sig(dat, names, j, br)
  }
  s1SBI[i,] <- c1 #Uncorrected
}    
SBI <- list(mB=mB, mA=mA)



#### BR-SBI####

Y1 <- dat[,1]
Y2 <- dat[,2][Y1==0]
dat11 <- cbind(1, dat[,names$U11])
dat14 <- dat[,names$U14]
dat2 <- cbind(1,dat[,names$U2])/sqrt(2)
dat13 <- 1
colnames(dat2)[1] <- "intercept"
dat2part <-  dat2[Y1==0,]
mf.sbi <- as.data.frame(dat2part)
mf.sbi$Y2 <- Y2
mB <- brglm(Y2~intercept+nuclear+ibf+sbf+milallia+
              milarm+fortrade+stalemat+dem_def-1, 
            data=mf.sbi, intercept=F,
            family=binomial(link="probit"),
            pl=TRUE)

p4hat <- as.numeric(pnorm(as.matrix(dat2) %*% mB$coef))
Z11hat <- p4hat *dat14
Z10hat <- (1-p4hat)*dat13
dat1 <- cbind(dat11, -Z10hat, -Z11hat)/sqrt(1+p4hat^2 + (1-p4hat)^2)
colnames(dat1) <- c("intercept_sq", names$U11, "itnerceptBD", names$U14)
mf.sbi <- as.data.frame(dat1)
mf.sbi$Y1 <- Y1
mA <- brglm(Y1~intercept_sq+tft+firmflex+dem_att+syear+itnerceptBD+nuclear+ibf+sbf+lbf+milallia+milarm-1,
            data=mf.sbi, intercept=F,
            family=binomial(link="probit"),
            pl=TRUE)



OmegaB.inv <- vcov(mA)

U11 <- dat11 %*% mA$coefficients[1:5]
U13 <- dat13 * mA$coefficients[6]
U14 <- dat14 %*% mA$coefficients[7:12]
normalize <- 1+p4hat^2 + (1-p4hat)^2
interior <- (U11-U13*(1-p4hat)-U14*p4hat)/sqrt(normalize)
chain <- (U13-U14)/sqrt(normalize) + ((1-2*p4hat)*((U11-U13*(1-p4hat)-U14*p4hat)))/(normalize^(3/2))
dLdP <- diag(drop(Y1*((dnorm(interior)/pnorm(interior))*( chain)) -
                    (1-Y1)*((dnorm(interior)/(pnorm(interior,lower=F)))*(chain))))
dLdB <-  drop(Y1*dnorm(interior)/pnorm(interior) - (1-Y1)*dnorm(interior)/pnorm(interior,lower=F)) *
  cbind(dat11,-dat13*(1-p4hat), -dat14*(p4hat))/sqrt(1+p4hat^2 + (1-p4hat)^2)


Omegap <- crossprod(dLdB, dLdP)
SIGMA <- drop(dnorm(as.matrix(dat2) %*% mB$coef)^2) *(dat2 %*% vcov(mB) %*% t(dat2))
SSE4 <- c(sqrt(diag(OmegaB.inv + OmegaB.inv %*% Omegap %*% SIGMA %*% t(Omegap) %*% OmegaB.inv)), sqrt(diag(vcov(mB))))
Scoefs4 <- c(mA$coef, mB$coef)
BR.SBI <- list(mB=mB, mA=mA)


#### BR-FIML ####
s2 <- try(egame12(outcome~tft+firmflex+dem_att+syear|
                    1|
                    nuclear+ibf+sbf+lbf+milallia+milarm-1|
                    nuclear+ibf+sbf+milallia+milarm+fortrade+stalemat+dem_def,
                  data = huth,
                  link = "probit", 
                  type = "private",
                  startvals="zero",
                  penalty="logF",
                  print.level=1,
                  iterlim=100, w=2,
                  method="BFGS"))





##Post estimation FIML test#
Y1 <- dat[,1]
Y2 <- dat[,2][Y1==0]
dat11 <- cbind(1, dat[,names$U11])
dat14 <- dat[,names$U14]
dat2 <- cbind(1,dat[,names$U2])/sqrt(2)
p4hat <- predict(s1, type="action")[,4]
Z11hat <- p4hat *dat14
Z10hat <- (1-p4hat)*dat13
dat1 <- cbind(dat11, -Z10hat, -Z11hat)
colnames(dat1) <- c("int1", names$U11, "int3", names$U14)
fullDesign <- cbind.data.frame(dat1, dat2); colnames(fullDesign) <- paste("V", 1:ncol(fullDesign), sep="")

total.lp0  <- detect_separation(x=fullDesign, ((1-dat[,1]) + ((1-dat[,1])*dat[,2]))==0,
                                family=binomial("probit"), intercept = F, 
                                control=list(linear_program = "primal", purpose="find"))
total.lp1  <- detect_separation(x=fullDesign, ((1-dat[,1]) + ((1-dat[,1])*dat[,2]))==1,
                                family=binomial("probit"), intercept = F, 
                                control=list(linear_program = "primal", purpose="find"))
total.lp2  <- detect_separation(x=fullDesign, ((1-dat[,1]) + ((1-dat[,1])*dat[,2]))==2,
                                family=binomial("probit"), intercept = F, 
                                control=list(linear_program = "primal", purpose="find"))



tab3 <- data.frame(Regressors = c("X[B]", "Z[SBI]", rep("Z[FIML] & X[B]", 3)),
                   Outcome = c("y[B] given y[A] = 1",
                               "y[A] = ",
                               "y[A] = 0",
                               'y[A] = y[B] = 0',
                               'y[A] = 0 & y[B] = 1'),
                   Results=c(Bsep.lp$separation,
                             Asep.lp$separation,
                             total.lp0$separation,
                             total.lp1$separation,
                             total.lp2$separation))
cat("Table 3 \n", file="../tables_and_figures/Table3.md")
cat(kable(tab3),
    sep="\n",
    file="../tables_and_figures/Table3.md",
    append=TRUE)

#### Check that BR-SBI doesn't change with Cauchy or logF penalty ####

brprobit <- function(b, X, y, penalty=c("Cauchy", "logF")){
  penalty <- match.arg(penalty)
  if(penalty=="Cauchy"){
    scale <- ifelse(str_detect(names(b), regex("Intercept", ignore_case=T)), 10,2.5)
    penalty <- sum(-log((1+ (b/scale)^2)))
  }else{
    penalty <- sum(b/2 -log(1+exp(b)))
  }
  LL <- sum(pnorm((2*y-1)*X%*%b, log.p=TRUE))+penalty
  return(LL)
}
probit <- function(b, X, y){
  LL <- sum(pnorm((2*y-1)*X%*%b, log.p=TRUE))
  return(-LL)
}
startB <- 0*mB$coef

startA <- 0*mA$coef


# cauchy penalty
Y1 <- dat[,1]
Y2 <- dat[,2][Y1==0]
dat11 <- cbind(1, dat[,names$U11])
dat14 <- dat[,names$U14]
dat2 <- cbind(1,dat[,names$U2])/sqrt(2)
dat13 <- 1
dat2part <-  dat2[Y1==0,]
mB <- maxLik(brprobit,start=startB, y=Y2, X=dat2part,method="BFGS")
vB <- solve(hessian(probit, mB$est,y=Y2, X=dat2part))
p4hat <- as.numeric(pnorm(as.matrix(dat2) %*% mB$est))
Z11hat <- p4hat *dat14
Z10hat <- (1-p4hat)*dat13
dat1 <- cbind(dat11, -Z10hat, -Z11hat)/sqrt(1+p4hat^2 + (1-p4hat)^2)
colnames(dat1) <- c("int1", names$U11, "int3", names$U14)
mA <- maxLik(brprobit, start=startA, y=Y1, X=dat1, method="BFGS")



OmegaB.inv <- vcov(mA)

U11 <- dat11 %*% mA$estimate[1:5]
U13 <- dat13 * mA$estimate[6]
U14 <- dat14 %*% mA$estimate[7:12]
normalize <- 1+p4hat^2 + (1-p4hat)^2
interior <- (U11-U13*(1-p4hat)-U14*p4hat)/sqrt(normalize)
chain <- (U13-U14)/sqrt(normalize) + ((1-2*p4hat)*((U11-U13*(1-p4hat)-U14*p4hat)))/(normalize^(3/2))
dLdP <- diag(drop(Y1*((dnorm(interior)/pnorm(interior))*( chain)) -
                    (1-Y1)*((dnorm(interior)/(pnorm(interior,lower=F)))*(chain))))
dLdB <-  drop(Y1*dnorm(interior)/pnorm(interior) - (1-Y1)*dnorm(interior)/pnorm(interior,lower=F)) * cbind(dat11,-dat13*(1-p4hat), -dat14*(p4hat))/sqrt(1+p4hat^2 + (1-p4hat)^2)


Omegap <- crossprod(dLdB, dLdP)
SIGMA <- drop(dnorm(as.matrix(dat2) %*% mB$estimate)^2) *(dat2 %*% vB %*% t(dat2))
SSE4c <- c(sqrt(diag(OmegaB.inv + OmegaB.inv %*% Omegap %*% SIGMA %*% t(Omegap) %*% OmegaB.inv)), sqrt(diag(vB)))
Scoefs4c <- c(mA$est, mB$est)

# log F

Y1 <- dat[,1]
Y2 <- dat[,2][Y1==0]
dat11 <- cbind(1, dat[,names$U11])
dat14 <- dat[,names$U14]
dat2 <- cbind(1,dat[,names$U2])/sqrt(2)
dat13 <- 1
dat2part <-  dat2[Y1==0,]
mB <- maxLik(brprobit,start=startB, y=Y2, X=dat2part,method="BFGS", penalty='logF')
vB <- solve(hessian(probit, mB$est,y=Y2, X=dat2part))
p4hat <- as.numeric(pnorm(as.matrix(dat2) %*% mB$est))
Z11hat <- p4hat *dat14
Z10hat <- (1-p4hat)*dat13
dat1 <- cbind(dat11, -Z10hat, -Z11hat)/sqrt(1+p4hat^2 + (1-p4hat)^2)
colnames(dat1) <- c("int1", names$U11, "int3", names$U14)
mA <- maxLik(brprobit, start=startA, y=Y1, X=dat1, method="BFGS", penalty='logF')



OmegaB.inv <- vcov(mA)

U11 <- dat11 %*% mA$estimate[1:5]
U13 <- dat13 * mA$estimate[6]
U14 <- dat14 %*% mA$estimate[7:12]
normalize <- 1+p4hat^2 + (1-p4hat)^2
interior <- (U11-U13*(1-p4hat)-U14*p4hat)/sqrt(normalize)
chain <- (U13-U14)/sqrt(normalize) + ((1-2*p4hat)*((U11-U13*(1-p4hat)-U14*p4hat)))/(normalize^(3/2))
dLdP <- diag(drop(Y1*((dnorm(interior)/pnorm(interior))*( chain)) -
                    (1-Y1)*((dnorm(interior)/(pnorm(interior,lower=F)))*(chain))))
dLdB <-  drop(Y1*dnorm(interior)/pnorm(interior) - (1-Y1)*dnorm(interior)/pnorm(interior,lower=F)) * cbind(dat11,-dat13*(1-p4hat), -dat14*(p4hat))/sqrt(1+p4hat^2 + (1-p4hat)^2)


Omegap <- crossprod(dLdB, dLdP)
SIGMA <- drop(dnorm(as.matrix(dat2) %*% mB$estimate)^2) *(dat2 %*% vB %*% t(dat2))
SSE4f <- c(sqrt(diag(OmegaB.inv + OmegaB.inv %*% Omegap %*% SIGMA %*% t(Omegap) %*% OmegaB.inv)), sqrt(diag(vB)))
Scoefs4f <- c(mA$est, mB$est)

# looking at hte SBI differences across penalties
cbind(Scoefs4, Scoefs4c,Scoefs4f)
cbind(SSE4, SSE4c,SSE4f)
# not too many


# #### BR-FIML Cauchy comparison for interested parties
# s3 <- try(egame12(outcome~tft+firmflex+dem_att+syear|
#                     1|
#                     nuclear+ibf+sbf+lbf+milallia+milarm-1|
#                     nuclear+ibf+sbf+milallia+milarm+fortrade+stalemat+dem_def,
#                   data = huth,
#                   link = "probit", 
#                   type = "private",
#                   startvals="zero",
#                   penalty="Cauchy",
#                   print.level=1,
#                   iterlim=100, w=2,
#                   method="BFGS"))
# 
# car::compareCoefs(s1,s2,s3)
# 

#### Table it ####

source("xtable2.R")

Scoefs <- list(s1$coef,
               apply(s1SBI, 2, median),
               s2$coef,
               Scoefs4)
SSE  <-  list(sqrt(diag(s1$vcov)),
              apply(s1SBI, 2, sd),
              sqrt(diag(s2$vcov)),
              SSE4)



names <- c("Const.",
           "Tit-for-Tat",
           "Firm-Flex",
           "Democratic Attacker",
           "Year",
           "Const.",
           "Nuclear",
           "Immediate Balance",
           "Short-term Balance",
           "Long-term Balance",
           "Military Alliance",
           "Arms Transfers",
           "Const.",
           "Nuclear",
           "Immediate Balance",
           "Short-term Balance",
           "Military Alliance",
           "Arms Transfers",
           "Foreign Trade",
           "Stalemate",
           "Democratic Defender")
names <- paste(c(rep("$U_A$(SQ):", 5),
                 "$U_A$(BD):",
                 rep("$U_A$(War):", 6),
                 rep("$U_B$(War):",9)), names) 


info.list <- list(UFMLE=list(logLik = logLik(s1),
                             Obs =  nrow(s1$model)),
                  SBI = list(logLik = "",
                             Obs =  nrow(s1$model)),
                  CFMLE=list(logLik = logLik(s2),
                             Obs =  nrow(s2$model)),
                  CSBI =list(logLik ="",
                             Obs =  nrow(s2$model)))

model.names <- c("FMLE",
                 "SBI",
                 "BR FMLE",
                 "BR SBI")

formattedCoef <- lapply(Scoefs, round, 2)
formattedSE<- lapply(SSE, function(x){paste0("(", formatC(x,format ="f",digits=2),")")})
tabOut <- matrix("", nrow=length(Scoefs[[1]])*2, ncol=4)
for(i in 1:4){tabOut[,i] <- matrix(rbind(formattedCoef[[i]], formattedSE[[i]]), ncol=1)}
tabOut <- cbind(matrix(rbind(names, ""), ncol=1), tabOut)
colnames(tabOut) <- c("", model.names)
tabOut <- rbind(tabOut, c("Obseravations", rep(nrow(huth), 4)))
cat("Table 4\n", file="../tables_and_figures/Table4.md")
cat(kable(tabOut),
    sep="\n",
    file="../tables_and_figures/Table4.md",
    append = TRUE)


save.image("sigTarar_output.rdata")




#### Plotting ####
mf.sbi <- mf.profile

ps1 <- profile(s1,which=c(17), dist=1.75, step=50)
ps1.2 <- profileModel(SBI$mB, which=5, stdn=1.5, gridsize=10, objective="ordinaryDeviance")
ps2 <- profile(s2,which=c(17), dist=1.75, step=50)
ps2.2 <- profile(BR.SBI$mB, which=5,  stdn=5, stepsize=3,gridsize=50)
names(ps1) <-  "Military Alliance"
ps1 <- lapply(ps1,
              function(x){
                colnames(x) <- c("Profiled Likelihood",
                                 names)
                return(x)
              })

p1temp <- ps1[[1]][,c(17+1,1)]

p2temp <- ps2[[1]][,c(17+1,1)]
p22temp <- data.frame(milallia = ps2.2$profilesBR$profiles$milallia[,"milallia"],
                      logLik = -ps2.2$profilesBR$profiles$milallia[,"Differences"] + logLik(SBI$mB))
p22temp <- p22temp[p22temp$logLik > -10,]

p1df <- data.frame(var = p1temp[,1],
                   LL = p1temp[,2])
p2df <- data.frame(var = p2temp[,1],
                   LL = p2temp[,2])
p22df <- data.frame(var = p22temp[,1],
                    LL = p22temp[,2])


plot.fiml <- ggplot(p1df)+
  geom_line(aes(x=var,
                y=LL))+
  geom_vline(aes(xintercept=s1$coef[17]), alpha=.5)+
  xlab(expression(U[B]:~"Military Alliance"))+
  ylab("Profiled log likelihood")+
  theme_bw(8)+
  ggtitle("Ordinary FIML")

plot.brfiml <- ggplot(p2df)+
  geom_line(aes(x=var,
                y=LL))+
  geom_vline(aes(xintercept=s2$coef[17]), alpha=.5)+
  xlab(expression(U[B]:~"Military Alliance"))+
  ylab("Profiled log likelihood")+
  theme_bw(8)+
  ggtitle("BR-FIML")

plot.brsbi <- ggplot(p22df)+
  geom_line(aes(x=var,
                y=LL))+
  geom_vline(aes(xintercept=Scoefs[[4]][17]), alpha=.5)+
  xlab(expression(U[B]:~"Military Alliance"))+
  ylab("Profiled log likelihood")+
  theme_bw(8)+
  ggtitle("BR-SBI")
pdf(width=3, height=5,# pointsize=10, 
    file="../tables_and_figures/figure3.pdf")
grid.arrange(plot.fiml, plot.brfiml, plot.brsbi, ncol=1)
dev.off()





save(list=c("s2", "dat", "dat2"), file="signorinotarar_brfit.rdata")
