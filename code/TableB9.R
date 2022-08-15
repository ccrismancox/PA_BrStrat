#######
# If needed
library(devtools)
install_github("ccrismancox/games2")
######
library(mc2d)
library(games2)
library(doParallel)
library(doRNG)
library(brglm)
library(detectseparation)
library(matrixStats)
library(readstata13)
library(data.table)
library(stringr)
rm(list=ls())

load("signorinotarar_brfit.rdata")
huth <- read.dta13("../data/huth.dta")
rmse <- function(bias, var){
  return(sqrt(sum(var)+ crossprod(bias)))
}

s2.fixed <- s2
dat2.fixed <- dat2
dat.fixed <- dat
s2.fixed$coefficients['u2(4):milallia'] <- 0

reject <- c(TRUE,TRUE, rep(FALSE,7))
workers <- makeCluster(25) 
registerDoParallel(workers)
set.seed(2)
out <- foreach(i =  1:500, .packages=c("games2", "brglm", "detectseparation", "mc2d"),
               .combine=rbind, .errorhandling = "remove")%dorng%{
                 
                 sep=FALSE;s2 <- list(coef=rep(200, 21))
                 #try to keep things from totally jumping the rails; stability condition for a small sample with shakey data
                 while( sep==FALSE | max(abs(s2$coef)) > 100){ 
                   probs <- predict(s2.fixed, type="outcome")
                   Y <- rmultinomial(n=58, size=1, as.matrix(probs))
                   
                   outcomeMC <- Y[,1] + Y[,2]*3 + Y[,3]*4
                   sep <- detect_separation(x=dat2.fixed[Y[,1]==0 ,],y=Y[Y[,1]==0 ,][,3], family=binomial("probit"), intercept = FALSE)$separation
                   
                   s2 <- try(egame12(outcomeMC~tft+firmflex+dem_att+syear|1|nuclear+ibf+sbf+lbf+milallia+milarm-1|
                                       nuclear+ibf+sbf+milallia+milarm+fortrade+stalemat+dem_def, 
                                     w=1,
                                     data = huth,
                                     link = "probit", 
                                     type = "private",
                                     startvals=s2.fixed$coefficients, 
                                     penalty="logF",
                                     iterlim=100, method="BFGS",
                                     print.level=1)) 
                   if(class(s2)[1]!="game" || s2$conv$code !=0){ s2 <- list(coef=rep(200, 21))}
                 }
                 
                 
                 
                 s1 <- try(egame12(outcomeMC~tft+firmflex+dem_att+syear|1|nuclear+ibf+sbf+lbf+milallia+milarm-1|nuclear+ibf+sbf+milallia+milarm+fortrade+stalemat+dem_def,
                                   data = huth,
                                   link = "probit", 
                                   type = "private",
                                   startvals=s2.fixed$coefficients, 
                                   iterlim=500, method="BFGS",
                                   print.level=0)) 
                 if(class(s1)[1]=="game"){
                   NCFIML <- s1$coefficients
                 }else{
                   NCFIML <- rep(NA, 21)
                 }
                 
                 if(class(s2)[1]=="game"){
                   WCFIML <- s2$coefficients
                 }else{
                   WCFIML <- rep(NA, 21)
                 }
                 
                 s3 <- try(egame12(outcomeMC~tft+firmflex+dem_att+syear|1|nuclear+ibf+sbf+lbf+milallia+milarm-1|nuclear+ibf+sbf+milallia+milarm+fortrade+stalemat+dem_def,
                                   data = huth,
                                   link = "probit", 
                                   type = "private",
                                   penalty="Cauchy",
                                   startvals=s2.fixed$coefficients, 
                                   iterlim=500, method="BFGS",
                                   print.level=0)) 
                 if(class(s3)[1]=="game"){
                   WCFIML1 <- s3$coefficients
                 }else{
                   WCFIML1 <- rep(NA, 21)
                 }
                 
                 dat <- as.data.frame(dat.fixed)
                 dat$outcomeMC <- outcomeMC
                 
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
                 dat$Y1 <- ifelse(dat$outcomeMC==1,
                                  1,
                                  0)
                 
                 dat$Y2 <- ifelse(dat$Y1==0 & dat$outcomeMC==4,
                                  1,
                                  0)
                 
                 
                 dat$outcomeMC <- NULL
                 dat <-  as.matrix(dat)
                 
                 Y1 <- dat[,1]
                 Y2 <- dat[,2][Y1==0]
                 dat11 <- cbind(1, dat[,names$U11])
                 dat14 <- dat[,names$U14]
                 dat2 <- cbind(1,dat[,names$U2])/sqrt(2)
                 dat13 <- 1
                 colnames(dat2)[1] <- "int"
                 dat2part <-  dat2[Y1==0,]
                 
                 mB <- try(glm(Y2~dat2part-1, 
                               x=TRUE, start=s2.fixed$coefficients[13:21],
                               family=binomial(link="probit")))
                 if("glm" %in% class(mB)){
                   p4hat <- as.numeric(pnorm(as.matrix(dat2) %*% mB$coef)) #sqrt 2 are factored in
                   Z11hat <- p4hat *dat14
                   Z10hat <- (1-p4hat)*dat13
                   dat1 <- cbind(dat11, -Z10hat, -Z11hat)/sqrt(1+p4hat^2 + (1-p4hat)^2)
                   colnames(dat1) <- c("int1", names$U11, "int3", names$U14)
                   mA <- try(glm(Y1~dat1-1,x=TRUE,start=s2.fixed$coefficients[1:12],
                                 family=binomial(link="probit")))
                   if("glm" %in% class(mA)){
                     NCSBI <- c(mA$coef, mB$coef)
                   }else{
                     NCSBI <- rep(NA, 21)
                   }
                 }else{
                   NCSBI <- rep(NA, 21)
                 }
                 
                 ## BR
                 Y1 <- dat[,1]
                 Y2 <- dat[,2][Y1==0]
                 dat11 <- cbind(1, dat[,names$U11])
                 dat14 <- dat[,names$U14]
                 dat2 <- cbind(1,dat[,names$U2])/sqrt(2)
                 dat13 <- 1
                 dat2part <-  dat2[Y1==0,]
                 mB2 <- try(brglm(Y2~dat2part-1, 
                                 x=TRUE, intercept=F,start=s2.fixed$coefficients[13:21],
                                 family=binomial(link="probit"),
                                 pl=TRUE))
                 if("glm" %in% class(mB2)){
                   p4hat <- as.numeric(pnorm(as.matrix(dat2) %*% mB2$coef))
                   Z11hat <- p4hat *dat14
                   Z10hat <- (1-p4hat)*dat13
                   dat1 <- cbind(dat11, -Z10hat, -Z11hat)/sqrt(1+p4hat^2 + (1-p4hat)^2)
                   colnames(dat1) <- c("int1", names$U11, "int3", names$U14)
                   mA <- try(brglm(Y1~dat1-1,x=TRUE,start=s2.fixed$coefficients[1:12],
                                   family=binomial(link="probit"), pl=TRUE))
                   if("glm" %in% class(mA)){
                     WCSBI <- c(mA$coef, mB2$coef)
                   }else{
                     WCSBI <- rep(NA, 21)
                   }
                   
                 }else{
                   WCSBI <- rep(NA, 21)
                 }
                 Htests.SBINC <- (summary(mB)$coef[,4] < 0.05) == reject
                 Htests.SBIWC <- (summary(mB2)$coef[,4] < 0.05) == reject
                 Htests.FIMLNC <- (summary(s1)$coef[,4][13:21] < 0.05) ==reject
                 Htests.FIMLWC <- (summary(s2)$coef[,4][13:21] < 0.05) ==reject
                 Htests.FIMLWC1 <- (summary(s3)$coef[,4][13:21] < 0.05) ==reject
                 
                 

                 if(any(abs(NCSBI) > 1000,na.rm=T)){
                   NCSBI <- rep(NA, length(NCSBI))
                   Htests.SBINC <- rep(NA,9)
                 }
                 if(any(abs(WCSBI) > 1000,na.rm=T)){
                   WCSBI <- rep(NA, length(WCSBI))
                   Htests.SBIWC <- rep(NA,9)
                 }
                 if(any(abs(NCFIML) > 1000,na.rm=T)){
                   NCFIML <- rep(NA, length(NCFIML))
                 }
                 if(any(abs(WCFIML) > 1000,na.rm=T)){
                   WCFIML <- rep(NA, length(WCFIML))
                 }
                 if(any(abs(WCFIML1) > 1000,na.rm=T)){
                   WCFIML1 <- rep(NA, length(WCFIML1))
                   Htests.FIMLWC1 <- rep(NA,9)
                 }
                 
                 
                 output <- c(NCSBI, WCSBI, NCFIML,WCFIML,WCFIML1,
                             Htests.SBINC, Htests.SBIWC,Htests.FIMLNC,
                             Htests.FIMLWC,Htests.FIMLWC1)
                 
                 output
               }
stopCluster(workers)
save.image("tableB9.rdata")



load("tableB9.rdata")
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


correctOutcomes <- matrix(colMeans(out[,106:ncol(out)], na.rm=T), ncol=5) #first 2 are power everything else is 1- type 1
rownames(correctOutcomes) <- names[str_detect(names, "U_B")]
colnames(correctOutcomes) <- c("Ordinary SBI",
                               "BR-SBI",
                               "Ordinary FIML",
                               "BR-FIML (Cauchy)",
                               "BR-FIML (Log-F)")
output <- rbind(correctOutcomes[2,,drop=FALSE], # power for Nuclear
                1-correctOutcomes[5,,drop=F]) #type 1 for alliance
cat("Table B.9\n",    file = "../tables_and_figures/TableB9.md")
cat(kable(output, digits=2), 
    sep="\n",
    file = "../tables_and_figures/TableB9.md",
    append = TRUE)

