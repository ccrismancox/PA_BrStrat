#######
# If needed
library(devtools)
#install_github("ccrismancox/games2")
######
library(games2)
library(brglm) 
library(detectseparation)
library(knitr)
rm(list=ls())
source("extraFunctions.r")
### Load the data (saved in the games package)
data(leblang2003)

## Fit the FIML
l1 <- egame12(outcome ~ capcont + lreserves + overval + creditgrow +
                USinterest + service + contagion + prioratt - 1 | 1 | 1 | unifgov + lexports
              + preelec + postelec + rightgov + realinterest + capcont + lreserves, data =
                leblang2003, link = "probit", type = "private")

dat <- l1$model # save the relevant data

names <-list(U1=c("capcont",
                  "lreserves",
                  "overval",
                  "creditgrow",
                  "USinterest",
                  "service",
                  "contagion",
                  "prioratt"),
             U2=c("unifgov",
                  "lexports",
                  "preelec",
                  "postelec",
                  "rightgov",
                  "realinterest",
                  "capcont",
                  "lreserves"))
Y1 <- ifelse(dat$outcome=="no attack",
             1,
             0)

Y2 <- ifelse(Y1==0 & dat$outcome=="defense",
             1,
             0)

dat <- cbind(Y1, Y2, dat)



dat$outcome <- NULL
dat <-  as.matrix(dat)

#### SBI - NC####
Y1 <- dat[,1]
Y2 <- dat[,2][Y1==0]
dat11 <- dat[,names$U1]
dat14 <- 1
dat2 <- cbind(1,dat[,names$U2])/sqrt(2)
colnames(dat2)[1] <- "int"
dat13 <- 1
dat2part <-  dat2[Y1==0,]
mB <- glm(Y2~dat2part-1, 
          family=binomial(link="probit"), x=TRUE)
p4hat <- as.numeric(pnorm(as.matrix(dat2) %*% mB$coef)) #sqrt(2) already in the data
Z11hat <- p4hat *dat14
Z10hat <- (1-p4hat)*dat13
dat1 <- cbind(dat11, -Z10hat, -Z11hat)/sqrt(1+p4hat^2 + (1-p4hat)^2)
colnames(dat1) <- c(names$U1, "int3", "int4")
mA <- glm(Y1~dat1-1,x=TRUE,
          family=binomial(link="probit"))
fullDesign <- cbind.data.frame(dat1, dat2); colnames(fullDesign) <- paste("V", 1:ncol(fullDesign), sep="")


Bsep.lp <- detect_separation(y=Y2, x=mB$x,
                             family=binomial("probit"), intercept = F, 
                             control=list(linear_program = "primal", purpose="find"))
Asep.lp <- detect_separation(y=Y1, x=mA$x,
                             family=binomial("probit"), intercept = F, 
                             control=list(linear_program = "primal", purpose="find"))

# Standard errors for the two-step procedure
U11 <- dat11 %*% mA$coefficients[1:8]
U13 <- dat13 * mA$coefficients[9]
U14 <- dat14 * mA$coefficients[10]
normalize <- 1+p4hat^2 + (1-p4hat)^2
interior <- (U11-U13*(1-p4hat)-U14*p4hat)/sqrt(normalize)
chain <- (U13-U14)/sqrt(normalize) + ((1-2*p4hat)*((U11-U13*(1-p4hat)-U14*p4hat)))/(normalize^(3/2))
dLdP <- diag(drop(Y1*((dnorm(interior)/pnorm(interior))*( chain)) -
                    (1-Y1)*((dnorm(interior)/(pnorm(interior,lower=F)))*(chain))))
dLdB <-  drop(Y1*dnorm(interior)/pnorm(interior) - (1-Y1)*dnorm(interior)/pnorm(interior,lower=F)) * cbind(dat11,-dat13*(1-p4hat), -dat14*(p4hat))/sqrt(1+p4hat^2 + (1-p4hat)^2)
OmegaB.inv <- vcov(mA)
Omegap <- crossprod(dLdB, dLdP)
SIGMA <- drop(dnorm(as.matrix(dat2) %*% mB$coef)^2) *(as.matrix(dat2) %*% vcov(mB) %*% t(as.matrix(dat2)))
SSE2 <- c(sqrt(diag(OmegaB.inv + OmegaB.inv %*% Omegap %*% SIGMA %*% t(Omegap) %*% OmegaB.inv)), sqrt(diag(vcov(mB))))
estNC <- c(mA$coefficients, mB$coefficients)
NC <- list(mA, mB)

#### BR- SBI ####
Y1 <- dat[,1]
Y2 <- dat[,2][Y1==0]
dat11 <- dat[,names$U1]
dat14 <- 1
dat2 <- cbind(1,dat[,names$U2])/sqrt(2)
colnames(dat2)[1] <- "int"
dat13 <- 1
dat2part <-  dat2[Y1==0,]
mB <- brglm(Y2~dat2part-1, 
            x=TRUE, intercept=F,
            family=binomial(link="probit"),
            pl=TRUE)
p4hat <- as.numeric(pnorm(as.matrix(dat2) %*% mB$coef/sqrt(2)))
Z11hat <- p4hat *dat14
Z10hat <- (1-p4hat)*dat13
dat1 <- cbind(dat11, -Z10hat, -Z11hat)/sqrt(1+p4hat^2 + (1-p4hat)^2)
colnames(dat1) <- c(names$U1, "int3", "int4")
mA <- brglm(Y1~dat1-1,x=TRUE,
            family=binomial(link="probit"), pl=TRUE)

fullDesign <- cbind.data.frame(dat1, dat2); colnames(fullDesign) <- paste("V", 1:ncol(fullDesign), sep="")

# standard errors # 
U11 <- dat11 %*% mA$coefficients[1:8]
U13 <- dat13 * mA$coefficients[9]
U14 <- dat14 * mA$coefficients[10]
normalize <- 1+p4hat^2 + (1-p4hat)^2
interior <- (U11-U13*(1-p4hat)-U14*p4hat)/sqrt(normalize)
chain <- (U13-U14)/sqrt(normalize) + ((1-2*p4hat)*((U11-U13*(1-p4hat)-U14*p4hat)))/(normalize^(3/2))

dLdP <- diag(drop(Y1*((dnorm(interior)/pnorm(interior))*( chain)) -
                    (1-Y1)*((dnorm(interior)/(pnorm(interior,lower=F)))*(chain))))
dLdB <-  drop(Y1*dnorm(interior)/pnorm(interior) - (1-Y1)*dnorm(interior)/pnorm(interior,lower=F)) * cbind(dat11,-dat13*(1-p4hat), -dat14*(p4hat))/sqrt(1+p4hat^2 + (1-p4hat)^2)

OmegaB.inv <- vcov(mA)
Omegap <- crossprod(dLdB, dLdP)
SIGMA <- drop(dnorm(as.matrix(dat2) %*% mB$coef)^2) *(as.matrix(dat2) %*% vcov(mB) %*% t(as.matrix(dat2)))
SSE4 <- c(sqrt(diag(OmegaB.inv + OmegaB.inv %*% Omegap %*% SIGMA %*% t(Omegap) %*% OmegaB.inv)), sqrt(diag(vcov(mB))))
estWC <- c(mA$coefficients, mB$coefficients)
WC <- list(mA, mB)






##Post estimation FIML test#
Y1 <- dat[,1]
Y2 <- dat[,2][Y1==0]
dat11 <- dat[,names$U1]
dat14 <- 1
dat2 <- cbind(1,dat[,names$U2])
colnames(dat2)[1] <-"int"
dat13 <- 1
dat2part <-  dat2[Y1==0,]
p4hat <- predict(l1, type="action")[,4]
Z11hat <- p4hat *dat14
Z10hat <- (1-p4hat)*dat13
dat1 <- cbind(dat11, -Z10hat, -Z11hat)
colnames(dat1) <- c(names$U1, "int3", "int4")
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



tabD1 <- data.frame(Regressors = c("X[B]", "Z[SBI]", rep("Z[FIML] & X[B]", 3)),
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
cat("Table D.1 \n", file="../tables_and_figures/TableD1.md")
cat(kable(tabD1),
    sep="\n",
    file="../tables_and_figures/TableD1.md",
    append=TRUE)
#### BR FIML ####
l2 <- egame12(outcome ~ capcont + lreserves + overval + creditgrow +
                USinterest + service + contagion + prioratt - 1 | 1 | 1 | unifgov + lexports
              + preelec + postelec + rightgov + realinterest + capcont + lreserves, data =
                leblang2003, link = "probit", type = "private", penalty="logF",
              start=l1$coefficients)

#### Convert to table ####

Lcoefs <- list(l1$coef, estNC,
               l2$coef, estWC)
LSE  <-  list(sqrt(diag(l1$vcov)), SSE2,
              sqrt(diag(l2$vcov)), SSE4)

names <- c("Capital Controls",
           "Log(Reserves)",
           "Overvalued",
           "Credit Growth",
           "U.S. Interest",
           "Service",
           "Contagion",
           "Prior Attack",
           "Devaluation",
           "Defense",
           "Const.",
           "Unified Gov't",
           "Log(Exports)",
           "Pre-election",
           "Post-election",
           "Right Gov't",
           "Real Interest",
           "Capital Control",
           "Log(Reserves)")


names <- paste(c(rep("$U_A$(SQ):", 8),
                 "$U_A$(Devalue):",
                 "$U_A$(Defend):",
                 rep("$U_B$(Defend):", 9)), names) 


info.list <- list(UFMLE=list(logLik = logLik(l1),
                             Obs =  nrow(l1$model)),
                  USBI =list(logLik = "",
                             Obs =  nrow(l1$model)),
                  CFMLE=list(logLik = logLik(l2),
                             Obs =  nrow(l2$model)),
                  CSBI =list(logLik ="",
                             Obs =  nrow(l2$model)))



model.names <- c("FMLE",
                 "SBI",
                 "BR FMLE",
                 "BR SBI")

formattedCoef <- lapply(Lcoefs, round, 2)
formattedSE<- lapply(LSE, function(x){paste0("(", formatC(x,format ="f",digits=2),")")})
tabOut <- matrix("", nrow=length(Lcoefs[[1]])*2, ncol=4)
for(i in 1:4){tabOut[,i] <- matrix(rbind(formattedCoef[[i]], formattedSE[[i]]), ncol=1)}
tabOut <- cbind(matrix(rbind(names, ""), ncol=1), tabOut)
colnames(tabOut) <- c("", model.names)
tabOut <- rbind(tabOut, c("Obseravations", sapply(info.list, function(x)return(x$Obs))))
cat("Table D2\n", file="../tables_and_figures/TableD2.md")
cat(kable(tabOut),
    sep="\n",
    file="../tables_and_figures/TableD2.md",
    append = TRUE)
save(list=c("Lcoefs", "LSE", "info.list"), file="Leblang_output.rdata")
