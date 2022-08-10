# ###### If needed
library(devtools)
install_github("ccrismancox/games2")
# #####
library(games2)
library(doParallel)
library(doRNG)
library(brglm)
library(detectseparation)
library(matrixStats)
library(Formula)
library(knitr)
rm(list=ls())

source("extraFunctions.r")
truth <- c(2, .25, 2.5, -.5,.25, 5)
N <- 500
b0 <- truth
b0[1] <- -b0[1]

workers <- makeCluster(37) 
registerDoParallel(workers)
set.seed(1)
out <- foreach(i =  1:5000, .packages=c("games2", "brglm", "detectseparation"),
               .combine=rbind, .errorhandling = "remove")%dorng%{
                 
                 bigX <- mvrnorm(N, mu=rep(0,3), Sigma = matrix(c(1,.7,-.7,
                                                                  .7, 1,-.1,
                                                                  -.7,-.1,1), byrow=T,nrow=3))         
                 xB <- bigX[,2]
                 xA <- bigX[,3]
                 X <- ifelse(bigX[,1]>0, 1,0)
                 
                 cor(cbind(xB,xA,X))
                 
                 yB <- truth[4] + truth[5]*xB+  truth[6]* X+ (rnorm(N)-rnorm(N))
                 yB <- ifelse(yB>0, 1,0)
                 p4 <- pnorm( (truth[4] + truth[5]*xB+  truth[6]* X)/sqrt(2))
                 
                 yA <-  truth[1] + drop(truth[2:3] %*% t(cbind(xA,X)*p4))+ (rnorm(N)-rnorm(N))
                 yA <- ifelse(yA>0, 1,0)
                 
                 Asep=checksep(X*p4,yA)
                 Bsep=checksep(X[yA==1],yB[yA==1])
                 
                 Bsep.lp <- detect_separation(x = data.frame(int=1,
                                                             X1=xB[yA==1], 
                                                             X2=X[yA==1]),
                                              y=yB[yA==1], 
                                              family=binomial("probit"), intercept = F)$separation
                 Asep.lp <- detect_separation(x=data.frame(int=1,X1=xA*p4,X2=X*p4),
                                              yA, 
                                              family=binomial("probit"), intercept = F)$separation
                 total.lp0  <- detect_separation(x=data.frame(int=1,
                                                              X1=xB, 
                                                              X2=X,
                                                              z1=xA*p4,
                                                              z2=X*p4), 
                                                 y=(yA + (yA*yB))==0, family=binomial("probit"),
                                                 intercept = F)$separation
                 total.lp1  <- detect_separation(x=data.frame(int=1,
                                                              X1=xB, 
                                                              X2=X,
                                                              z1=xA*p4, 
                                                              z2=X*p4), 
                                                 y=(yA + (yA*yB))==1, 
                                                 family=binomial("probit"), intercept = F)$separation
                 total.lp2  <- detect_separation(x=data.frame(int=1,
                                                              X1=xB, 
                                                              X2=X,
                                                              z1=xA*p4, 
                                                              z2=X*p4),
                                                 y=(yA + (yA*yB))==2,
                                                 family=binomial("probit"), intercept = F)$separation
                 
                 mf = data.frame(yA,yB,xA,xB,X)
                 
                 mf.sbi <- mf
                 mf.sbi$xB <- mf$xB/sqrt(2)
                 mf.sbi$Xroot2 <- mf$X/sqrt(2)
                 mf.sbi$const <- 1/sqrt(2)
                 b0 <- truth
                 b0[1] <- -b0[1]
                 trueSBI.se <- standardErrors.sbi(yA + yB ~ 1 | 0 | xA+X-1 | xB+X, mf =mf ,
                                                  par=truth)$SE
                 trueFIML.se <- standardErrors.fiml(yA + yB ~ 1 | 0 | xA+X-1 | xB+X, mf = mf,
                                                    par=b0)
                 
                 step1WC <- try(brglm(yB~const+xB-1,subset=yA==1, data=mf.sbi, family=binomial(link="probit"), pl=TRUE))
                 
                 if(class(step1WC)[1]=="try-error"){
                   step1WC <- step2WC <- list(coef=rep(NA,3))
                 }else{
                   p4hat <- pnorm(with(mf.sbi, cbind(const, xB)) %*% step1WC$coef)
                   mf.sbi$Zhat1 <- xA*p4hat/sqrt(2)
                   mf.sbi$Zhat2 <- X*p4hat/sqrt(2)
                   step2WC <-  try(brglm(yA~const + Zhat1-1, data=mf.sbi, family=binomial(link="probit"), pl=TRUE))

                   if(class(step2WC)[1]=="try-error"){ 
                     step2WC <- list(coef=rep(NA,3))
                     SBIWC <- c(step2WC$coef,0,step1WC$coef,0)
                     SBIWC.se <- rep(NA,6)
                     pB.WC.bias.SBI <-NA
                     pB.WC.bias.SBI.xB <-NA
                     coverage.WC.SBI <-rep(NA,6)
                     power.WC.SBI <-  rep(NA,6)
                   }else{
                     SBIWC <- c(step2WC$coef,0,step1WC$coef,0)
                     
                     SBIWC.se <- standardErrors.sbi(yA + yB ~ 1 | 0 | xA-1 | xB,
                                                    mf = data.frame(yA,yB,xA,xB,X),
                                                    par=c(step2WC$coef,step1WC$coef))$SE
                     SBIWC.se <- c(SBIWC.se[1:2], NA,
                                   SBIWC.se[3:4], NA)
                     pB.WC.bias.SBI <- mean(p4hat-p4)
                     pB.WC.bias.SBI.xB <- mean(p4hat[xB==1]-p4[xB==1])
                     
                     coverage.WC.SBI <- ( truth < SBIWC + 1.96*SBIWC.se & truth > SBIWC - 1.96*SBIWC.se )
                     power.WC.SBI <-  1-( 0 < SBIWC + 1.96*SBIWC.se & 0 > SBIWC - 1.96*SBIWC.se )
                     coverage.WC.SBI[c(3,6)] <- power.WC.SBI[c(3,6)] <-NA
                   }
                 }    
                 
                 step1NC <- try(glm(yB~const+xB-1,subset=yA==1, data=mf.sbi, family=binomial(link="probit")))
                 if(any(class(step1NC)=="try-error")|| step1NC$converged==FALSE||anyNA(step1NC$coef)){
                   step1NC <- step2NC <- list(coef=rep(NA,3))
                 }else{
                   p4hat <- pnorm(with(mf.sbi, cbind(const, xB)) %*% step1NC$coef)
                   mf.sbi$Zhat1 <- xA*p4hat/sqrt(2)
                   mf.sbi$Zhat2 <- X*p4hat/sqrt(2)
                   
                 
                   step2NC <- try(glm(yA~const+Zhat1-1, data=mf.sbi, family=binomial(link="probit")))
                   if(any(class(step2NC)=="try-error") || step2NC$converged==FALSE){
                     step2NC <- list(coef=rep(NA,3))
                     SBINC <- c(step2NC$coef,step1NC$coef)
                     
                     SBINC.se <-rep(NA,6)
                     pB.NC.bias.SBI <-pB.NC.bias.SBI.xB <- NA
                     coverage.NC.SBI <- rep(NA,6)
                     power.NC.SBI <- rep(NA,6)
                     
                     
                   }else{
                     SBINC <- c(step2NC$coef,0,
                                step1NC$coef,0)
                     
                     SBINC.se <- standardErrors.sbi(yA + yB ~ 1 | 0 | xA-1 | xB,
                                                    mf = data.frame(yA,yB,xA,xB,X),
                                                    par=c(step2NC$coef,step1NC$coef))$SE
                     SBINC.se <- c(SBINC.se[1:2], NA,
                                   SBINC.se[3:4], NA)
                     pB.NC.bias.SBI <- mean(p4hat-p4)
                     pB.NC.bias.SBI.xB<- mean(p4hat[xB==1]-p4[xB==1])
                     
                     coverage.NC.SBI <- ( truth < SBINC + 1.96*SBINC.se & truth > SBINC - 1.96*SBINC.se )
                     power.NC.SBI <- 1- ( 0 < SBINC + 1.96*SBINC.se & 0 > SBINC - 1.96*SBINC.se )
                     coverage.NC.SBI[c(3,6)] <- power.NC.SBI[c(3,6)] <-NA
                     
                   }
                 }    

                 
                 
                 NC <- try(egame12(yA + yB ~ 1 | 0 | xA-1 | xB, start='zero', link="probit"))
                 if(any(class(NC)=="try-error" )|| NC$convergence$code != 0||NC$localID==FALSE){
                   NC <- NC.se <- NC.coverage <- NC.power <- rep(NA,6)
                   pB.NC.bias <- pB.NC.bias.xB <- NA
                 }else{
                   NC.se <- sqrt(diag(vcov(NC)))
                   NC.se <- c(NC.se[1:2], NA,
                              NC.se[3:4],NA)
                   phat <- pnorm(with(mf, cbind(1,xB) %*% NC$coef[3:4]/sqrt(2)))
                   pB.NC.bias <- mean(phat-p4)
                   pB.NC.bias.xB <-  mean(phat[xB==1]-p4[xB==1])
                   NC <- c(NC$coef[1:2] , 0,
                           NC$coef[3:4], 0)
                   coverage.NC <- ( b0 < NC + 1.96*NC.se & b0 > NC - 1.96*NC.se )
                   
                   power.NC <-  1-( 0 < NC + 1.96*NC.se & 0 > NC - 1.96*NC.se )
                   NC[1] <- -NC[1]
                   coverage.NC[c(3,6)] <- power.NC[c(3,6)] <-NA
                 }
                 
                 
                 WC <- try(egame12(yA + yB ~ 1 | 0 | xA-1 | xB, start='zero', link="probit"))
                 if(any(class(WC)=="try-error" )|| WC$convergence$code != 0||WC$localID==FALSE){
                   WC <- WC.se <- WC.coverage <- WC.power <- rep(NA,6)
                   pB.WC.bias <- pB.WC.bias.xB <- NA
                 }else{
                   WC.se <- sqrt(diag(vcov(WC)))
                   WC.se <- c(WC.se[1:2], NA,
                              WC.se[3:4], NA)
                   phat <- pnorm(with(mf, cbind(1,xB) %*% WC$coef[3:4]/sqrt(2)))
                   pB.WC.bias <- mean(phat-p4)
                   pB.WC.bias.xB <-  mean(phat[xB==1]-p4[xB==1])
                   WC <- c(WC$coef[1:2] ,0,
                           WC$coef[3:4],0)
                   coverage.WC <- ( b0 < WC + 1.96*WC.se & b0 > WC - 1.96*WC.se )
                   power.WC <-  1-( 0 < WC + 1.96*WC.se & 0 > WC - 1.96*WC.se )
                   coverage.WC[c(3,6)] <- power.WC[c(3,6)] <-NA
                   
                   WC[1] <- -WC[1]
                 }
                 
                 
                 WC2 <- try(egame12(yA + yB ~ 1 | 0 | xA-1 | xB, start='zero', link="probit"))
                 if(any(class(WC2)=="try-error" )|| WC2$convergence$code != 0||WC2$localID==FALSE){
                   WC2 <- WC2.se <- WC2.coverage <- WC2.power <- rep(NA,6)
                   pB.WC2.bias <- pB.WC2.bias.xB <- NA
                 }else{
                   WC2.se <- sqrt(diag(vcov(WC2)))
                   WC2.se <- c(WC2.se[1:2], NA,
                               WC2.se[3:4], NA)
                   phat <- pnorm(with(mf, cbind(1,xB) %*% WC2$coef[3:4]/sqrt(2)))
                   pB.WC2.bias <- mean(phat-p4)
                   pB.WC2.bias.xB <-  mean(phat[xB==1]-p4[xB==1])
                   WC2 <- c(WC2$coef[1:2] ,0,
                            WC2$coef[3:4],0)
                   coverage.WC2 <- ( b0 < WC2 + 1.96*WC2.se & b0 > WC2 - 1.96*WC2.se )
                   power.WC2 <-  1-( 0 < WC2 + 1.96*WC2.se & 0 > WC2 - 1.96*WC2.se )
                   coverage.WC2[c(3,6)] <- power.WC2[c(3,6)] <-NA
                   
                   WC2[1] <- -WC2[1]
                 } 
                 
                 WC3 <- try(egame12(yA + yB ~ 1 | 0 | xA-1 | xB, start='zero', link="probit"))
                 if(any(class(WC3)=="try-error" )|| WC3$convergence$code != 0||WC3$localID==FALSE){
                   WC3 <- WC3.se <- WC3.coverage <- WC3.power <- rep(NA,6)
                   pB.WC3.bias <- pB.WC3.bias.xB <- NA
                 }else{
                   WC3.se <- sqrt(diag(vcov(WC3)))
                   WC3.se <- c(WC3.se[1:2], NA,
                               WC3.se[3:4], NA)
                   phat <- pnorm(with(mf, cbind(1,xB) %*% WC3$coef[3:4]/sqrt(2)))
                   pB.WC3.bias <- mean(phat-p4)
                   pB.WC3.bias.xB <-  mean(phat[xB==1]-p4[xB==1])
                   WC3 <- c(WC3$coef[1:2] ,0,
                            WC3$coef[3:4],0)
                   coverage.WC3 <- ( b0 < WC3 + 1.96*WC3.se & b0 > WC3 - 1.96*WC3.se )
                   
                   power.WC3 <-  1-( 0 < WC3 + 1.96*WC3.se & 0 > WC3 - 1.96*WC3.se )
                   coverage.WC3[c(3,6)] <- power.WC3[c(3,6)] <-NA
                   
                   WC3[1] <- -WC3[1]
                 }
                 
        
                 output <- c(SBINC, 
                             SBIWC,  
                             NC,
                             WC,  
                             WC2, 
                             WC3,  
                             Asep, Bsep, 
                             Asep.lp, Bsep.lp,
                             total.lp0,total.lp1,total.lp2, 
                             SBINC.se, 
                             SBIWC.se, 
                             NC.se, 
                             WC.se, 
                             WC2.se, 
                             WC3.se, 
                             pB.NC.bias.SBI, 
                             pB.WC.bias.SBI, 
                             pB.NC.bias, 
                             pB.WC.bias, 
                             pB.WC2.bias,
                             pB.WC3.bias,
                             pB.NC.bias.SBI.xB, 
                             pB.WC.bias.SBI.xB, 
                             pB.NC.bias.xB, 
                             pB.WC.bias.xB, 
                             pB.WC2.bias.xB,
                             pB.WC3.bias.xB,
                             coverage.NC.SBI,
                             coverage.WC.SBI, 
                             coverage.NC, 
                             coverage.WC, 
                             coverage.WC2, 
                             coverage.WC3, 
                             power.NC.SBI, 
                             power.WC.SBI, 
                             power.NC,
                             power.WC, 
                             power.WC2, 
                             power.WC3, 
                             trueSBI.se, 
                             trueFIML.se 
                 )
                 
                 output
               }
stopCluster(workers)
save.image("tableB6.rdata")



load("tableB6.rdata")
k <- length(truth)
Index <- list(coefs=matrix(rep(1:(k*6)), 6, byrow=TRUE),
              separation=(k*6+1):(k*6+7),
              standardErrors=matrix((k*6+8):(k*6+7 + k*6), 6, byrow=TRUE),
              p = (k*6+7 + k*6+1):(k*6+7 + k*6+12),
              coverage=matrix((k*6+7 + k*6+13):(k*6+7 + k*6+13+k*6-1), 6, byrow=TRUE),
              power=matrix((k*6+7 + k*6+13+k*6):(k*6+7 + k*6+13+k*6+(k)*6-1), 6, byrow=TRUE),
              trueSE=( k*6+7 + k*6+13+k*6+k*6):(k*6+7 + k*6+13+k*6+k*6+2*k-1))



out <- out[out[,Index$separation[6]]==1 & (out[,Index$separation[5]]==1 | out[,Index$separation[4]]==1),]

biasSBI.nc <- rowMeans((t(out)[Index$coefs[1,], ] - truth ), na.rm=T)
biasSBI.wc <- rowMeans((t(out)[Index$coefs[2,], ] - truth ), na.rm=T)
biasFML.nc <- rowMeans((t(out)[Index$coefs[3,], ] - truth ), na.rm=T)
biasFML.wc <- rowMeans((t(out)[Index$coefs[4,], ] - truth ), na.rm=T)
biasFML.wc2 <- rowMeans((t(out)[Index$coefs[5,], ] - truth ), na.rm=T)
biasFML.wc3 <- rowMeans((t(out)[Index$coefs[6,], ] - truth ), na.rm=T)

estSBI.nc <- colMeans(out[,Index$coefs[1,]], na.rm=T)
estSBI.wc <- colMeans(out[,Index$coefs[2,]], na.rm=T)
estFML.nc <- colMeans(out[,Index$coefs[3,]], na.rm=T)
estFML.wc <- colMeans(out[,Index$coefs[4,]], na.rm=T)
estFML.wc2 <- colMeans(out[,Index$coefs[5,]], na.rm=T)
estFML.wc3 <- colMeans(out[,Index$coefs[6,]], na.rm=T)


seSBI.nc <- colVars(out[,Index$coefs[1,]], na.rm=T)
seSBI.wc <- colVars(out[,Index$coefs[2,]], na.rm=T)
seFML.nc <- colVars(out[,Index$coefs[3,]], na.rm=T)
seFML.wc <- colVars(out[,Index$coefs[4,]], na.rm=T)
seFML.wc2 <- colVars(out[,Index$coefs[5,]], na.rm=T)
seFML.wc3 <- colVars(out[,Index$coefs[6,]], na.rm=T)

rmseSBI.nc <- rmse(biasSBI.nc,seSBI.nc)
rmseSBI.wc <- rmse(biasSBI.wc,seSBI.wc)
rmseFML.nc <- rmse(biasFML.nc,seFML.nc)
rmseFML.wc <- rmse(biasFML.wc,seFML.wc)
rmseFML.wc2 <- rmse(biasFML.wc2,seFML.wc2)
rmseFML.wc3 <- rmse(biasFML.wc3,seFML.wc3)



sdtrueSBI <- colMeans(out[,Index$trueSE[1:6]])
sdtrueFIML <- colMeans(out[,Index$trueSE[7:12]],na.rm=T)



output <- list(rmse=c(rmseSBI.nc, rmseSBI.wc, rmseFML.nc, rmseFML.wc,rmseFML.wc2,rmseFML.wc3),
                sd=sqrt(c(seSBI.nc, seSBI.wc, seFML.nc, seFML.wc,seFML.wc2,seFML.wc3)),
                est=c(estSBI.nc, estSBI.wc, estFML.nc, estFML.wc, estFML.wc2, estFML.wc3))
output$est[seq(1,k*6,by=6)] <- output$est[seq(1,k*6,by=6)]*-1
output$coverage <- c(colMeans(out[,Index$coverage[1,]], na.rm=T),
                      colMeans(out[,Index$coverage[2,]], na.rm=T),
                      colMeans(out[,Index$coverage[3,]], na.rm=T),
                      colMeans(out[,Index$coverage[4,]], na.rm=T),
                      colMeans(out[,Index$coverage[5,]], na.rm=T),
                      colMeans(out[,Index$coverage[6,]], na.rm=T))
output$power <-  c(colMeans(out[,Index$power[1,]], na.rm=T),
                    colMeans(out[,Index$power[2,]], na.rm=T),
                    colMeans(out[,Index$power[3,]], na.rm=T),
                    colMeans(out[,Index$power[4,]], na.rm=T),
                    colMeans(out[,Index$power[5,]], na.rm=T),
                    colMeans(out[,Index$power[6,]], na.rm=T))
output$SE <-  c(colMeans(out[,Index$standardErrors[1,]], na.rm=T),
                 colMeans(out[,Index$standardErrors[2,]], na.rm=T),
                 colMeans(out[,Index$standardErrors[3,]], na.rm=T),
                 colMeans(out[,Index$standardErrors[4,]], na.rm=T),
                 colMeans(out[,Index$standardErrors[5,]], na.rm=T),
                 colMeans(out[,Index$standardErrors[6,]], na.rm=T))
TableB6 <- matrix(t(cbind(matrix(output$est,nrow=6, byrow=T),
                         matrix(output$sd, nrow=6, byrow=T),
                         matrix(output$SE, nrow=6, byrow=T),
                         matrix(output$power, nrow=6, byrow=T),
                         matrix(output$coverage, nrow=6, byrow=T)
)), 
ncol=6, byrow=T)
TableB6 <- round(TableB6,2)
TableB6 <- as.data.frame(cbind(TableB6, c(rbind(output$rmse, matrix(NA,4,6)))))
TableB6 <- rbind(TableB6, c(b0,NA), c(sdtrueSBI, NA), c(sdtrueFIML, NA))
colnames(TableB6) <- c("$\\alpha_0$", "$\\alpha_1$","$\\alpha_2$",
                      "$\\beta_0$", "$\\beta_1$", "$\\beta_2$", "RMSE")
TableB6 <- cbind.data.frame(data.frame(Estimator=c("Ordinary SBI", rep(NA,4),
                                                  "BR-SBI", rep(NA,4),
                                                  "Ordinary FIML", rep(NA,4),
                                                  "BR-FIML (Firth)", rep(NA,4),
                                                  "BR-FIML (Cauchy)", rep(NA,4),
                                                  "BR-FIML (log-$F$)", rep(NA,4),
                                                  "Truth",rep(NA,2)),
                                      Quantity=c( rep(c("Est","St.Dev.", "St.Err.", "Power", "Coverage"), 6),
                                                  "Parameters", 
                                                  "St. Err. (SBI)",
                                                  "St. Err. (FIML)")),
                           TableB6)

cat("Table B.6\n",    file = "../tables_and_figures/TableB6.md")
cat(kable(TableB6, digits=2), 
    sep="\n",
    file = "../tables_and_figures/TableB6.md",
    append = TRUE)

