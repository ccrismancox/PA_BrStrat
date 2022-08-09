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
library(xtable)
rm(list=ls())
load('table1.rdata')

k <- length(truth)
Index <- list(coefs=matrix(rep(1:(k*6)), 6, byrow=TRUE),
              separation=(k*6+1):(k*6+7),
              standardErrors=matrix((k*6+8):(k*6+7 + k*6), 6, byrow=TRUE),
              p = (k*6+7 + k*6+1):(k*6+7 + k*6+24),
              coverage=matrix((k*6+7 + k*6+25):(k*6+7 + k*6+25+k*6-1), 6, byrow=TRUE),
              power=matrix((k*6+7 + k*6+25+k*6):(k*6+7 + k*6+25+k*6+(k)*6-1), 6, byrow=TRUE),
              fixedp=matrix((k*6+7 + k*6+25+k*6+(k)*6):(k*6+7 + k*6+25+k*6+k*6+2*k*4-1),2, byrow=TRUE),
              trueSE=( k*6+7 + k*6+25+k*6+k*6+2*k*4):(k*6+7 + k*6+25+k*6+k*6+2*k*4+7))



table(out[,Index$separation[5]]) # Check all covariates against SQ 0
table(out[,Index$separation[6]]) # Check all covariates against BD 4762 (good)
table(out[,Index$separation[7]]) # Check all covariates against SF 205 (interesting)

# only keep the separation cases 
out <- out[out[,Index$separation[6]]==1,]

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



sdtrueSBI <- colMeans(out[,Index$trueSE[1:4]])
sdtrueFIML <- colMeans(out[,Index$trueSE[5:8]],na.rm=T)



output2 <- list(rmse=c(rmseSBI.nc, rmseSBI.wc, rmseFML.nc, rmseFML.wc,rmseFML.wc2,rmseFML.wc3),
                sd=sqrt(c(seSBI.nc, seSBI.wc, seFML.nc, seFML.wc,seFML.wc2,seFML.wc3)),
                est=c(estSBI.nc, estSBI.wc, estFML.nc, estFML.wc, estFML.wc2, estFML.wc3))
output2$est[seq(1,k*6,by=4)] <- output2$est[seq(1,k*6,by=4)]*-1 #convert from SBI to regular form 
output2$coverage <- c(colMeans(out[,Index$coverage[1,]], na.rm=T),
                      colMeans(out[,Index$coverage[2,]], na.rm=T),
                      colMeans(out[,Index$coverage[3,]], na.rm=T),
                      colMeans(out[,Index$coverage[4,]], na.rm=T),
                      colMeans(out[,Index$coverage[5,]], na.rm=T),
                      colMeans(out[,Index$coverage[6,]], na.rm=T))
output2$power <-  c(colMeans(out[,Index$power[1,]], na.rm=T),
                    colMeans(out[,Index$power[2,]], na.rm=T),
                    colMeans(out[,Index$power[3,]], na.rm=T),
                    colMeans(out[,Index$power[4,]], na.rm=T),
                    colMeans(out[,Index$power[5,]], na.rm=T),
                    colMeans(out[,Index$power[6,]], na.rm=T))
output2$SE <-  c(colMeans(out[,Index$standardErrors[1,]], na.rm=T),
                 colMeans(out[,Index$standardErrors[2,]], na.rm=T),
                 colMeans(out[,Index$standardErrors[3,]], na.rm=T),
                 colMeans(out[,Index$standardErrors[4,]], na.rm=T),
                 colMeans(out[,Index$standardErrors[5,]], na.rm=T),
                 colMeans(out[,Index$standardErrors[6,]], na.rm=T))
Table1 <- matrix(t(cbind(matrix(output2$est,nrow=6, byrow=T),
                         matrix(output2$sd, nrow=6, byrow=T),
                         matrix(output2$SE, nrow=6, byrow=T),
                         matrix(output2$power, nrow=6, byrow=T),
                         matrix(output2$coverage, nrow=6, byrow=T)
)), 
ncol=4, byrow=T)
Table1 <- round(Table1,2)
Table1 <- as.data.frame(cbind(Table1, c(rbind(output2$rmse, matrix(NA,4,6)))))
Table1 <- rbind(Table1, c(b0,NA), c(sdtrueSBI, NA), c(sdtrueFIML, NA))
colnames(Table1) <- c("$\\alpha_0$", "$\\alpha_1$", "$\\beta_0$", "$\\beta_1$", "RMSE")
Table1 <- cbind.data.frame(data.frame(Estimator=c("Ordinary SBI", rep(NA,4),
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
                           Table1)
options(knitr.kable.NA="")
cat("Table 1\n", file="../tables_and_figures/Table1.md")
cat(kable(Table1, digits=2), 
    file = "../tables_and_figures/Table1.md",
    append = TRUE)



############## what about pB (Table 2)#########

pB.analysis <- cbind(
  colMeans(out[,Index$p[1:6]], na.rm=TRUE),
  sqrt(colMeans(out[,Index$p[1:6]]^2, na.rm=TRUE)), 
  colMeans(out[,Index$p[7:12]], na.rm=TRUE),
  sqrt(colMeans(out[,Index$p[7:12]]^2, na.rm=TRUE)),
  matrix(colMeans(out[,Index$p[13:24]], na.rm=TRUE),nrow=6, byrow=T))
pB.analysis[,ncol(pB.analysis)] <- sqrt(pB.analysis[,ncol(pB.analysis)])


colnames(pB.analysis) <- c("Bias ($X_B=0$)", "RMSE ($X_B=0$)",
                           "Bias ($X_B=1$)", "RMSE ($X_B=1$)",
                           "Bias (total)", "RMSE (total)")
rownames(pB.analysis) <- c("Ordinary SBI",
                           "BR-SBI",
                           "Ordinary FIML",
                           "BR-FIML (Firth)",
                           "BR-FIML (Cauchy)",
                           "BR-FIML (log-$F$)")
pB.analysis <- pB.analysis
cat("Table 2\n", file="../tables_and_figures/Table2.md")
cat(kable(pB.analysis, digits=3), 
    file = "../tables_and_figures/Table2.md",
    append = TRUE)



############## what about a fixed pB (Table B.1)#########
fixed.p <-rbind.data.frame(colMeans(out[,Index$fixedp[1,]][,1:2]),
                           colMeans(out[,Index$fixedp[2,]][,1:2]),
                           #sampling sd
                           colSds(out[,Index$fixedp[1,]][,1:2], na.rm=T),
                           colSds(out[,Index$fixedp[2,]][,1:2], na.rm=T),
                           
                           #sd err
                           colMeans(out[,Index$fixedp[1,]][,5:6], na.rm=T),
                           colMeans(out[,Index$fixedp[2,]][,5:6], na.rm=T),
                           
                           #coverage
                           colMeans(out[,Index$fixedp[1,]][,9:10], na.rm=T),
                           colMeans(out[,Index$fixedp[2,]][,9:10], na.rm=T),
                           
                           #power
                           colMeans(out[,Index$fixedp[1,]][,13:14], na.rm=T),
                           colMeans(out[,Index$fixedp[2,]][,13:14], na.rm=T))
fixed.p <- fixed.p[c(seq(1, 10, by=2), seq(2,10,by=2)),]
fixed.p <- cbind.data.frame(c("Ordinary SBI",  rep(NA, 4),
                              "BR-SBI",rep(NA,4)),
                            rep(c("Expected value",
                                  "St.\ Dev.",  
                                  "St.\ Err.",
                                  "Coverage", 
                                  "Power"),2),
                            fixed.p)

colnames(fixed.p) <- c("Estimator", "Statistic", "$\\alpha_0$", "$\\alpha_1$")
fixed.p[c(1,6), 3:4] <- fixed.p[c(1,6), 3:4]*-1
cat("Table B.1\n", file="../tables_and_figures/TableB1.md")
cat(kable(fixed.p, digits=2), 
    file = "../tables_and_figures/TableB1.md",
    append = TRUE)


