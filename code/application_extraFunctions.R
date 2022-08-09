
# sepTest <- function(Y, X){
#   require(stringr)
#   catch <- tryCatch(glm.fit(x=X, y=Y, family=binomial()), warning=function(w){w})
#   if(class(catch)[1] == "list"){
#     return(FALSE)
#   }else{
#     if(str_detect(as.character(catch),"fitted probabilities numerically 0 or 1 occurred")){
#       return(TRUE)
#     }else{
#       cat("Unknown warning detected:", paste("'", catch, "'", sep=""))
#       return(FALSE)
#     }
#   }
# }

convTest <- function(Y,X, fn){
  require(stringr)
  catch <- tryCatch(fn(x=X, y=Y, family=binomial(link="probit")), warning=function(w){w})
  if(class(catch)[1] == "list"){
    return(catch$coef)
  }else{
    if(str_detect(as.character(catch),"fitted probabilities numerically 0 or 1 occurred")  | 
       str_detect(as.character(catch),"#successes in a binomial glm!") ){
      return(fn(x=X, y=Y, family=binomial(link="probit"))$coef)
    }else{
      cat("Unknown warning detected:", paste("'", catch, "'\n", sep=""))
      return(NA)
    }
  }
  
}
# 
# CorrectedLL <- function(b, X, Y){
#   XB <- X %*% b
#   p <- plogis(XB)
#   dW <- p * (1-p)
#   W <- diag(as.numeric(dW))
#   Ib <- t(X) %*% W %*% X
#   LL <- (1-Y) * (-XB) - log(1 + exp(-XB)) 
#   LL <- LL +
#     1/2 * log(det(Ib))/length(Y)
#   return(sum(LL))
# }
# 
# logitLL <- function(b, x, y){
#   xb=x%*%b
#   logD=ifelse(-xb>708,-xb,log(1+exp(-xb)))
#   llik=(1-y)*(-xb)-logD
#   sum(llik)		
# }
# 
# 
# mse <- function(est, true, se){
#   return(sum(se^2) + norm(est-true, "2")^2)
# }
# mse2 <- function(est, true, se){
#   bias2 <- (est-true)^2
#   var <- se^2
#   return(mean(var + bias2))
# }
# 
# SBI2 <- function(data, i, cols, A=0.2){
#   BS <- data[i,]
#   y1 <- BS[,1]
#   y2 <- BS[,2]
#   X11 <- BS[,3:(2+cols[1]), drop=FALSE]/sqrt(2)
#   X13 <- BS[,(2+cols[1]+1):(2+cols[1]+cols[2]), drop=FALSE]/sqrt(2)
#   X14 <- BS[,(2+cols[1]+cols[2]+1):(2+cols[1]+cols[2]+cols[3]), drop=FALSE]/sqrt(2)
#   X2 <- BS[,(2+cols[1]+cols[2]+cols[3]+1):(ncol(BS)), drop=FALSE]/sqrt(2)
#   
#   
#   mB <- glm.fit(y=y2[y1==1],x=(X2[y1==1,]), family=binomial())
#   p4hat <- plogis(X2 %*% A*mB$coef/sqrt(2))
#   p3hat <- 1-p4hat
#   Z11hat <- as.numeric(p4hat)*X14
#   Z10hat <- as.numeric(p3hat)*(X13)
#   mA <- glm.fit(y=y1,x=cbind(-X11, Z10hat, Z11hat), family=binomial())
#   na.omit(c(mA$coef, mB$coef))
# }
# SBI2C <- function(data, i, cols, A=0.2){
#   BS <- data[i,]
#   y1 <- BS[,1]
#   y2 <- BS[,2]
#   X11 <- BS[,3:(2+cols[1]), drop=FALSE]/sqrt(2)
#   X13 <- BS[,(2+cols[1]+1):(2+cols[1]+cols[2]), drop=FALSE]/sqrt(2)
#   X14 <- BS[,(2+cols[1]+cols[2]+1):(2+cols[1]+cols[2]+cols[3]), drop=FALSE]/sqrt(2)
#   X2 <- BS[,(2+cols[1]+cols[2]+cols[3]+1):(ncol(BS)), drop=FALSE]/sqrt(2)
#   
#   
#   mB <- brglm.fit(y=y2[y1==1],x=(X2[y1==1,]), family=binomial())
#   p4hat <- plogis(X2 %*% A*mB$coef/sqrt(2))
#   p3hat <- 1-p4hat
#   Z11hat <- as.numeric(p4hat)*X14
#   Z10hat <- as.numeric(p3hat)*(X13)
#   mA <- glm.fit(y=y1,x=cbind(-X11, Z10hat, Z11hat), family=binomial())
#   na.omit(c(mA$coef, mB$coef))
# }
# 
# SBI1C <- function(data, i, cols, A=0.2){
#   BS <- data[i,]
#   y1 <- BS[,1]
#   y2 <- BS[,2]
#   X11 <- BS[,3:(2+cols[1]), drop=FALSE]/sqrt(2)
#   X13 <- BS[,(2+cols[1]+1):(2+cols[1]+cols[2]), drop=FALSE]/sqrt(2)
#   X14 <- BS[,(2+cols[1]+cols[2]+1):(2+cols[1]+cols[2]+cols[3]), drop=FALSE]/sqrt(2)
#   X2 <- BS[,(2+cols[1]+cols[2]+cols[3]+1):(ncol(BS)), drop=FALSE]/sqrt(2)
#   
#   
#   mB <- glm.fit(y=y2[y1==1],x=(X2[y1==1,]), family=binomial())
#   p4hat <- plogis(X2 %*% A*mB$coef/sqrt(2))
#   p3hat <- 1-p4hat
#   Z11hat <- as.numeric(p4hat)*X14
#   Z10hat <- as.numeric(p3hat)*(X13)
#   mA <- brglm.fit(y=y1,x=cbind(-X11, Z10hat, Z11hat), family=binomial())
#   na.omit(c(mA$coef, mB$coef))
# }
# 
# 
# 
# 
# SBI12C <- function(data, i, cols, A=0.2){
#   BS <- data[i,]
#   y1 <- BS[,1]
#   y2 <- BS[,2]
#   X11 <- BS[,3:(2+cols[1]), drop=FALSE]/sqrt(2)
#   X13 <- BS[,(2+cols[1]+1):(2+cols[1]+cols[2]), drop=FALSE]/sqrt(2)
#   X14 <- BS[,(2+cols[1]+cols[2]+1):(2+cols[1]+cols[2]+cols[3]), drop=FALSE]/sqrt(2)
#   X2 <- BS[,(2+cols[1]+cols[2]+cols[3]+1):(ncol(BS)), drop=FALSE]/sqrt(2)
#   
#   
#   mB <- brglm.fit(y=y2[y1==1],x=(X2[y1==1,]), family=binomial())
#   p4hat <- plogis(X2 %*% A*mB$coef/sqrt(2))
#   p3hat <- 1-p4hat
#   Z11hat <- as.numeric(p4hat)*X14
#   Z10hat <- as.numeric(p3hat)*(X13)
#   mA <- brglm.fit(y=y1,x=cbind(-X11, Z10hat, Z11hat), family=binomial())
#   na.omit(c(mA$coef, mB$coef))
# }
# 
# 
# num2str <- function(x, digits=2){formatC(x, digits=digits, format='f')}
