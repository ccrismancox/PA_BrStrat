checksep=function(x,y){
  x0=x[y==0]
  x1=x[y==1]
  minx0=min(x0)
  maxx0=max(x0)
  minx1=min(x1)
  maxx1=max(x1)
  sep=(minx0 >= minx1 & minx0 >= maxx1) | (minx1 >= minx0 & minx1 >=maxx0)
  return(sep)
}

probitHessian <- function(theta, X, y){
  mu <- drop(X %*% theta)
  d <- dnorm(mu)
  p <- pnorm(mu)
  p0 <- pnorm(mu, lower=FALSE)
  # scalar <- -d * ((y*((d+mu*p)/(p^2))) +(1-y)*((d-mu*(p0))/(p0^2)))
  scalar <- -(d^2)/(p*p0)
  Xstar <- X*scalar
  H <- t(Xstar) %*% X
  return(H)
}

probit.br <- function(theta, X, y, penalty=c("Cauchy", "logF")){
  penalty=match.arg(penalty)
  mu <- drop(X %*% theta)
  ll <- ifelse(y==1,
               pnorm(mu,log=TRUE),
               pnorm(mu, lower=FALSE, log=TRUE))
  if(penalty=="Cauchy"){
    scale <- ifelse(str_detect(names(theta), "Intercept"), 10,2.5)
    ll <- sum(ll,sum(-log((1+ (theta/scale)^2))))
  }
  if(penalty=="logF"){
    ll <- sum(ll,sum(theta/2 -log(1+exp(theta))))
    
  }
  
  return(-sum(ll))
}

standardErrors.sbi <- function(formulas,mf,par){
  
  regr <- list()
  formulas <- as.Formula(formulas)
  for (i in seq_len(length(formulas)[2])){
    regr[[i]] <- model.matrix(formulas, data =mf, rhs = i)
  }
  rcols <- sapply(regr, ncol)
  
  b0 <- par
  alpha <- b0[1:sum(rcols[1:3])]
  beta <-  b0[(sum(rcols[1:3])+1):sum(rcols)]
  
  
  yf <- model.part(formulas, mf, lhs = 1, drop = TRUE)
  yf <- games2:::makeResponse12(yf)
  y <- as.numeric(yf)
  yB <- ifelse(y==1,NA, ifelse(y==2, 0,1))
  yA <- ifelse(y==1,0, 1)
  datB <- na.omit(cbind.data.frame(yB,regr[[4]]/sqrt(2)))
  datA <- cbind.data.frame(yA, do.call(cbind,regr[1:3])/sqrt(2))
  XB <- as.matrix(datB[,-1])
  XA <- as.matrix(datA[,-1])
  VB <- solve(-probitHessian(beta, X=XB, y=datB$yB))
  se.beta <- sqrt(diag(VB))
  
  pB <- pnorm(drop(regr[[4]] %*% beta/sqrt(2)))
  Z <- cbind(XA[,1:rcols[1]], #empty alpha_bd is hard coded
             XA[,(rcols[1]+1):sum(rcols[1:3])]*pB) 
  OmegaA.inv <- solve(-probitHessian(alpha, X=Z, y=datA$yA))
  
  dat11 <- Z[,1:rcols[1],drop=F]
  dat13 <- matrix(nrow=nrow(dat11), ncol=0)
  dat14 <-  Z[,(rcols[1]+1):sum(rcols[1:3]), drop=FALSE]
  dat2 <- regr[[4]]/sqrt(2)
  U11 <- dat11 %*% alpha[1:rcols[1]]
  U13 <- 0
  U14 <- dat14 %*% alpha[(rcols[1]+1):sum(rcols[1:3])]
  interior <- (U11+U13*(1-pB)+U14*pB)
  
  dLdP <- diag(drop(yA*((dnorm(interior)/(pnorm(interior)))*( (U14-U13))) -
                      (1-yA)*((dnorm(interior)/(pnorm(interior,lower=F)))*((U14-U13)))))
  dLdB <-  drop(yA*dnorm(interior)/(pnorm(interior)) - (1-yA)*dnorm(interior)/(pnorm(interior,lower=F))) *
    cbind(dat11,dat13*(1-pB), dat14*(pB))
  
  
  Omegap <- crossprod(dLdB, dLdP)
  SIGMA <- (drop(dnorm(dat2 %*% beta))*dat2) %*%VB %*% t((drop(dnorm(dat2 %*% beta))*dat2))
  VA <- OmegaA.inv + OmegaA.inv %*% Omegap %*% SIGMA %*% t(Omegap) %*% OmegaA.inv
  se.2step <- c(sqrt(diag(VA)),se.beta)
  return(list(SE=se.2step, VB=VB, VA=VA))
}
standardErrors.fiml <-  function(formulas, mf, par, link="probit", type="agent"){
  b0 <- par
  
  regr <- list()
  formulas <- as.Formula(formulas)
  for (i in seq_len(length(formulas)[2])){
    regr[[i]] <- model.matrix(formulas, data =mf, rhs = i)
  }
  rcols <- sapply(regr, ncol)
  yf <- model.part(formulas, mf, lhs = 1, drop = TRUE)
  yf <- games2:::makeResponse12(yf)
  y <- as.numeric(yf)
  ## generate true standard errors fiml
  u <-  games2:::makeUtils(b0, regr, nutils = 4,
                           unames = c("u11", "u13", "u14", "u24"))
  
  names(regr) <- character(length(regr))
  names(regr)[1:4] <- c("X1", "X3", "X4", "Z")
  rcols <- sapply(regr, ncol)
  n <- nrow(regr$Z)
  
  FirthExtra <- list()
  FirthExtra$regr2 <- regr
  FirthExtra$regr2$Z <- cbind(matrix(0,
                                     nrow=n, ncol= sum(rcols[1:3])),
                              regr$Z)
  FirthExtra$regr2$X1 <- cbind(regr$X1,
                               matrix(0,
                                      nrow=n, ncol= sum(rcols[2:4])))
  FirthExtra$regr2$X3 <- cbind(matrix(0,
                                      nrow=n, ncol= sum(rcols[1])),
                               regr$X3,
                               matrix(0,
                                      nrow=n, ncol= sum(rcols[3:4])))
  FirthExtra$regr2$X4 <- cbind(matrix(0,
                                      nrow=n, ncol= sum(rcols[1:2])),
                               regr$X4,
                               matrix(0,
                                      nrow=n, ncol= sum(rcols[4])))
  probs <- games2:::makeProbs12(b0, regr,link=link, type=type)
  
  A <- -games2:::hessian12(b0, y=y, regr=regr,link=link, type=type, FirthExtra=FirthExtra, p=probs, u=u)
  se.fiml.true <- sqrt(diag(solve(A)))
  return(se.fiml.true)
}

br.sbi <- function(formulas,mf,par, penalty=c("Firth", "Cauchy", "logF")){
  penalty <- match.arg(penalty)
  g <- switch(penalty,
              "Firth"= function(theta){.5*determinant(probitHessian(theta ,X=X, y=y))$modulus},
              "Cauchy"= function(theta){sum(-log((1+ (theta/2.5)^2)))},
              "logF"=function(theta){sum(b/2 -log(1+exp(b)))})
  
  regr <- list()
  formulas <- as.Formula(formulas)
  for (i in seq_len(length(formulas)[2])){
    regr[[i]] <- model.matrix(formulas, data =mf, rhs = i)
  }
  rcols <- sapply(regr, ncol)
  
  b0 <- par
  alpha <- b0[1:sum(rcols[1:3])]
  beta <-  b0[(sum(rcols[1:3])+1):sum(rcols)]
  
  
  yf <- model.part(formulas, mf, lhs = 1, drop = TRUE)
  yf <- games2:::makeResponse12(yf)
  y <- as.numeric(yf)
  yB <- ifelse(y==1,NA, ifelse(y==2, 0,1))
  yA <- ifelse(y==1,0, 1)
  datB <- na.omit(cbind.data.frame(yB,regr[[4]]/sqrt(2)))
  datA <- cbind.data.frame(yA, do.call(cbind,regr[1:3])/sqrt(2))
  XB <- as.matrix(datB[,-1])
  XA <- as.matrix(datA[,-1])
  
  brprobit.LL <- function(theta, X,y,g){
    ll <- sum(ifelse(y==1,
                     pnorm(X%*%theta, log.p=TRUE), 
                     pnorm(X%*%theta,log.p=TRUE, lower=FALSE)))+
      g(theta)
    return(ll)
  }
  mB <- optim(beta, brprobit.LL, X=XB, y=datB$yB)
  beta.hat <- mB$par
  pB.hat <- pnorm(regr[[4]]%*%beta.hat/sqrt(2))
  
  Z <- cbind(XA[,1:rcols[1]], #empty alpha_bd is hard coded
             XA[,(rcols[1]+1):sum(rcols[1:3])]*pB.hat) 
  mA <- optim(alpha, brprobit.LL, X=Z, y=yA)
  alpha.hat <- mA$par
  return(list(coef=c(alpha.hat,beta.hat), p4hat=pB.hat))
}


rmse <- function(bias, var){
  return(sqrt(sum(var)+ crossprod(bias)))
}
