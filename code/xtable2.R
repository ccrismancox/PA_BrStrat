# ======================================
#      Project: Tabling functions
#      Casey Crisman-Cox
#      Useful Functions
#      Tabling.R
# ======================================
#######Dependencies#####
if("xtable" %in% installed.packages()){
  library(xtable)
}else{
  install.packages("xtable")
  library(xtable)
}
if("stringr" %in% installed.packages()){
  library(stringr)
}else{
  install.packages("stringr")
  library(stringr)
}


num2str <- function(x, digits=2){formatC(x, digits=digits, format='f')}

mod.sum<-function(mod.list, 
                  apsr=FALSE, 
                  stars=c("default", "all", "CI", "none"),
                  seNote="Standard Errors in Parenthesis",
                  bootCI.level=.95,
                  model.numbers=NULL,
                  model.names=NULL,
                  dep.varnames=NULL,
                  k=1,
                  length=NULL){  ##Function to add in LogLik and N to xtable
  if(is.null(model.numbers)){
    model.numbers <- paste("\\multicolumn{1}{c}{Model ",
                           k:(k+length(mod.list)-1),
                           "}",
                           sep="")
  }
  if(is.null(dep.varnames)){
    dep.varnames <- rep("%", length(mod.list))
  }
  if(is.null(model.names)){
    model.names <- rep("%", length(mod.list))
  }
  
  out<-NULL
  if(apsr==FALSE){
    mod<-mod.list[[1]]
    ##Works on maxLik, lm, glm, game, and surv
    if("maxLik" %in% class(mod)){ ##Check if the class is maxLik
      out<-paste("\\hline ", 
                 "Log $L$ & \\multicolumn{4}{c}{",num2str(logLik(mod)), "}\\\\\n",
                 "$N$ & \\multicolumn{4}{c}{", nrow(mod$gradientObs), "}\\\\\n", 
                 sep="") ##Creates two rows which we'll add to xtable
    }
    if("lm" == class(mod)[1]){
      out<-paste("\\hline ",
                 "adj. R$^2$& \\multicolumn{4}{c}{", num2str(summary(mod)$adj.r.squared),
                 "}\\\\\n",
                 "$N$ & \\multicolumn{4}{c}{", nrow(mod$model), "}\\\\\n",
                 sep="") ##Same as above
    }
    
    if("glm" %in% class(mod) | "game" %in% class(mod) 
       |"polr" %in% class(mod)){
      out<-paste("\\hline ",
                 "Log $L$ & \\multicolumn{4}{c}{",num2str(logLik(mod)), "}\\\\\n",
                 "$N$ & \\multicolumn{4}{c}{", nrow(mod$model), "}\\\\\n", 
                 sep="") ##Creates two rows which we'll add to xtable
    }
    if("polywog" == class(mod)[1]){
      out<-paste("\\hline ",
                 "Lambda & \\multicolumn{4}{c}{",num2str(mod$lambda), "}\\\\\n",
                 "$N$ & \\multicolumn{4}{c}{", mod$nobs, "}\\\\\n", 
                 sep="") ##Creates two rows which we'll add to xtable
    }
    if("list" == class(mod)[1]){
      out<-paste("\\hline ",
                 "Log $L$ & \\multicolumn{4}{c}{",num2str(mod$logLik), "}\\\\\n",
                 "$N$ & \\multicolumn{4}{c}{", num2str(mod$Obs), "}\\\\\n", 
                 sep="") ##Creates two rows which we'll add to xtable
    }
    out<-list(out=out, 
              length=length(unique(names(coef(mod)))))
  }else{
    out.L<-c("\\hline Log $L$")
    out.A <- c("$\\alpha$")
    out.Lam<-c("$\\lambda$")
    out.Rho<-c("$\\rho$")
    out.R<-c("adj. R$^2$")
    out.N<-c("$N$")
    modNums <- " "
    if(all(dep.varnames=="%")){ 
      depNames <- "%"
    }else{
      depNames <- " "
    }
    if(all(model.names=="%")){ 
      modNames <- "%"
    }else{
      modNames <- " "
    }
    
    
    if(length(dep.varnames)< length(mod.list)){
      dep.varnames <- c(dep.varnames, rep("%", length(mod.list)-  length(dep.varnames)))
    }
    if(length(model.names)< length(mod.list)){
      model.names <- c(model.names, rep("%", length(mod.list)-  length(model.names)))
    }
    
    
    for(i in 1:length(mod.list)){
      mod.i<-mod.list[[i]]
      ##Works on maxLik, lm, glm, game, and surv
      if("maxLik" %in% class(mod.i)){ ##Check if the class is maxLik
        out.L<-str_c(out.L,
                     paste(" & ",
                           "\\multicolumn{1}{c}{",
                           num2str(logLik(mod.i)),
                           "}"))
        out.A<-str_c(out.A, " & ")
        out.Lam<-str_c(out.Lam, " & ")
        out.Rho<-str_c(out.Rho, " & ")
        out.R<-str_c(out.R, " & ")
        out.N<- str_c(out.N,
                      paste(" & ",
                            "\\multicolumn{1}{c}{", 
                            nrow(mod.i$gradientObs),
                            "}"))
      }
      if("lm" == class(mod.i)[1]){
        out.L<-str_c(out.L, " & ")
        out.A<-str_c(out.A, " & ")
        out.Lam<-str_c(out.Lam, " & ")
        out.Rho<-str_c(out.Rho, " & ")
        out.R<-str_c(out.R,
                     paste(" & ",
                           "\\multicolumn{1}{c}{", 
                           num2str(summary(mod.i)$adj.r.squared)),
                     "}")
        out.N<-str_c(out.N,
                     paste("& ", 
                           "\\multicolumn{1}{c}{",
                           nrow(mod.i$model),
                           "}"))
      }
      if("ivreg" == class(mod.i)[1] ){
        out.L<-str_c(out.L, " & ")
        out.A<-str_c(out.A, " & ")
        out.Lam<-str_c(out.Lam, " & ")
        out.Rho<-str_c(out.Rho, " & ")
        out.R<-str_c(out.R, " & ")
        out.N<- str_c(out.N,
                      paste(" & ", 
                            "\\multicolumn{1}{c}{",
                            nrow(mod.i$model),
                            "}"))
      }
      
      if("glm" == class(mod.i)[1]  | "game" %in% class(mod.i) 
         |"polr" %in% class(mod.i)){
        out.L<-str_c(out.L,
                     paste(" & ",
                           "\\multicolumn{1}{c}{",
                           num2str(as.numeric(logLik(mod.i))),
                           "}"))
        out.A<-str_c(out.A, " & ")
        out.Lam<-str_c(out.Lam, " & ")
        out.Rho<-str_c(out.Rho, " & ")
        out.R<-str_c(out.R, " & ")
        out.N<- str_c(out.N,
                      paste(" & ",   
                            "\\multicolumn{1}{c}{",
                            nrow(mod.i$model),
                            "}"))
      }
      if("negbin" == class(mod.i)[1]){
        out.L<-str_c(out.L,
                     paste(" & ",
                           "\\multicolumn{1}{c}{",
                           num2str(mod.i$twologlik/2),
                           "}"))                
        out.A <- str_c(out.A, 
                       paste(" & ",    
                             "\\multicolumn{1}{c}{",
                             num2str(1/mod.i$theta),
                             "}"))
        out.Lam<-str_c(out.Lam, " & ")
        out.Rho<-str_c(out.Rho, " & ")
        out.R<-str_c(out.R, " & ")
        out.N<- str_c(out.N,
                      paste(" & ",
                            "\\multicolumn{1}{c}{", 
                            nrow(mod.i$model),
                            "}"))
      }
      if("polywog" == class(mod.i)[1]){
        out.L<-str_c(out.L, " & ")
        out.A<-str_c(out.A, " & ")
        out.Lam<-str_c(out.Lam, 
                       paste(" & ",                                       
                             "\\multicolumn{1}{c}{", 
                             num2str(mod.i$lambda),
                             "}"))
        out.Rho<-str_c(out.Rho, " & ")
        out.R<-str_c(out.R, " & ")
        out.N<- str_c(out.N,
                      paste(" & ",
                            "\\multicolumn{1}{c}{", 
                            mod.i$nobs,
                            "}"))
      }
      if("spautolm" == class(mod.i)[1]){
        out.L<-str_c(out.L,
                     paste(" & ",
                           "\\multicolumn{1}{c}{",
                           num2str(mod.i$LL),
                           "}"))
        out.A<-str_c(out.A, " & ")
        out.Lam<-str_c(out.Lam, " & ")
        out.Rho<-str_c(out.Rho, 
                       paste(" & ",                                       
                             "\\multicolumn{1}{c}{", 
                             num2str(mod.i$lambda),
                             "}"))
        out.R<-str_c(out.R, " & ")
        out.N<- str_c(out.N,
                      paste(" & ",                                       
                            "\\multicolumn{1}{c}{", 
                            nrow(mod.i$X),
                            "}"))
      }
      if("list" == class(mod.i)[1]){
        out.L<-str_c(out.L,
                     paste(" & ",                                       
                           "\\multicolumn{1}{c}{",
                           num2str(mod.i$logLik),
                           "}"))
        out.A<-str_c(out.A, " & ")
        out.Lam<-str_c(out.Lam, " & ")
        out.Rho<-str_c(out.Rho, " & ")
        out.R<-str_c(out.R, " & ")
        out.N<- str_c(out.N,
                      paste(" & ",                                       
                            "\\multicolumn{1}{c}{", 
                            mod.i$Obs,
                            "}"))
      }
      modNames <- str_c(modNames, "& ", model.names[i])
      modNums <- str_c(modNums, "& ", model.numbers[i])
      depNames <- str_c(depNames, " & ", dep.varnames[i])
    }   
    
    if(all("lm" %in% lapply(mod.list, class))){
      out<-str_c(out.R, "\\\\\n", 
                 out.N,   "\\\\ \\hline  \n")
    }else{
      ifelse(any("negbin" %in% sapply(mod.list, class)),
             ifelse(any("spautolm" %in% lapply(mod.list, class)),
                    ifelse(any("polywog" %in% lapply(mod.list, class)),
                           ifelse(any("lm" %in% lapply(mod.list, class)),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.A, "\\\\\n",
                                             out.R, "\\\\\n", 
                                             out.Lam, "\\\\\n", 
                                             out.Rho, "\\\\\n",
                                             out.N,   "\\\\ \\hline  \n"),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.A, "\\\\\n",
                                             out.Lam, "\\\\\n", 
                                             out.Rho, "\\\\\n",
                                             out.N,   "\\\\ \\hline  \n")),
                           ifelse(any("lm" %in% lapply(mod.list, class)),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.A, "\\\\\n",
                                             out.Rho, "\\\\\n",
                                             out.R, "\\\\\n", 
                                             out.N,   "\\\\ \\hline  \n"),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.A, "\\\\\n",
                                             out.Rho, "\\\\\n",
                                             out.N,   "\\\\ \\hline  \n"))),
                    ifelse(any("polywog" %in% lapply(mod.list, class)),
                           ifelse(any("lm" %in% lapply(mod.list, class)),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.A, "\\\\\n",
                                             out.R, "\\\\\n", 
                                             out.Lam, "\\\\\n", 
                                             out.N,   "\\\\ \\hline  \n"),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.A, "\\\\\n",
                                             out.Lam, "\\\\\n",
                                             out.N,   "\\\\ \\hline  \n")),
                           ifelse(any("lm" %in% lapply(mod.list, class)),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.A, "\\\\\n",
                                             out.R, "\\\\\n", 
                                             out.N,   "\\\\ \\hline  \n"),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.A, "\\\\\n",
                                             out.N,   "\\\\ \\hline  \n")))),
             ifelse(any("spautolm" %in% lapply(mod.list, class)),
                    ifelse(any("polywog" %in% lapply(mod.list, class)),
                           ifelse(any("lm" %in% lapply(mod.list, class)),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.R, "\\\\\n", 
                                             out.Lam, "\\\\\n", 
                                             out.Rho, "\\\\\n",
                                             out.N,   "\\\\ \\hline  \n"),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.Lam, "\\\\\n", 
                                             out.Rho, "\\\\\n",
                                             out.N,   "\\\\ \\hline  \n")),
                           ifelse(any("lm" %in% lapply(mod.list, class)),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.Rho, "\\\\\n",
                                             out.R, "\\\\\n", 
                                             out.N,   "\\\\ \\hline  \n"),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.Rho, "\\\\\n",
                                             out.N,   "\\\\ \\hline  \n"))),
                    ifelse(any("polywog" %in% lapply(mod.list, class)),
                           ifelse(any("lm" %in% lapply(mod.list, class)),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.R, "\\\\\n", 
                                             out.Lam, "\\\\\n", 
                                             out.N,   "\\\\ \\hline  \n"),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.Lam, "\\\\\n",
                                             out.N,   "\\\\ \\hline  \n")),
                           ifelse(any("lm" %in% lapply(mod.list, class)),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.R, "\\\\\n", 
                                             out.N,   "\\\\ \\hline  \n"),
                                  out<-str_c(out.L, "\\\\\n", 
                                             out.N,   "\\\\ \\hline  \n")))))
    } 
    
    stars <- ifelse(stars=="all",
                    paste("\\multicolumn{", 
                          1+length(mod.list),
                          "}{l}{\\footnotesize{\\emph{Notes:} $^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^\\dagger p<0.1$ }} \\\\ \\multicolumn{",
                          1+length(mod.list),
                          "}{l}{\\footnotesize{", 
                          seNote,
                          "}}%",
                          sep=""),
                    ifelse(stars=="default",
                           paste("\\multicolumn{", 
                                 1+length(mod.list),
                                 "}{l}{\\footnotesize{\\emph{Notes:} $^{*}p<0.05$}}\\\\ \\multicolumn{",
                                 1+length(mod.list),
                                 "}{l}{\\footnotesize{", 
                                 seNote,
                                 "}}%",
                                 sep=""),
                           ifelse(stars=="none",
                                  paste("\\multicolumn{", 
                                        1+length(mod.list),
                                        "}{l}{\\footnotesize{", 
                                        seNote,
                                        "}}%",
                                        sep=""),
                                  paste("\\multicolumn{", 
                                        1+length(mod.list),
                                        "}{l}{\\footnotesize{\\emph{Notes:}  $^{*}$",
                                        bootCI.level*100,
                                        "\\% CI does not contain 0}}\\\\ \\multicolumn{",
                                        1+length(mod.list),
                                        "}{l}{\\footnotesize{", 
                                        seNote,
                                        "}}%",
                                        sep=""))))


header <- str_c(depNames, " \\\\ \n ", 
                modNames, " \\\\ \n",
                modNums, " \\\\ \n")



out<-list(header=header,
          out=str_c(out, stars), 
          length=length)
if(!is.null(length)){
  out$length <- c(0, length)
}
  }
  if(is.null(out)){stop("Object class not supported")}
  return(out)
}

##handrolled Coeftest
coeftest2<-function(beta, sd){
  z<-beta/sd
  p.value<-pnorm(abs(z), lower.tail=FALSE)*2
  results<-cbind(beta, sd,z, p.value)
  colnames(results)<-c("Estimate",
                       "Std. Error",
                       "Z Value",
                       "Pr(|z|)")
  return(results)
}

##A wrapper for xtable to allow stackign multiple models and inputing 
##Bootstrapped standard errors
xtable2<-function(mod.list, 
                  se=NULL, 
                  info.list=NULL, #Only if coef=TRUE
                  coef=FALSE,
                  apsr=FALSE, 
                  bootCI=FALSE,
                  bootCI.level = .95,
                  stars=c("default", "all", "CI", "none"),
                  caption=NULL,
                  label=NULL, 
                  align=NULL, 
                  digits=2,
                  order=NULL,
                  seNote="Standard Errors in Parenthesis",
                  covariate.labels=NULL,
                  model.numbers=NULL,
                  model.names=NULL,
                  dep.varnames=NULL,         
                  k=1,
                  length=NULL,
                  print.xtable.options=list()
){
  
  stars <- match.arg(stars)
  
  
  if(coef){
    if(length(dep.varnames)==length(info.list)){
      dep.varnames <- paste("\\multicolumn{1}{c}{", dep.varnames, "}")
    }
    if(length(model.names)==length(info.list)){
      model.names <- paste("\\multicolumn{1}{c}{", model.names, "}")
    }
    if(length(model.numbers)==length(info.list)){
      model.numbers <- paste("\\multicolumn{1}{c}{", model.numbers, "}")
    }
    
    info <- mod.sum(info.list, apsr, stars, seNote, bootCI.level, model.numbers, model.names, dep.varnames,k, length)
  }else{
    if(length(dep.varnames)==length(mod.list)){
      dep.varnames <- paste("\\multicolumn{1}{c}{", dep.varnames, "}")
    }
    if(length(model.names)==length(mod.list)){
      model.names <- paste("\\multicolumn{1}{c}{", model.names, "}")
    }
    if(length(model.numbers)==length(mod.list)){
      model.numbers <- paste("\\multicolumn{1}{c}{", model.numbers, "}")
    }
    info <- mod.sum(mod.list, apsr, stars, seNote,bootCI.level, model.numbers, model.names, dep.varnames,k, length)
  }
  
  if(apsr){
    print.xtable.default <- list(booktabs=TRUE,
                                 sanitize.text.function=function(x){x},
                                 include.rownames=FALSE,
                                 caption.placement="top",
                                 table.placement="h!",
                                 include.colnames =FALSE,
                                 add.to.row=list(pos=list(info[[3]][1],info[[3]][2]),
                                                 command=c(info[[1]], info[[2]]))
    )
  }else{
    print.xtable.default <-list(booktabs=TRUE,
                                sanitize.text.function=function(x){x},
                                include.rownames=FALSE,
                                caption.placement="top",
                                table.placement="h!",
                                include.colnames =FALSE,
                                add.to.row=list(pos=list(info[[2]]),
                                                command=info[[1]])
    )
    
  }
  
  print.xtable.options <- modifyList(print.xtable.default, print.xtable.options)
  
  
  
  if(!coef){
    coef.list<-lapply(mod.list, coef)
  }else{
    if(bootCI){
      coef.list<-lapply(mod.list, apply, 2, mean, na.rm=TRUE)
    }else{
      coef.list<-mod.list
    }
  }
  if(is.null(order)){
    order <- unique(unlist(lapply(coef.list, names)))
    
  }else{
    order <- unique(c(order, unlist(lapply(coef.list, names))))
  }
  order <- rep(order, each=2)
  order[seq(2, length(order), by=2)] <- paste(order[seq(2, length(order), by=2)], 
                                              "_SE", 
                                              sep="")
  
  output<-data.frame()
  bootCI.tail <- (1-bootCI.level)/2
  
  
  for(i in 1:nrow(summary(mod.list))){
    
    if(!is.null(se)){
      table<-coeftest2(coef.list[[i]], se[[i]])
    }else{
      if(bootCI){
        table <- cbind(coef.list[[i]], 
                       t(apply(mod.list[[i]],
                               2, 
                               quantile,
                               c(bootCI.tail, 1-bootCI.tail),
                               na.rm=TRUE)))
      }else{
        table<-coeftest2(coef.list[[i]], 
                         sqrt(diag(vcov(mod.list[[i]]))))
      }
    }
    
    if(apsr){
      a.table<-data.frame()
      for(j in 1:nrow(table)){
        if(!bootCI){
          temp<-table[j, 1:2]
          temp<-t(t(temp))
          temp<-num2str(temp, digits=digits)
          temp<-matrix(str_trim(temp), nrow=nrow(temp))
          temp[2,]<-paste("(",
                          temp[2,],
                          ")", sep="")
        }else{
          temp <- table[j, c(1:3)]
          temp <- num2str(temp, digits=digits)
          temp[2] <- paste("\\myalign{,}{(",
                           str_c(temp[2:3], 
                                 collapse=", "),
                           ")}",
                           sep="")
          temp <- temp[-3]
          temp <- t(t(temp))
        }
        if(stars=="default"){
          if(table[j, 4]<0.05){
            temp[1,]<- paste(temp[1,] ,"$^*$", sep="")
          }
        }else{
          if(stars=="all"){
            if(table[j, 4]<0.1 & table[j, 4]>0.05){
              temp[1,]<- paste(temp[1,] 
                               ,"$^\\dagger$", sep="")
            }
            if(table[j, 4]<0.05 & table[j, 4]>0.01){
              temp[1,]<- paste(temp[1,],
                               "$^*$", sep="")
            }
            if(table[j, 4]<0.01& table[j, 4]>0.001){
              temp[1,]<- paste(temp[1,],
                               "$^{**}$", sep="")
            }
            if(table[j, 4]<0.001){
              temp[1,]<- paste(temp[1,],
                               "$^{***}$", sep="")
            }
          }else{
            if(stars=="CI"){
              
              if(sign(table[j,2])==sign(table[j,3])){
                temp[1,] <- paste(temp[1,],
                                  "^{*}",
                                  sep="")
              }
              # if(sign(table[j,2])==sign(table[j,3])){
              #   temp[1,] <- paste(temp[1,],
              #                     "^{***}",
              #                     sep="")
              # }else{
              #   if(sign(table[j,4])==sign(table[j,5])){
              #     temp[1,] <- paste(temp[1,],
              #                       "^{**}", 
              #                       sep="")
              #   }else{
              #     if(sign(table[j,6])==sign(table[j,7])){
              #       temp[1,] <- paste(temp[1,],
              #                         "^{*}", 
              #                         sep="")
              #     }
              #   }
              # }
            }
          }
        }
        row.names(temp)<-NULL
        var.name<-row.names(table)[j]
        colnames(temp)<-paste("Model.", i, sep="")
        row.names(temp)<-c(var.name, paste(var.name, "_SE", sep=""))
        a.table<-rbind(a.table, temp)  
      }
      if(i==1){
        a.output<-a.table
        ` `<-row.names(a.output)
        row.names(a.output)<-NULL
        a.output<-cbind(` `, a.output)
        if(is.null(order)){
          a.output$` ` <- factor(a.output$` `, 
                                 levels=as.character(a.output$` `))
        }else{
          a.output$` ` <- factor(a.output$` `, 
                                 levels=order)
        }
      }else{
        ` `<-row.names(a.table)
        row.names(a.table)<-NULL
        a.table<-cbind(` `, a.table) 
        if(is.null(order)){
          a.table$` ` <- factor(a.table$` `, 
                                levels=as.character(a.table$` `))
        }else{
          a.table$` ` <- factor(a.table$` `, 
                                levels=order)
        }
        
        a.output<-merge(a.output, 
                        a.table,
                        by= " ",
                        all=TRUE,
                        sort=TRUE)        
        a.output <- a.output[order(a.output[,1]),]
      }
    }else{
      ` `<-row.names(table)
      row.names(table)<-NULL
      table<-cbind.data.frame(` `, table)
      output<-rbind(output, table)
      int <- output[1,]
      output <- output[2:nrow(output),]
      output <- rbind(output, int)
    }
  }
  if(apsr){
    a.output$` ` <- as.character(a.output$` `)
    a.output$` `[str_detect(a.output$` `, "_SE")]<-""
    output<-a.output
    const <- str_detect(output[,1], fixed("intercept", ignore_case = TRUE)) | 
      str_detect(output[,1], fixed("constant", ignore_case = TRUE))
    const[is.na(const)] <- FALSE
    idx <- which(const)
    if(length(idx)>0){
      int <- output[c(idx, idx+1),]
      output <- output[-c(idx, idx+1),]
      output <- rbind(output, int)
      colnames(output) <- rep("", ncol(output)) 
    }
    
  }
  if(class(covariate.labels)=="character"){
    output[output[,1]!="",1] <- covariate.labels
  }
  align<-ifelse(is.null(align),
                ifelse(apsr==FALSE, 
                       str_c(rep("r", ncol(output)+1), collapse=""),
                       ifelse(!bootCI,
                              str_c(c("rr", rep("S", ncol(output)-1)), collapse=""),
                              str_c(c("rr", rep("d", ncol(output)-1)), collapse=""))),
                align)
  #      return(output)
  suppressWarnings(
    x <- xtable(output,
                caption=caption,
                label=label, 
                align=align, 
                digits=digits))
  print.xtable.options <- modifyList(print.xtable.options, list(x=x))
  do.call(print, print.xtable.options)
  
  if(bootCI){
    cat("Reminder: Add the following to your TeX preamble:\n")
    cat("\\usepackage{dcolumn}\n
\\newcolumntype{d}{D{.}{.}{-1}}\n
\\newcolumntype{,}{D{,}{,\\,}{-1}}\n
\\newcommand*{\\myalign}[2]{\\multicolumn{1}{#1}{#2}}\n")
  }
  if(str_detect(align, "S")){
    cat("Reminder: Add the following to your TeX preamble:\n")
    cat("\\usepackage{siunitx}
\\sisetup{\n
     input-symbols=(),\n
	table-align-text-post = false,\n
	group-digits=false,\n
} ")
  }
}
###Note this function requires 


#\usepackage{siunitx}
#Be in the preamble 

