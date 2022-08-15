install.packages("devtools") #any version should be fine

library(devtools)
install_github("ccrismancox/games2")

install_version("doParallel", "1.0.16", upgrade = "always")
install_version("doRNG", "1.8.2", upgrade = "always")
install_version("brglm", "0.7.2", upgrade = "always")
install_version("detectseparation", "0.2", upgrade = "always")
install_version("matrixStats", "0.60.1", upgrade = "always")
install_version("Formula", "1.2-4", upgrade = "always")
install_version("knitr", "1.33", upgrade = "always")
install_version("numDeriv", "2016.8-1.1", upgrade = "always")
install_version("readstata13", "0.10.0", upgrade = "always")
install_version("ggplot2", "3.3.6", upgrade = "always")
install_version("gridExtra", "2.3", upgrade = "always")
install_version("stringr", "1.4.0", upgrade = "always")
install_version("mc2d", "0.1-21", upgrade = "always")


library(games2)
library(doRNG)
library(brglm)
library(detectseparation)
library(matrixStats)
library(Formula)
library(knitr)
library(numDeriv)
library(readstata13)
library(ggplot2)
library(gridExtra)
library(stringr)
library(mc2d)

print(sessionInfo())