library(ageage)
library(Rcpp)
ageage<-Rcpp::Module("ageage", PACKAGE="ageage")

model<-new(ageage$AgeAgeModel)
model$nobs<-10
model$age<-c(1,2,3,4,5,6,7,8,9,10)
model$ape<-c(1,2,3,4,5,6,7,8,9,10)
model$n<-c(1,2,3,4,5,6,7,8,9,10)
argv<-c("")
model$Run(argv)
