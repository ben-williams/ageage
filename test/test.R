library(Rcpp)
library(ageage)
ageage <- Rcpp::Module("ageage", PACKAGE = "ageage")

interface<-new(ageage$AgeAgeInterface)
interface$sigma1<-10
interface$sigma2<-2