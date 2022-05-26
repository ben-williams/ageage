#ifndef RINTERFACE
#define RINTERFACE
#include <Rcpp.h>

class AgeAgeInterface{
public:
//data section
 int nobs;
 Rcpp::NumericVector age;
 Rcpp::NumericVector ape;
 Rcpp::NumericVector n;

 //parameter section
 //Rcpp::IntegerVector integer_control_flags;
 // dvector double_control_flags;

  double sigma1; //param_init_number
  double sigma2; //param_init_number
//   param_vector Perc_Corr;
//   param_vector Perc_Corr1;
//   param_vector Perc_Corr2;
//   param_vector Phat;
//   param_stddev_vector sigma_a;
//   double sigma_inc;//param_number
//    Rcpp::NumericVector RSS;
//   //double prior_function_value;//
//  // double likelihood_function_value;//param_number

};

RCPP_EXPOSED_CLASS(AgeAgeInterface)
RCPP_MODULE(ageage) {
    Rcpp::class_<AgeAgeInterface>("AgeAgeInterface")
            .constructor()
            .field("sigma1", &AgeAgeInterface::sigma1)
            .field("sigma2", &AgeAgeInterface::sigma2);

};

#endif