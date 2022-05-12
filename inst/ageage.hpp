#if !defined(_ageage_)
#  define _ageage_

class model_data : public ad_comm{
  data_int nobs;
  data_vector age;
  data_vector ape;
  data_vector n;
  ~model_data();
  model_data(int argc,char * argv[]);
  friend class model_parameters;
};

class model_parameters : public model_data ,
  public function_minimizer
{
public:
  ~model_parameters();
  void preliminary_calculations(void);
  void set_runtime(void);
  virtual void * mycast(void) {return (void*)this;}
  static int mc_phase(void)
  {
    return initial_params::mc_phase;
  }
  static int mceval_phase(void)
  {
    return initial_params::mceval_phase;
  }
  static int sd_phase(void)
  {
    return initial_params::sd_phase;
  }
  static int current_phase(void)
  {
    return initial_params::current_phase;
  }
  static int last_phase(void)
  {
    return (initial_params::current_phase
      >=initial_params::max_number_phases);
  }
  static prevariable current_feval(void)
  {
    return *objective_function_value::pobjfun;
  }
private:
  ivector integer_control_flags;
  dvector double_control_flags;
  param_init_number sigma1;
  param_init_number sigma2;
  param_vector Perc_Corr;
  param_vector Perc_Corr1;
  param_vector Perc_Corr2;
  param_vector Phat;
  param_stddev_vector sigma_a;
  param_number sigma_inc;
  param_vector RSS;
  param_number prior_function_value;
  param_number likelihood_function_value;
  objective_function_value f;
public:
  virtual void userfunction(void);
  virtual void report(const dvector& gradients);
  virtual void final_calcs(void);
  model_parameters(int sz,int argc, char * argv[]);
  virtual void initializationfunction(void);
  void get_A_SD_est(void);
  void evaluate_the_objective_function(void);

};
#endif
