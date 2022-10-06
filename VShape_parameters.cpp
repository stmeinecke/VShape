struct parameters{
  //class to collect all parameters of the system.
  
  //define parameters here:
  double alpha_Q;
  double alpha_G;
  double kappa;
  double delta;
  double J_Q;
  double J_G;
  double gamma_Q;
  double gamma_G;
  double r_s;
  double T;
  double T01;
  double T12;
  double TR;
  double omega;
  
  double Q_0;
  
  //noise
  double nStr;
  double SqrtNoiseStr;
  
  //delay
  double tau;
  double k;
  double C;
  
  parameters(){}
}; 
