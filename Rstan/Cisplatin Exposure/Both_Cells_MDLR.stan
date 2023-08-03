functions {
real[] cpt(real t, real[] y1, real[] theta, real[] x_r, int[] x_i) {
  
  // Define the state variables
  
  real Dose = y1[1];
  real KidneyPt = y1[2];
  real AccuPt = y1[3];
  real DD = y1[4];
  real Rep = y1[5];
  
  real ks_DD_CPT = theta[1];
  real ks_DDkid_CPT = theta[2];
  real kd_DDrep_CPT = theta[3];
  real kd_DD_CPT = theta[4];
  real ks_rep_CPT = theta[5];
  real ks_repDD_CPT = theta[6];
  real kd_rep_CPT = theta[7];
  real km_dd_CPT = theta[8];
  
  /** Below you can fill in the differental equations of the model*/

  real dDose= y1[2] * 0.022052503115 - y1[1] * (1.13252628565 + 0.1813857857015);
  real dKidneyPt = y1[3] * + 0.0003284364143375 + y1[1] * 1.13252628565  - y1[2] * (0.022052503115 + 0.00181720824448 + 0.001673191758);
  real dAccuPt = y1[2] * 0.001673191758 - y1[3] * 0.0003284364143375;
  real dDD = theta[1] + (theta[2] * (y1[2] + y1[3]))/( theta[8] + (y1[2] + y1[3])) - y1[5] * y1[4] * theta[3] - y1[4]*theta[4];
  real dRep = theta[5] + theta[6]* y1[4] - y1[5] * theta[7];

  
  return {dDose, dKidneyPt,dAccuPt,dDD,dRep};
}

real[] ppt(real t, real[] y2, real[] theta2, real[] x_r, int[] x_i) {
  
  // Define the state variables
  real Dose = y2[1];
  real KidneyPt = y2[2];
  real AccuPt = y2[3];
  real DD = y2[4];
  real Rep = y2[5];
  
  real ks_DD_CPT = theta2[1];
  real ks_DDkid_CPT = theta2[2];
  real kd_DDrep_CPT = theta2[3];
  real kd_DD_CPT = theta2[4];
  real ks_rep_CPT = theta2[5];
  real ks_repDD_CPT = theta2[6];
  real kd_rep_CPT = theta2[7];
  real km_dd_PPT = theta2[8];
  
  /** Below you can fill in the differental equations of the model*/

  real dDose= y2[2] * 0.022052503115 - y2[1] * (1.13252628565 + 0.1813857857015);
  real dKidneyPt = y2[3] * + 0.0003284364143375 + y2[1] * 1.13252628565  - y2[2] * (0.022052503115 + 0.00181720824448 + 0.001673191758);
  real dAccuPt = y2[2] * 0.001673191758 - y2[3] * 0.0003284364143375;
  real dDD = theta2[1] + (theta2[2] * (y2[2] + y2[3]))/( theta2[8] + (y2[2] + y2[3])) - y2[5] * y2[4] * theta2[3] - y2[4]*theta2[4];
  real dRep = theta2[5] + theta2[6]* y2[4] - y2[5] * theta2[7];

  
  return {dDose, dKidneyPt,dAccuPt,dDD,dRep};
}
}

data {
  
  //CPT
  
  int N_CPT;           // Number of time points
  real t0_CPT;        
  real ts_CPT[N_CPT];                // Time points
  real y1[N_CPT,2];              // Observed state
  real sigma1[N_CPT,2]; // Known standard deviations of the observations
  real y0_CPT[5];
  
  //PPT
  
  int N_PPT;           // Number of time points
  real t0_PPT;        
  real ts_PPT[N_PPT];                // Time points
  real y2[N_PPT,2];              // Observed state
  real sigma2[N_PPT,2]; // Known standard deviations of the observations
  real y0_PPT[5];
}


transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {

      real <lower =0>ks_DD_CPT;
      real <lower =ks_DD_CPT>ks_DDkid_CPT;
      real <lower =0>kd_DDrep_CPT;
      real <lower =0>kd_DD_CPT;
      real <lower =0>ks_rep_CPT;
      real <lower =0>ks_repDD_CPT;
      real <lower =0>kd_rep_CPT;
      real <lower =0>km_dd_CPT;
      
      real <lower =0>ks_DD_PPT;
      real <lower =ks_DD_PPT>ks_DDkid_PPT;
      real <lower =0>kd_DDrep_PPT;
      real <lower =0>kd_DD_PPT;
      real <lower =0>ks_rep_PPT;
      real <lower =0>ks_repDD_PPT;
      real <lower =0>kd_rep_PPT;
      real <lower =0>km_dd_PPT;
  }

transformed parameters {
  
  
  real y1_hat[N_CPT,5];
  real y2_hat[N_PPT,5];
  
  {   real theta[8] = {ks_DD_CPT,ks_DDkid_CPT,kd_DDrep_CPT,kd_DD_CPT,ks_rep_CPT,ks_repDD_CPT,kd_rep_CPT,km_dd_CPT};
    y1_hat = integrate_ode_rk45(cpt, y0_CPT, t0_CPT, ts_CPT, theta, x_r, x_i);
    }
    
  {   real theta2[8] = {ks_DD_PPT,ks_DDkid_PPT,kd_DDrep_PPT,kd_DD_PPT,ks_rep_PPT,ks_repDD_PPT,kd_rep_PPT,km_dd_PPT};
    y2_hat = integrate_ode_rk45(ppt, y0_PPT, t0_PPT, ts_PPT, theta2, x_r, x_i);
    }
  
}


model {
      // Initial Priors
  
      //ks_DD_CPT ~ exponential(1);
      //ks_DDkid_CPT ~ exponential(1);
      //kd_DDrep_CPT ~ exponential(1);
      //kd_DD_CPT ~ exponential(1);
      //ks_rep_CPT ~ exponential(1);
      //ks_repDD_CPT ~ exponential(1);
      //kd_rep_CPT ~ exponential(1);
      //km_dd_CPT ~ exponential(1);
      
      //ks_DD_PPT ~ exponential(1);
      //ks_DDkid_PPT ~ exponential(1);
      //kd_DDrep_PPT ~ exponential(1);
      //kd_DD_PPT ~ exponential(1);
      //ks_rep_PPT ~ exponential(1);
      //ks_repDD_PPT ~ exponential(1);
      //kd_rep_PPT ~ exponential(1);
      //km_dd_PPT ~ exponential(1);
      
      // Estimate After 500 iterations of inference 
      
      //ks_DD_CPT ~ normal(0.00287497815547,0.002051015);
      //ks_DDkid_CPT ~ normal(0.012580947823,0.0127486);
      //kd_DDrep_CPT ~ normal(0.203068029,0.203323);
      //kd_DD_CPT ~ normal(0.0566460704,0.05506955);
      //ks_rep_CPT ~ normal(0.00191226279153,0.00199301);
      //ks_repDD_CPT ~ normal(0.0081144773591,0.007481125);
      //kd_rep_CPT ~ normal(0.00549823287,0.005421455);
      //km_dd_CPT ~ normal(3.50556669336,3.32082);
      //ks_DD_PPT ~ normal(2.68378947233e-05,1.736e-05);
      //ks_DDkid_PPT ~ normal(0.01111632744,0.01094395);
      //kd_DDrep_PPT ~ normal(0.00036681907214,0.0002559995);
      //kd_DD_PPT ~ normal(0.0209860358,0.02065895);
      //ks_rep_PPT ~ normal(0.00347996437597,0.00138043);
      //ks_repDD_PPT ~ normal(1.3750359777,1.116715);
      //kd_rep_PPT ~ normal(0.471837293,0.3805365);
      //km_dd_PPT ~ normal(37.3319374,36.98925);
      
ks_DD_CPT ~ normal(0.0495799306,0.0505882);
      ks_DDkid_CPT ~ normal(0.0777890694,0.0755196);
      kd_DDrep_CPT ~ normal(1.52498794,1.507955);
      kd_DD_CPT ~ normal(0.0821785272,0.0815538);
      ks_rep_CPT ~ normal(0.00100544757,0.00100731);
      ks_repDD_CPT ~ normal(0.0001926558375,0.0001852185);
      kd_rep_CPT ~ normal(0.00830756177,0.008314);
      km_dd_CPT ~ normal(0.51481286667,0.4362555);
      ks_DD_PPT ~ normal(0.000131400635479,0.0001192335);
      ks_DDkid_PPT ~ normal(0.253051961,0.251054);
      kd_DDrep_PPT ~ normal(0.00154795117116,0.00140219);
      kd_DD_PPT ~ normal(0.0159841378,0.01581185);
      ks_rep_PPT ~ normal(0.00080801904993,0.000687473);
      ks_repDD_PPT ~ normal(0.09109991482,0.08792275);
      kd_rep_PPT ~ normal(0.7617531017,0.741532);
      km_dd_PPT ~ normal(307.502638,302.921);

  // Likelihood
  for (n in 1:N_CPT) {
    
    
    // Codes to investigate values per iteration in case of errors
    
    //print("y_hat1 for time point ", n, ": ", y_hat[n,1]);
    //print("y_hat2 for time point ", n, ": ", y_hat[n,2]);
    //print("y_hat3 for time point ", n, ": ", y_hat[n,3]);
    //print("y_hat_sum for time point ", n, ": ", y_hat_sum);
    //print("theta for time point ", n, ": ", {ks_DD,ks_DDkid,kd_DDrep,kd_DD,ks_rep,ks_repDD,km_rep,kd_rep,n_rep});
    //print("Likelihood for y_hat[", n, "4] = ", normal_lpdf( y[n,1] | y_hat[n,4], sigma[n,1]));
    //print("Likelihood for y_hat[", n, "5] = ", normal_lpdf( y[n,2] | y_hat[n,5], sigma[n,2]));

    
   y1[n,1] ~ normal(y1_hat[n,4],sigma1[n,1]);
   y1[n,2] ~ normal(y1_hat[n,5], sigma1[n,2]);

  }
  
    for (n in 1:N_PPT) {

   y2[n,1] ~ normal(y2_hat[n,4],sigma2[n,1]);
   y2[n,2] ~ normal(y2_hat[n,5], sigma2[n,2]);

  }
  
}

