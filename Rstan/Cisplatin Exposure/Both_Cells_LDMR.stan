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
  real km_rep_CPT = theta[7];
  real kd_rep_CPT = theta[8];
  
  /** Below you can fill in the differental equations of the model*/

  real dDose= y1[2] * 0.022052503115 - y1[1] * (1.13252628565 + 0.1813857857015);
  real dKidneyPt = y1[3] * + 0.0003284364143375 + y1[1] * 1.13252628565  - y1[2] * (0.022052503115 + 0.00181720824448 + 0.001673191758);
  real dAccuPt = y1[2] * 0.001673191758 - y1[3] * 0.0003284364143375;
  real dDD = theta[1] + theta[2] * (y1[2] + y1[3]) - y1[5] * y1[4] * theta[3] - y1[4]*theta[4];
  real dRep = theta[5] + (theta[6]* y1[4])/(theta[7] + y1[4]) - y1[5] * theta[8];

  
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
  real km_rep_CPT = theta2[7];
  real kd_rep_CPT = theta2[8];
  
  /** Below you can fill in the differental equations of the model*/

  real dDose= y2[2] * 0.022052503115 - y2[1] * (1.13252628565 + 0.1813857857015);
  real dKidneyPt = y2[3] * + 0.0003284364143375 + y2[1] * 1.13252628565  - y2[2] * (0.022052503115 + 0.00181720824448 + 0.001673191758);
  real dAccuPt = y2[2] * 0.001673191758 - y2[3] * 0.0003284364143375;
  real dDD = theta2[1] + theta2[2] * (y2[2] + y2[3]) - y2[5] * y2[4] * theta2[3] - y2[4]*theta2[4];
  real dRep = theta2[5] + (theta2[6]* y2[4])/(theta2[7] + y2[4]) - y2[5] * theta2[8];

  
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
      real <lower = ks_DD_CPT>ks_DDkid_CPT;
      real <lower =0>kd_DDrep_CPT;
      real <lower =0>kd_DD_CPT;
      real <lower =0>ks_rep_CPT;
      real <lower =0>ks_repDD_CPT;
      real <lower =0>km_rep_CPT;
      real <lower =0>kd_rep_CPT;
      
      real <lower =0>ks_DD_PPT;
      real <lower = ks_DD_PPT>ks_DDkid_PPT;
      real <lower =0>kd_DDrep_PPT;
      real <lower =0>kd_DD_PPT;
      real <lower =0>ks_rep_PPT;
      real <lower =0>ks_repDD_PPT;
      real <lower =0>km_rep_PPT;
      real <lower =0>kd_rep_PPT;
  }

transformed parameters {
  
  
  real y1_hat[N_CPT,5];
  real y2_hat[N_PPT,5];
  
  {   real theta[8] = {ks_DD_CPT,ks_DDkid_CPT,kd_DDrep_CPT,kd_DD_CPT,ks_rep_CPT,ks_repDD_CPT,km_rep_CPT,kd_rep_CPT};
    y1_hat = integrate_ode_rk45(cpt, y0_CPT, t0_CPT, ts_CPT, theta, x_r, x_i);
    }
    
  {   real theta2[8] = {ks_DD_PPT,ks_DDkid_PPT,kd_DDrep_PPT,kd_DD_PPT,ks_rep_PPT,ks_repDD_PPT,km_rep_PPT,kd_rep_PPT};
    y2_hat = integrate_ode_rk45(ppt, y0_PPT, t0_PPT, ts_PPT, theta2, x_r, x_i);
    }
  
}


model {
  // Priors
  
      //ks_DD_CPT ~ exponential(1);
      //ks_DDkid_CPT ~ exponential(1);
      //kd_DDrep_CPT ~ exponential(1);
      //kd_DD_CPT ~ exponential(1);
      //ks_rep_CPT ~ exponential(1);
      //ks_repDD_CPT ~ exponential(1);
      //km_rep_CPT ~ exponential(1);
      //kd_rep_CPT ~ exponential(1);
      
      //ks_DD_PPT ~ exponential(1);
      //ks_DDkid_PPT ~ exponential(1);
      //kd_DDrep_PPT ~ exponential(1);
      //kd_DD_PPT ~ exponential(1);
      //ks_rep_PPT ~ exponential(1);
      //ks_repDD_PPT ~ exponential(1);
      //km_rep_PPT ~ exponential(1);
      //kd_rep_PPT ~ exponential(1);
      
      
      // Estimate After 500 iterations of inference 
      
      //ks_DD_CPT ~ normal(0.005215722449,0.005121025);
      //ks_DDkid_CPT ~ normal(0.0003156098954,0.0003118735);
      //kd_DDrep_CPT ~ normal(0.044981422949855,0.04307045);
      //kd_DD_CPT ~ normal(0.089454983305,0.08805995);
      //ks_rep_CPT ~ normal(0.003033013669323,0.003091405);
      //ks_repDD_CPT ~ normal(0.01066388463075,0.01032065);
      //km_rep_CPT ~ normal(2.555986846015,2.463875);
      //kd_rep_CPT ~ normal(0.007662119844,0.00761897);
      //ks_DD_PPT ~ normal(0.00058358581,0.0005835195);
      //ks_DDkid_PPT ~ normal(0.00014845004685,0.0001479905);
      //kd_DDrep_PPT ~ normal(0.0120857502825,0.01209365);
      //kd_DD_PPT ~ normal(0.013408715016,0.0133661);
      //ks_rep_PPT ~ normal(0.00102103997648575,0.0007401595);
      //ks_repDD_PPT ~ normal(1.3884391493,0.966107);
      //km_rep_PPT ~ normal(6.52645063975,6.398135);
      //kd_rep_PPT ~ normal(0.085593543605,0.0493742);
      
      ks_DD_CPT ~ normal(0.000704274126,0.0007050185);
      ks_DDkid_CPT ~ normal(0.000759480635,0.0007538585);
      kd_DDrep_CPT ~ normal(0.00098652523882,0.000709777);
      kd_DD_CPT ~ normal(0.026796613,0.0265124);
      ks_rep_CPT ~ normal(0.00034141419094,0.00036242);
      ks_repDD_CPT ~ normal(0.0002054594352826,0.000161667);
      km_rep_CPT ~ normal(1.656964177578,1.619985);
      kd_rep_CPT ~ normal(0.002014444994,0.001987295);
      ks_DD_PPT ~ normal(0.00015526913077,0.000137243);
      ks_DDkid_PPT ~ normal(0.000724523571,0.0007198015);
      kd_DDrep_PPT ~ normal(0.0022777745016,0.00204878);
      kd_DD_PPT ~ normal(0.0146287286,0.01453385);
      ks_rep_PPT ~ normal(0.00153598861868,0.00135716);
      ks_repDD_PPT ~ normal(3.360715368,3.30083);
      km_rep_PPT ~ normal(19.94343563,18.77315);
      kd_rep_PPT ~ normal(1.5137038249,1.430065);
  

  // Likelihood
  for (n in 1:N_CPT) {
    // real y_hat_sum;
    // y_hat_sum = (y_hat[n,2]+y_hat[n,3])*scale;
    
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
  
     
      
   //target += -0.1*abs(ks_DD_CPT - ks_DD_PPT)^2;
   //target += -0.1*abs(ks_DDkid_CPT - ks_DDkid_PPT)^2;
   //target += -0.1*abs(kd_DDrep_CPT - kd_DDrep_PPT)^2;
   //target += -0.1*abs(kd_DD_CPT - kd_DD_PPT)^2;
   //target += -0.1*abs(ks_rep_CPT - ks_rep_PPT)^2;
   //target += -0.1*abs(ks_repDD_CPT - ks_repDD_PPT)^2;
   //target += -0.1*abs(km_rep_CPT - km_rep_PPT)^2;
   //target += -0.1*abs(kd_rep_CPT - kd_rep_PPT)^2;

}

