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
  real km_rep_CPT = theta[9];
  
  /** Below you can fill in the differental equations of the model*/

  real dDose= y1[2] * 0.022052503115 - y1[1] * (1.13252628565 + 0.1813857857015);
  real dKidneyPt = y1[3] * + 0.0003284364143375 + y1[1] * 1.13252628565  - y1[2] * (0.022052503115 + 0.00181720824448 + 0.001673191758);
  real dAccuPt = y1[2] * 0.001673191758 - y1[3] * 0.0003284364143375;
  real dDD = theta[1] + (theta[2] * (y1[2] + y1[3]))/( theta[8] + (y1[2] + y1[3])) - y1[5] * y1[4] * theta[3] - y1[4]*theta[4];
  real dRep = theta[5] + (theta[6]* y1[4])/(theta[9] + y1[4]) - y1[5] * theta[7];

  
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
  real km_rep_PPT = theta2[9];
  
  /** Below you can fill in the differental equations of the model*/

  real dDose= y2[2] * 0.022052503115 - y2[1] * (1.13252628565 + 0.1813857857015);
  real dKidneyPt = y2[3] * + 0.0003284364143375 + y2[1] * 1.13252628565  - y2[2] * (0.022052503115 + 0.00181720824448 + 0.001673191758);
  real dAccuPt = y2[2] * 0.001673191758 - y2[3] * 0.0003284364143375;
  real dDD = theta2[1] + (theta2[2] * (y2[2] + y2[3]))/( theta2[8] + (y2[2] + y2[3])) - y2[5] * y2[4] * theta2[3] - y2[4]*theta2[4];
  real dRep = theta2[5] + (theta2[6]* y2[4])/(theta2[9] + y2[4]) - y2[5] * theta2[7];

  
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
      real <lower =0>km_rep_CPT;
      
      real <lower =0>ks_DD_PPT;
      real <lower =ks_DD_PPT>ks_DDkid_PPT;
      real <lower =0>kd_DDrep_PPT;
      real <lower =0>kd_DD_PPT;
      real <lower =0>ks_rep_PPT;
      real <lower =0>ks_repDD_PPT;
      real <lower =0>kd_rep_PPT;
      real <lower =0>km_dd_PPT;
      real <lower =0>km_rep_PPT;
  }

transformed parameters {
  
  
  real y1_hat[N_CPT,5];
  real y2_hat[N_PPT,5];
  
  {   real theta[9] = {ks_DD_CPT,ks_DDkid_CPT,kd_DDrep_CPT,kd_DD_CPT,ks_rep_CPT,ks_repDD_CPT,kd_rep_CPT,km_dd_CPT,km_rep_CPT};
    y1_hat = integrate_ode_rk45(cpt, y0_CPT, t0_CPT, ts_CPT, theta, x_r, x_i);
    }
    
  {   real theta2[9] = {ks_DD_PPT,ks_DDkid_PPT,kd_DDrep_PPT,kd_DD_PPT,ks_rep_PPT,ks_repDD_PPT,kd_rep_PPT,km_dd_PPT,km_rep_PPT};
    y2_hat = integrate_ode_rk45(ppt, y0_PPT, t0_PPT, ts_PPT, theta2, x_r, x_i);
    }
  
}


model {
  // Priors
  
      //ks_DD_CPT ~ exponential(1);
      //ks_DDkid_CPT ~ exponential(1);
    //  kd_DDrep_CPT ~ exponential(1);
      //kd_DD_CPT ~ exponential(1);
      //ks_rep_CPT ~ exponential(1);
      //ks_repDD_CPT ~ exponential(1);
      //kd_rep_CPT ~ exponential(1);
      //km_dd_CPT ~ exponential(1);
      //km_rep_CPT ~ exponential(1);
      
      //ks_DD_PPT ~ exponential(1);
      //ks_DDkid_PPT ~ exponential(1);
      //kd_DDrep_PPT ~ exponential(1);
      //kd_DD_PPT ~ exponential(1);
      //ks_rep_PPT ~ exponential(1);
      //ks_repDD_PPT ~ exponential(1);
      //kd_rep_PPT ~ exponential(1);
      //km_dd_PPT ~ exponential(1);
      //km_rep_PPT ~ exponential(1);
      
      
      // Estimate After 500 iterations of inference 
      
      //ks_DD_CPT ~ normal(0.00291289335,0.002956995);
      //ks_DDkid_CPT ~ normal(0.01493226385,0.01428785);
      //kd_DDrep_CPT ~ normal(0.185274107,0.185052);
      //kd_DD_CPT ~ normal(0.0657925491,0.06193685);
      //ks_rep_CPT ~ normal(0.00198922894,0.001679615);
      //ks_repDD_CPT ~ normal(0.01450072882,0.01489885);
      //kd_rep_CPT ~ normal(0.00623109639,0.005805755);
      //km_dd_CPT ~ normal(6.5572490186,6.822505);
      //km_rep_CPT ~ normal(1.48862689,1.355035);
      //ks_DD_PPT ~ normal(2.114788998e-05,1.919905e-05);
      //ks_DDkid_PPT ~ normal(0.0243772413,0.02446685);
      //kd_DDrep_PPT ~ normal(0.0004006536167,0.0003452485);
      //kd_DD_PPT ~ normal(0.0180947661,0.01784555);
      //ks_rep_PPT ~ normal(0.000788126916,0.000908736);
      //ks_repDD_PPT ~ normal(10.54172152,11.18905);
      //kd_rep_PPT ~ normal(0.0212288078,0.02197885);
      //km_dd_PPT ~ normal(112.193088,114.61);
      //km_rep_PPT ~ normal(1112.034377,1106.425);
      
ks_DD_CPT ~ normal(0.04670520274,0.04827715);
      ks_DDkid_CPT ~ normal(0.0799981729,0.07837555);
      kd_DDrep_CPT ~ normal(1.51886506,1.49916);
      kd_DD_CPT ~ normal(0.0815001023,0.08080075);
      ks_rep_CPT ~ normal(0.000852225971,0.000880784);
      ks_repDD_CPT ~ normal(0.00075950766622,0.000731869);
      kd_rep_CPT ~ normal(0.00825063161,0.00826385);
      km_dd_CPT ~ normal(0.48480164747,0.413111);
      km_rep_CPT ~ normal(1.41520015206,1.33467);
      ks_DD_PPT ~ normal(0.00013239622052,0.000119158);
      ks_DDkid_PPT ~ normal(0.265409031,0.260441);
      kd_DDrep_PPT ~ normal(0.00170609995478,0.001537985);
      kd_DD_PPT ~ normal(0.0164544809,0.016369);
      ks_rep_PPT ~ normal(0.00042872055691,0.0003714145);
      ks_repDD_PPT ~ normal(0.5603652843,0.535639);
      kd_rep_PPT ~ normal(0.2847765068,0.2626275);
      km_dd_PPT ~ normal(311.906617,303.066);
      km_rep_PPT ~ normal(17.74124034,16.61145);
  

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
  
     
      
   //target += 0.1*abs(ks_DD_CPT - ks_DD_PPT)^2;
   //target += 0.1*abs(ks_DDkid_CPT - ks_DDkid_PPT)^2;
   //target += 0.1*abs(kd_DDrep_CPT - kd_DDrep_PPT)^2;
   //target += 0.1*abs(kd_DD_CPT - kd_DD_PPT)^2;
   //target += 0.1*abs(ks_rep_CPT - ks_rep_PPT)^2;
   //target += 0.1*abs(ks_repDD_CPT - ks_repDD_PPT)^2;
   //target += 0.1*abs(km_rep_CPT - km_rep_PPT)^2;
   //target += 0.1*abs(kd_rep_CPT - kd_rep_PPT)^2;

}

