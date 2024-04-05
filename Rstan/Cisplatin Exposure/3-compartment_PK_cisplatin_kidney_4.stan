functions {
real[] ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
  // Define the state variables
  real Dose = y[1];
  real KidneyPt = y[2];
  real AccuPt = y[3];

  /** Below you can fill in the differental equations of the model*/

  real dDose= y[2] * theta[3] - y[1] * (theta[1] + theta[2]);
  real dKidneyPt = y[3] * + theta[6] + y[1] * theta[1]  - y[2] * (theta[3] + theta[4] + theta[5]);
  real dAccuPt = y[2] * theta[5] - y[3] * theta[6];
  
  // theta[1] = k1
  // theta[2] = ke_dose
  // theta[3] = k-1
  // theta[4] = ke_kidney
    
  return {dDose, dKidneyPt,dAccuPt};
}
}

data {
  int<lower=1> N;           // Number of time points
  real t0;        
  real ts[N];                // Time points
  real y[N,2];              // Observed state
  real<lower=0> sigma[N,2]; // Known standard deviations of the observations
  real y0[3];
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real<lower=0> theta[6];   // Parameters to be inferred
  real scale;
  
  }

transformed parameters {
  
    real y_hat[N,3];
    y_hat = integrate_ode_bdf(ode, y0, t0, ts, theta, x_r, x_i);
  
}


model {
  // Priors
  theta[1] ~ normal(1.132526e+00,0.0521909772);
  theta[2] ~ normal(1.813858e-01,0.0492411816	);
  theta[3] ~ normal(2.205250e-02,0.0014384697	);
  theta[4] ~ normal(1.817208e-03,0.0008922234);
  theta[5] ~ normal(1.673192e-03,0.0001968134);
  theta[6] ~ normal(3.284364e-04,0.0001026295);
  scale ~ normal(4.313859e+00,0.2256496508);

  // Likelihood
  for (n in 1:N) {
    real y_hat_sum;
    y_hat_sum = (y_hat[n,2]+y_hat[n,3])*scale;
    
    //print("y_hat1 for time point ", n, ": ", y_hat[n,1]);
    //print("y_hat2 for time point ", n, ": ", y_hat[n,2]);
    //print("y_hat3 for time point ", n, ": ", y_hat[n,3]);
    //print("y_hat_sum for time point ", n, ": ", y_hat_sum);
    //print("theta for time point ", n, ": ", theta);
    //print("Likelihood for y_hat[", n, "1] = ", normal_lpdf( y[n,1] | y_hat[n,1], sigma[n,1]));
    //print("Likelihood for y_hat[", n, "1] = ", normal_lpdf( y[n,2] | y_hat[n,2]+y_hat[n,3], sigma[n,2]));
    
   y[n,1] ~ normal(y_hat[n,1], sigma[n,1]);
   y[n,2] ~ normal(y_hat_sum, sigma[n,2]);
  }
}

