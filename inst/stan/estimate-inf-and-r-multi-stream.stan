functions {
  #include "functions/convolve_with_delay.stan"
  #include "functions/renewal.stan"
  #include "functions/geometric_random_walk.stan"
}

data {
  int n;                // number of days
  int I0;               // number initially infected

  // shared latent process
  int gen_time_max;     // maximum generation time
  array[gen_time_max] real gen_time_pmf; // pmf of generation time distribution

  // stream 1: cases (infection-to-report delay, ascertainment)
  array[n] int cases;
  int<lower = 1> case_delay_max;
  array[case_delay_max + 1] real case_delay_pmf;

  // stream 2: deaths (infection-to-death delay, IFR)
  array[n] int deaths;
  int<lower = 1> death_delay_max;
  array[death_delay_max + 1] real death_delay_pmf;

  // stream 3: wastewater (infection-to-shedding delay, scaling)
  array[n] real ww;     // log-scale wastewater concentration
  int<lower = 1> ww_delay_max;
  array[ww_delay_max + 1] real ww_delay_pmf;
}

parameters {
  real<lower = 0> init_R;          // initial reproduction number
  array[n-1] real rw_noise;        // random walk noise
  real<lower = 0> rw_sd;           // random walk standard deviation
  real<lower = 0, upper = 1> ascertainment; // proportion of infections reported as cases
  real<lower = 0, upper = 1> ifr;  // infection fatality ratio
  real<lower = 0> ww_scale;        // wastewater scaling (signal per infection)
  real<lower = 0> ww_sigma;        // wastewater observation standard deviation (log scale)
}

transformed parameters {
  // one shared infection / Rt process feeds every stream
  array[n] real R = geometric_random_walk(init_R, rw_noise, rw_sd);
  array[n] real infections = renewal(I0, R, gen_time_pmf);

  // each stream is a convolution of the SAME infections with its own delay
  array[n] real exp_cases;
  array[n] real exp_deaths;
  array[n] real exp_ww;
  {
    array[n] real case_conv = convolve_with_delay(infections, case_delay_pmf);
    array[n] real death_conv = convolve_with_delay(infections, death_delay_pmf);
    array[n] real ww_conv = convolve_with_delay(infections, ww_delay_pmf);
    for (i in 1:n) {
      exp_cases[i] = ascertainment * case_conv[i];
      exp_deaths[i] = ifr * death_conv[i];
      exp_ww[i] = ww_scale * ww_conv[i];
    }
  }
}

model {
  // priors
  init_R ~ normal(1, 0.5) T[0, ];
  rw_noise ~ std_normal();
  rw_sd ~ normal(0, 0.05) T[0, ];
  ascertainment ~ beta(2, 2);
  ifr ~ beta(1, 50);
  ww_scale ~ normal(1, 1) T[0, ];
  ww_sigma ~ normal(0, 0.5) T[0, ];

  // joint likelihood: each stream contributes its own term off the shared process
  cases ~ poisson(exp_cases);
  deaths ~ poisson(exp_deaths);
  ww ~ normal(log(exp_ww), ww_sigma);
}
