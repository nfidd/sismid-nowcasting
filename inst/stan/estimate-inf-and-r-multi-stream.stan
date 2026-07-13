functions {
  #include "functions/convolve_with_delay.stan"
  #include "functions/renewal_seeded.stan"
  #include "functions/geometric_random_walk.stan"
}

data {
  int n;                // number of days
  int I0;               // number initially infected

  // shared latent process
  int gen_time_max;     // maximum generation time
  array[gen_time_max] real gen_time_pmf; // pmf of generation time distribution

  // per-stream switches: set to 1 to fit a stream, 0 to leave it out
  int<lower = 0, upper = 1> use_cases;
  int<lower = 0, upper = 1> use_deaths;
  int<lower = 0, upper = 1> use_ww;

  // switch on a random walk on the IFR rather than a constant IFR
  int<lower = 0, upper = 1> tv_death_scale;

  // switch cases to a negative binomial observation model
  int<lower = 0, upper = 1> use_nb_cases;

  // priors for each stream's scaling, as (mean, sd) pairs
  array[2] real ascertainment_p; // ascertainment prior, truncated to [0, 1]
  array[2] real ifr_p;           // IFR prior, truncated to [0, 1]
  array[2] real ww_scale_p;      // wastewater scaling prior, truncated to [0, ]

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
  real<lower = 0> seed_base;       // infections at the start of the window
  real initial_growth;             // daily (log) growth rate over the seed
  real<lower = 0> init_R;          // initial reproduction number
  array[n-1] real rw_noise;        // random walk noise
  real<lower = 0> rw_sd;           // random walk standard deviation
  real<lower = 0, upper = 1> ascertainment; // proportion of infections reported
  real<lower = 0, upper = 1> ifr;  // infection fatality ratio (initial level)
  real<lower = 0> ww_scale;        // wastewater scaling (signal per infection)
  real<lower = 0> ww_sigma;        // wastewater obs sd (log scale)
  // reciprocal overdispersion for the cases, 1 / sqrt(phi)
  real<lower = 0> cases_overdispersion;
  // random walk on the death scaling, sized to zero when not used
  array[tv_death_scale ? n - 1 : 0] real ifr_rw_noise;
  array[tv_death_scale ? 1 : 0] real<lower = 0> ifr_rw_sd;
}

transformed parameters {
  // one shared infection / Rt process feeds every stream
  array[n] real R = geometric_random_walk(init_R, rw_noise, rw_sd);
  // initial history over gen_time_max days, grown from seed_base, so the
  // renewal has a complete generation time window on the first observed day
  array[gen_time_max] real seed_infections;
  for (t in 1:gen_time_max) {
    seed_infections[t] =
      seed_base * exp(initial_growth * (t - gen_time_max));
  }
  array[n] real infections =
    renewal_seeded(seed_infections, R, gen_time_pmf);

  // death scaling: constant ifr, or a random walk starting from it
  array[n] real death_scale;
  if (tv_death_scale) {
    death_scale = geometric_random_walk(ifr, ifr_rw_noise, ifr_rw_sd[1]);
  } else {
    for (i in 1:n) death_scale[i] = ifr;
  }

  // each stream convolves the same infections with its own delay
  array[n] real exp_cases;
  array[n] real exp_deaths;
  array[n] real exp_ww;
  {
    array[n] real case_conv = convolve_with_delay(infections, case_delay_pmf);
    array[n] real death_conv = convolve_with_delay(infections, death_delay_pmf);
    array[n] real ww_conv = convolve_with_delay(infections, ww_delay_pmf);
    for (i in 1:n) {
      exp_cases[i] = ascertainment * case_conv[i];
      exp_deaths[i] = death_scale[i] * death_conv[i];
      exp_ww[i] = ww_scale * ww_conv[i];
    }
  }
}

model {
  // priors
  seed_base ~ normal(I0, I0) T[0, ];
  initial_growth ~ normal(0, 0.2);
  init_R ~ normal(1, 0.5) T[0, ];
  rw_noise ~ std_normal();
  rw_sd ~ normal(0, 0.05) T[0, ];
  ascertainment ~ normal(ascertainment_p[1], ascertainment_p[2]) T[0, 1];
  ifr ~ normal(ifr_p[1], ifr_p[2]) T[0, 1];
  ww_scale ~ normal(ww_scale_p[1], ww_scale_p[2]) T[0, ];
  ww_sigma ~ normal(0, 0.5) T[0, ];
  cases_overdispersion ~ normal(0, 1) T[0, ];
  if (tv_death_scale) {
    ifr_rw_noise ~ std_normal();
    ifr_rw_sd ~ normal(0, 0.1);   // half-normal via <lower = 0> on parameter
  }

  // the streams are conditionally independent given the infections, so the
  // joint log-likelihood is the sum of the terms we switch on
  if (use_cases) {
    if (use_nb_cases) {
      // variance = mu + mu^2 * cases_overdispersion^2
      cases ~ neg_binomial_2(exp_cases, inv_square(cases_overdispersion));
    } else {
      cases ~ poisson(exp_cases);
    }
  }
  if (use_deaths) {
    deaths ~ poisson(exp_deaths);
  }
  if (use_ww) {
    ww ~ normal(log(exp_ww), ww_sigma);
  }
}

generated quantities {
  // posterior predictive draws, with observation error, for each stream
  array[n] int pp_cases;
  array[n] int pp_deaths;
  array[n] real pp_ww;
  for (i in 1:n) {
    if (use_nb_cases) {
      pp_cases[i] = neg_binomial_2_rng(
        exp_cases[i], inv_square(cases_overdispersion)
      );
    } else {
      pp_cases[i] = poisson_rng(exp_cases[i]);
    }
    pp_deaths[i] = poisson_rng(exp_deaths[i]);
    pp_ww[i] = normal_rng(log(exp_ww[i]), ww_sigma);
  }
}
