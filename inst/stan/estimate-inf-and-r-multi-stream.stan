functions {
  #include "functions/convolve_with_delay.stan"
  #include "functions/renewal.stan"
  #include "functions/geometric_random_walk.stan"

  // Renewal equation with a seeded initial history. The shared renewal()
  // starts from a single scalar I0, which is fine at the start of an epidemic
  // but collapses to near zero when (as here) we start mid-outbreak: the first
  // few days only "see" that single seed through a truncated generation-time
  // window. This variant instead takes a full initial history (`seed`, one
  // value per day) so the convolution window is complete from the first
  // observed day and the estimated infections start at the right level. The
  // recursion is otherwise identical to renewal().
  array[] real renewal_seeded(array[] real seed, array[] real R,
                              array[] real gen_time) {
    int seed_n = num_elements(seed);
    int n = num_elements(R);
    int max_gen_time = num_elements(gen_time); // gen_time starts at day 1
    array[seed_n + n] real I;
    I[1:seed_n] = seed;
    for (i in 1:n) {
      int t = seed_n + i;                       // current day in full series
      int first = max(1, t - max_gen_time);
      int len = t - first;                      // past days contributing
      array[len] real segment = I[first:(t - 1)];
      array[len] real gen_pmf = reverse(gen_time[1:len]);
      I[t] = dot_product(segment, gen_pmf) * R[i];
    }
    return I[(seed_n + 1):(seed_n + n)];        // observed window only
  }
}

data {
  int n;                // number of days
  int I0;               // number initially infected

  // shared latent process
  int gen_time_max;     // maximum generation time
  array[gen_time_max] real gen_time_pmf; // pmf of generation time distribution

  // per-stream switches: set to 1 to fit a stream, 0 to leave it out.
  // This lets the SAME model fit one stream on its own, two streams linked
  // together, or all three jointly, so we can build the model up in parts.
  int<lower = 0, upper = 1> use_cases;
  int<lower = 0, upper = 1> use_deaths;
  int<lower = 0, upper = 1> use_ww;

  // switch on a time-varying death scaling (a random walk on the IFR) to let
  // the model explain a drifting severity rather than forcing a compromise.
  int<lower = 0, upper = 1> tv_death_scale;

  // use an overdispersed (negative binomial) observation model for cases;
  // deaths stay Poisson and wastewater log-normal. Cases are the timely but
  // noisiest stream, so overdispersion matters most there.
  int<lower = 0, upper = 1> use_nb_cases;

  // priors for each stream's scaling, as (mean, sd) pairs. A caller can pass a
  // TIGHT prior near the truth to a single-stream diagnostic fit (which only
  // identifies infections x scaling, so needs the scale informed to recover the
  // level) or a RELAXED prior to the joint fit (where several streams anchor
  // the level between them).
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
  // seed for the initial infection history: a level near the start of the
  // window and an exponential growth rate, used to build a full history over
  // gen_time_max days so the renewal starts at the right level, not near zero.
  real<lower = 0> seed_base;       // infections at the start of the window
  real initial_growth;             // daily (log) growth rate over the seed
  real<lower = 0> init_R;          // initial reproduction number
  array[n-1] real rw_noise;        // random walk noise
  real<lower = 0> rw_sd;           // random walk standard deviation
  real<lower = 0, upper = 1> ascertainment; // proportion of infections reported
  real<lower = 0, upper = 1> ifr;  // infection fatality ratio (initial level)
  real<lower = 0> ww_scale;        // wastewater scaling (signal per infection)
  real<lower = 0> ww_sigma;        // wastewater obs sd (log scale)
  // reciprocal overdispersion for the negative-binomial cases, 1 / sqrt(phi):
  // near 0 this recovers the Poisson, larger values allow more overdispersion.
  real<lower = 0> cases_overdispersion;
  // optional random walk on the (log) death scaling, only used when
  // tv_death_scale = 1. Sized to zero otherwise so it costs nothing.
  array[tv_death_scale ? n - 1 : 0] real ifr_rw_noise;
  array[tv_death_scale ? 1 : 0] real<lower = 0> ifr_rw_sd;
}

transformed parameters {
  // one shared infection / Rt process feeds every stream
  array[n] real R = geometric_random_walk(init_R, rw_noise, rw_sd);
  // seed a full gen_time_max-day initial history by exponential growth, so the
  // renewal starts at the right level instead of collapsing from a single I0.
  array[gen_time_max] real seed_infections;
  for (t in 1:gen_time_max) {
    seed_infections[t] =
      seed_base * exp(initial_growth * (t - gen_time_max));
  }
  array[n] real infections =
    renewal_seeded(seed_infections, R, gen_time_pmf);

  // death scaling: either constant (ifr) or a geometric random walk starting
  // from ifr, which lets a drifting severity be absorbed rather than forced
  // into the shared infections.
  array[n] real death_scale;
  if (tv_death_scale) {
    death_scale = geometric_random_walk(ifr, ifr_rw_noise, ifr_rw_sd[1]);
  } else {
    for (i in 1:n) death_scale[i] = ifr;
  }

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
      exp_deaths[i] = death_scale[i] * death_conv[i];
      exp_ww[i] = ww_scale * ww_conv[i];
    }
  }
}

model {
  // priors
  // seed the initial history near the observed starting level I0, with a
  // weak prior on the initial growth rate.
  seed_base ~ normal(I0, I0) T[0, ];
  initial_growth ~ normal(0, 0.2);
  init_R ~ normal(1, 0.5) T[0, ];
  rw_noise ~ std_normal();
  rw_sd ~ normal(0, 0.05) T[0, ];
  // configurable scaling priors: a tight (mean, sd) near the truth anchors the
  // level in a single-stream diagnostic fit; a relaxed one lets the streams
  // pin the level between them in the joint fit.
  ascertainment ~ normal(ascertainment_p[1], ascertainment_p[2]) T[0, 1];
  ifr ~ normal(ifr_p[1], ifr_p[2]) T[0, 1];
  ww_scale ~ normal(ww_scale_p[1], ww_scale_p[2]) T[0, ];
  ww_sigma ~ normal(0, 0.5) T[0, ];
  cases_overdispersion ~ normal(0, 1) T[0, ];
  if (tv_death_scale) {
    ifr_rw_noise ~ std_normal();
    ifr_rw_sd ~ normal(0, 0.1);   // half-normal via <lower = 0> on parameter
  }

  // joint likelihood: each stream contributes its own term off infections,
  // and only the streams we switch on are added. Because the streams are
  // conditionally independent given the infections, the joint log-likelihood
  // is simply the sum of the per-stream terms.
  if (use_cases) {
    if (use_nb_cases) {
      // negative binomial: variance = mu + mu^2 * cases_overdispersion^2
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
  // posterior predictive draws WITH observation error, from each stream's own
  // likelihood. These are the full predictive distribution (not the expected
  // signal exp_*), so posterior predictive checks compare like with like:
  // cases as negative-binomial (or Poisson) counts, deaths as Poisson counts,
  // wastewater on the log scale as normal draws.
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
