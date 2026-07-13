// renewal equation started from a full initial history rather than a single
// I0, so the generation time window is complete on the first observed day
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
