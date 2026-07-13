# nfidd.nowcasting 1.3.0

- Split the combined SISMID course into a standalone nowcasting course and renamed the package from `nfidd` to `nfidd.nowcasting`. Forecasting sessions, slides, datasets, and Stan-free helpers now live in the companion `nfidd.forecasting` package.
- Corrected the expectation model in the complex reporting session. `enw_expectation()` takes a formula for the growth rate, so `~ rw(day)` put a random walk on the growth rate rather than the geometric random walk on expected counts the session described. The daily models now use a random effect by day and the weekly models a random effect by week.
- Made the simulated cases in the joint session informative. The negative binomial observation model was correct, but the simulated cases were too noisy to constrain the latent infections, so the random walk prior on $R_t$ dominated and the fit drifted above the truth.
- Reordered the sessions to match the timetable, so joint fitting comes before complex reporting, and renumbered the session ordering.
- Extracted `renewal_seeded()` into its own Stan functions file.
- Added learning objectives for the complex reporting session.

# nfidd 1.2.0

- Localised the course material to the SISMID course

# nfidd 1.1.2

- added the R version of `condition_onsets_by_report()`

# nfidd 1.1.1

- adapted `nfidd_cmdstan_model()` to work with an include path option, and a model file name argument

# nfidd 1.1.0

- fixed a bug in convolution function which affected the earliest part of convoluted time series #475.
- renamed `target_day` to `origin_day` for clarity #465
- added `nffid_sample()` function to speed up default inference #457
- replaced `vapply` with for loop in `convolve_with_delay` for clarity #433
- streamlined the use of logged and natural R #424
- added the `summarise_lognormal()` function for mean/sd summarises #406

# nfidd 1.0.0

In development version of the package and teaching material for teaching in Bangkok in November 2024.

This included a complete redevelopment of the package where what previously were snippets are now functions.

# nfidd 0.1.0

Initial release of the `nfidd` package and teaching material for teaching in June 2024 in Stockholm.
