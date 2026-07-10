## Saves two figures for the closing slides, sourced from the
## multiple-data-streams and forecasting-nowcasting sessions (simulation only,
## no Stan fitting). Run from the repo root.
library("nfidd.nowcasting")
library("dplyr")
library("tidyr")
library("ggplot2")
set.seed(123)
theme_set(theme_minimal(base_size = 16))

figdir <- "sessions/slides/figures"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

## ---- three streams from one trajectory (multiple-data-streams slides) ----
infections <- make_daily_infections(infection_times)$infections
n <- length(infections)
case_delay  <- censored_delay_pmf(rgamma, max = 10, shape = 3,  rate = 1)
death_delay <- censored_delay_pmf(rgamma, max = 21, shape = 16, rate = 1)
ww_delay    <- censored_delay_pmf(rgamma, max = 10, shape = 2,  rate = 1)
obs <- tibble(
  day = seq_len(n),
  infections = infections,
  cases  = rpois(n, 0.4  * convolve_with_delay(infections, case_delay)),
  deaths = rpois(n, 0.05 * convolve_with_delay(infections, death_delay)),
  ww     = rnorm(n, log(2 * convolve_with_delay(infections, ww_delay)), 0.2)
)
p_streams <- obs |>
  pivot_longer(c(infections, cases, deaths, ww),
               names_to = "stream", values_to = "value") |>
  mutate(stream = factor(stream,
    levels = c("infections", "cases", "deaths", "ww"),
    labels = c("Infections (latent)", "Cases", "Deaths", "Wastewater (log)"))) |>
  ggplot(aes(day, value)) +
  geom_line(linewidth = 0.8, colour = "#447099") +
  facet_wrap(~stream, scales = "free_y", nrow = 1) +
  labs(x = "Day", y = NULL)
ggsave(file.path(figdir, "multiple-data-streams.png"), p_streams,
       width = 9, height = 4.2, dpi = 150, bg = "white")

## ---- different views of the data (forecasting-nowcasting session) ----
gen_time_pmf <- make_gen_time_pmf()
ip_pmf <- make_ip_pmf()
onset_df <- simulate_onsets(make_daily_infections(infection_times), ip_pmf)
cutoff <- 61
forecast_horizon <- 14
reporting_delay_pmf <- censored_delay_pmf(rlnorm, max = 15, meanlog = 1, sdlog = 0.5)
reporting_triangle <- onset_df |>
  filter(day <= cutoff + forecast_horizon) |>
  mutate(reporting_delay = list(tibble(d = 0:15, reporting_delay = reporting_delay_pmf))) |>
  unnest(reporting_delay) |>
  mutate(reported_onsets = rpois(n(), onsets * reporting_delay)) |>
  mutate(reported_day = day + d)
filtered_reporting_triangle <- reporting_triangle |> filter(reported_day <= cutoff)
available_onsets <- filtered_reporting_triangle |>
  summarise(available_onsets = sum(reported_onsets), .by = day)
complete_at_horizon <- reporting_triangle |>
  filter(day <= cutoff + forecast_horizon) |>
  summarise(complete_onsets = sum(reported_onsets), .by = day)
complete_threshold <- 14
complete_data <- available_onsets |> filter(day <= cutoff - complete_threshold)
p_views <- ggplot() +
  geom_line(data = complete_at_horizon,
            aes(x = day, y = complete_onsets, colour = "Complete data"), linewidth = 1) +
  geom_line(data = available_onsets,
            aes(x = day, y = available_onsets, colour = "Available now"), linewidth = 1) +
  geom_line(data = complete_data,
            aes(x = day, y = available_onsets, colour = "Complete as of now"), linewidth = 1.2) +
  geom_vline(xintercept = cutoff, linetype = "dotted") +
  geom_vline(xintercept = cutoff - complete_threshold, linetype = "dotted", alpha = 0.5) +
  scale_colour_manual(values = c(
    "Complete data" = "black", "Available now" = "red", "Complete as of now" = "green")) +
  labs(x = "Day", y = "Onsets", title = "Different views of the data",
       colour = "Data type") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(figdir, "nowcasting-forecasting-data.png"), p_views,
       width = 8, height = 5, dpi = 150, bg = "white")

cat("Saved figures to", figdir, "\n")
