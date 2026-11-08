## 02_make_synthetic_data.R

## Generating a synthetic dataset with the same columns as data/processed/daily,csv, to hide my own personal data.

## Output: data/synthetic/daily_synthetic.csv  (+ a short README in that folder)

## How it works
##  1. Predictors (sleep, heart rate, steps, workouts) are drawn from a
##     multivariate normal using only SUMMARY statistics of my real data
##     (means, SDs, correlations). No real daily values are stored in this script.
##  2. Biocharge is then built from a KNOWN formula (below), so any model can be
##     tested on whether it recovers a formula that is known to be true.
##  3. A gap of missing data is added (1 to 14 May) like the real data.
##
## IMPORTANT: the Biocharge formula below is INVENTED for testing. It is NOT
## Zepp's formula. Its coefficients were set to give similar behaviour to my
## real data (about 63 waking, about 19 at end of day, floor at 5).
##
## Needs: tidyverse, lubridate, here. MASS ships with R.

library(tidyverse)
library(lubridate)
library(here)

set.seed(2026)

clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

make_psd <- function(R) {          # rounding can make a correlation matrix very slightly non-PSD
  e  <- eigen(R, symmetric = TRUE)
  R2 <- e$vectors %*% diag(pmax(e$values, 1e-3)) %*% t(e$vectors)
  D  <- diag(1 / sqrt(diag(R2)))
  D %*% R2 %*% D
}

## ---- 1. Summary statistics of the predictors -------------------------------
vars <- c("sleep_total_min", "deep_min", "rem_min", "awake_min", "sleep_onset_hr",
          "sleep_hr_p5", "hr_day_p5", "hr_day_mean", "steps", "workout_latent")
mu   <- c(407.4, 100.7, 87.8, 17.1, 24.77, 56.2, 58.8, 82.2, 9734, 0)
sdv  <- c( 77.2,  23.2, 26.1, 17.2,  1.56,  6.5,  6.7,  5.6, 3509, 1)   # workout_latent is standard normal

R <- matrix(c(
  1.00,  0.56,  0.72,  0.25, -0.72, -0.24, -0.59, -0.62, -0.01,  0.25,
  0.56,  1.00,  0.22,  0.19, -0.42, -0.10, -0.28, -0.25, -0.02,  0.14,
  0.72,  0.22,  1.00,  0.09, -0.47, -0.29, -0.49, -0.52, -0.02,  0.21,
  0.25,  0.19,  0.09,  1.00, -0.33,  0.05, -0.05, -0.06, -0.13, -0.01,
  -0.72, -0.42, -0.47, -0.33,  1.00,  0.15,  0.52,  0.51,  0.00, -0.17,
  -0.24, -0.10, -0.29,  0.05,  0.15,  1.00,  0.29,  0.51,  0.18,  0.12,
  -0.59, -0.28, -0.49, -0.05,  0.52,  0.29,  1.00,  0.78, -0.02, -0.32,
  -0.62, -0.25, -0.52, -0.06,  0.51,  0.51,  0.78,  1.00,  0.21, -0.06,
  -0.01, -0.02, -0.02, -0.13,  0.00,  0.18, -0.02,  0.21,  1.00,  0.49,
  0.25,  0.14,  0.21, -0.01, -0.17,  0.12, -0.32, -0.06,  0.49,  1.00),
  nrow = 10, byrow = TRUE, dimnames = list(vars, vars))

## ---- 2. Simulate the predictors --------------------------------------------
N     <- 270
dates <- seq(as.Date("2026-02-01"), by = "day", length.out = N)

Z <- MASS::mvrnorm(N, mu = rep(0, length(vars)), Sigma = make_psd(R))
colnames(Z) <- vars
X <- as_tibble(sweep(sweep(Z, 2, sdv, "*"), 2, mu, "+"))

daily <- tibble(date = dates) |>
  bind_cols(X) |>
  mutate(
    ## sleep stages (total is rebuilt so that total = deep + light + rem always holds)
    deep_min        = round(clip(deep_min, 30, 185)),
    rem_min         = round(clip(rem_min, 10, 175)),
    awake_min       = round(clip(awake_min, 0, 120)),
    sleep_total_min = round(clip(sleep_total_min, 150, 660)),
    light_min       = pmax(40, sleep_total_min - deep_min - rem_min),
    sleep_total_min = deep_min + light_min + rem_min,
    time_in_bed_min = sleep_total_min + awake_min,
    sleep_onset_hr  = round(clip(sleep_onset_hr, 19, 31), 2),
    sleep_offset_hr = round(sleep_onset_hr - 24 + time_in_bed_min / 60, 2),
    sleep_start     = round_date(as.POSIXct(date, tz = "UTC") + (sleep_onset_hr - 24) * 3600, "minute"),
    sleep_stop      = sleep_start + time_in_bed_min * 60,
    n_naps          = rbinom(N, 1, 0.06),
    
    ## heart rate
    sleep_hr_p5     = round(clip(sleep_hr_p5, 43, 80)),
    sleep_hr_mean   = round(sleep_hr_p5 + pmax(2, rnorm(N, 7.8, 2.9)), 1),
    sleep_hr_sd     = round(clip(rnorm(N, 5.8, 2), 2, 15), 1),
    sleep_hr_n      = round(sleep_total_min * 0.9 + rnorm(N, 0, 15)),
    hr_day_p5       = round(clip(hr_day_p5, 45, 85)),
    hr_day_mean     = round(pmax(clip(hr_day_mean, 66, 105), hr_day_p5 + 8), 1),
    hr_day_max      = round(hr_day_mean + clip(rnorm(N, 60.7, 16.5), 25, 110)),
    hr_day_sd       = round(clip(rnorm(N, 17.5, 3.6), 8, 30), 1),
    hr_awake_mean   = round(hr_day_mean + rnorm(N, 8.1, 3.2), 1),
    hr_n            = round(clip(rnorm(N, 1346, 60), 1200, 1440)),
    
    ## activity
    steps           = round(clip(steps, 1500, 24000)),
    distance_m      = round(steps * rnorm(N, 0.735, 0.029)),
    run_distance_m  = round(distance_m * rnorm(N, 0.784, 0.023)),
    active_records  = round(steps * 0.0234 + rnorm(N, 0, 40)),
    brisk_records   = pmax(0, round(rnorm(N, 0.128 * active_records, 8))),
    peak_steps_record = round(rnorm(N, 119, 6)),
    
    ## workouts (zero on about 40% of days)
    workout_day     = workout_latent > qnorm(0.40),
    workout_min     = if_else(workout_day, round(clip(55 + 50 * workout_latent, 15, 236)), 0),
    n_workouts      = if_else(workout_day, 1L + rpois(N, 0.9), 0L),
    workout_kcal    = round(workout_min * rnorm(N, 8.6, 1.2)),
    workout_dist_m  = if_else(workout_day & runif(N) < 0.26, round(workout_min * runif(N, 60, 150)), 0),
    calories_kcal   = round(pmax(250, 27 + 0.035 * steps + 7.1 * workout_min + rnorm(N, 0, 100)))
  )

## ---- 3. Weight: slow drift, first weigh-in on 12 Feb ------------------------
w <- numeric(N); w[1] <- 87.4
for (t in 2:N) w[t] <- 87.4 + 0.98 * (w[t - 1] - 87.4) + rnorm(1, 0, 0.25)
daily$weight_kg <- round(w, 1)
daily$weight_kg[1:11] <- NA

## ---- 4. Biocharge from a KNOWN formula (invented, see header) --------------
##   waking_t = 7.6 + 0.109*sleep_total + 0.063*deep - 0.203*sleep_hr_p5
##              + 0.845*end_of_day_{t-1} + noise_t
##              (noise is AR(1): phi = 0.37, overall SD 11.5)
##   end_t    = max(5, -37.7 + 0.905*waking_t - 0.465*steps/1000 - 0.017*workout_min
##              + noise(SD 10.1))        # the 5 is the app's floor (censoring)
phi <- 0.37
waking <- end <- numeric(N)
end_prev <- 20
eps <- 0
for (t in seq_len(N)) {
  eps <- phi * eps + rnorm(1, 0, 11.5 * sqrt(1 - phi^2))
  w_raw <- 7.6 + 0.109 * daily$sleep_total_min[t] + 0.063 * daily$deep_min[t] -
    0.203 * daily$sleep_hr_p5[t] + 0.845 * end_prev + eps
  waking[t] <- clip(round(w_raw), 5, 96)
  e_raw <- -37.7 + 0.905 * waking[t] - 0.465 * daily$steps[t] / 1000 -
    0.017 * daily$workout_min[t] + rnorm(1, 0, 10.1)
  end[t]   <- clip(round(e_raw), 5, 55)
  end_prev <- end[t]
}
daily$biocharge_waking <- waking
daily$biocharge_end    <- end

## ---- 5. Missing data like the real file ------------------------------------
gap_idx <- which(dates >= as.Date("2026-05-01") & dates <= as.Date("2026-05-14"))
sleep_cols <- c("deep_min", "light_min", "rem_min", "awake_min", "sleep_total_min", "sleep_onset_hr",
                "sleep_offset_hr", "time_in_bed_min", "sleep_hr_mean", "sleep_hr_p5", "sleep_hr_sd", "sleep_hr_n")
hr_cols    <- c("active_records", "brisk_records", "peak_steps_record", "hr_day_mean", "hr_day_p5",
                "hr_day_max", "hr_day_sd", "hr_awake_mean", "hr_n")

for (cc in c(sleep_cols, hr_cols, "biocharge_waking", "biocharge_end")) daily[[cc]][gap_idx] <- NA
for (cc in c("steps", "distance_m", "run_distance_m", "calories_kcal")) daily[[cc]][gap_idx] <- 0
daily$sleep_start[gap_idx] <- as.POSIXct(daily$date[gap_idx], tz = "UTC") + 23 * 3600   # placeholder, as in the export
daily$sleep_stop[gap_idx]  <- daily$sleep_start[gap_idx]

rest <- setdiff(seq_len(N), gap_idx)
daily$biocharge_waking[sample(rest, 2)] <- NA          # a couple of extra missing scores
miss_end <- sample(rest, 2)
daily$biocharge_end[miss_end] <- NA
for (cc in c("sleep_hr_mean", "sleep_hr_p5", "sleep_hr_sd", "sleep_hr_n")) daily[[cc]][sample(rest, 5)] <- NA

daily$biocharge_end_at_floor <- daily$biocharge_end == 5

## ---- 6. Same columns, same order as daily.csv, and save --------------------
daily <- daily |>
  select(date, biocharge_waking, biocharge_end, biocharge_end_at_floor, sleep_start, sleep_stop,
         deep_min, light_min, rem_min, awake_min, n_naps, sleep_total_min, sleep_onset_hr,
         sleep_offset_hr, time_in_bed_min, sleep_hr_mean, sleep_hr_p5, sleep_hr_sd, sleep_hr_n,
         steps, distance_m, run_distance_m, calories_kcal, active_records, brisk_records,
         peak_steps_record, hr_day_mean, hr_day_p5, hr_day_max, hr_day_sd, hr_awake_mean, hr_n,
         n_workouts, workout_min, workout_kcal, workout_dist_m, weight_kg)

dir.create(here("data", "synthetic"), recursive = TRUE, showWarnings = FALSE)
write_csv(daily, here("data", "synthetic", "daily_synthetic.csv"))

writeLines(c(
  "# Synthetic data",
  "",
  "`daily_synthetic.csv` is simulated by `R/02_make_synthetic_data.R` (seed 2026).",
  "It has the same columns as the real processed file, but contains no real observations.",
  "",
  "- Predictors are drawn from a multivariate normal using summary statistics of the real data.",
  "- Biocharge is built from a KNOWN, invented formula, so models can be checked on whether they recover it.",
  "- The formula is NOT Zepp's. It exists only to test the analysis pipeline.",
  "- 1 to 14 May are missing, as in the real data."
), here("data", "synthetic", "README.md"))

message("Synthetic data: ", nrow(daily), " rows, ", ncol(daily), " columns")
message("Mean waking Biocharge: ", round(mean(daily$biocharge_waking, na.rm = TRUE), 1),
        " | mean end-of-day: ", round(mean(daily$biocharge_end, na.rm = TRUE), 1),
        " | share at floor of 5: ", round(mean(daily$biocharge_end_at_floor, na.rm = TRUE), 2))