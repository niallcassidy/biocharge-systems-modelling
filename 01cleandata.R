
  
  
  
##  01 clean data · R
## 01_clean_data.R
## Reads the raw Zepp exports in data/raw/ plus the hand-entered Biocharge
## scores, and builds ONE ROW PER DAY.
##
## Output: data/processed/daily.csv
##         figures/biocharge_over_time.png
##
## Needs: tidyverse, lubridate, here   (dplyr >= 1.1.0 for join_by)
## Run from the project root (open the .Rproj in RStudio).

library(tidyverse)
library(lubridate)
library(here)

## ---- Settings --------------------------------------------------------------
raw_dir     <- here("data", "raw")
BRISK_STEPS <- 100   # steps in one activity record to count as "brisk" (arbitrary; change freely)

## ---- Helpers ---------------------------------------------------------------
# Find exactly one file in data/raw matching a pattern (export names carry a timestamp)
find_file <- function(pattern) {
  f <- list.files(raw_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(f) != 1) {
    stop("Expected exactly one file matching '", pattern, "' in ", raw_dir,
         " but found ", length(f), ": ", paste(basename(f), collapse = ", "))
  }
  f
}

# Build a timestamp from a date and an "HH:MM" string.
# Zepp writes midnight at the end of a day as "24:00", so work in minutes.
make_ts <- function(date, time) {
  mins <- as.integer(substr(time, 1, 2)) * 60L + as.integer(substr(time, 4, 5))
  as.POSIXct(as.character(date), tz = "UTC") + mins * 60
}

## ---- 1. Daily activity (steps, distance, calories) -------------------------
activity <- read_csv(find_file("^ACTIVITY_[0-9]+\\.csv$"), show_col_types = FALSE) |>
  transmute(date = as_date(date),
            steps,
            distance_m     = distance,      # assumed metres
            run_distance_m = runDistance,
            calories_kcal  = calories)

## ---- 2. Sleep (one row per night, dated by the wake-up day) ----------------
# The 'naps' column holds JSON with backslash-escaped quotes, so tell readr.
sleep <- read_delim(find_file("^SLEEP_[0-9]+\\.csv$"), delim = ",",
                    escape_backslash = TRUE, escape_double = FALSE,
                    col_types = cols(start = col_character(),   # keep as text: readr would otherwise turn
                                     stop  = col_character(),   # these into datetimes and drop exact-midnight times
                                     naps  = col_character()),
                    show_col_types = FALSE) |>
  transmute(date        = as_date(date),
            sleep_start = ymd_hms(substr(start, 1, 19), tz = "UTC"),
            sleep_stop  = ymd_hms(substr(stop,  1, 19), tz = "UTC"),
            deep_min    = deepSleepTime,
            light_min   = shallowSleepTime,
            rem_min     = REMTime,
            awake_min   = wakeTime,
            n_naps      = str_count(coalesce(naps, ""), "start")) |>
  mutate(sleep_total_min = deep_min + light_min + rem_min,
         no_sleep_data   = sleep_total_min == 0,   # a night of all zeros is missing data, not zero sleep
         # clock times as hours after midnight; evening starts (before 12:00 = after midnight) get +24
         sleep_onset_hr  = hour(sleep_start) + minute(sleep_start) / 60,
         sleep_onset_hr  = if_else(sleep_onset_hr < 12, sleep_onset_hr + 24, sleep_onset_hr),
         sleep_offset_hr = hour(sleep_stop) + minute(sleep_stop) / 60,
         time_in_bed_min = as.numeric(difftime(sleep_stop, sleep_start, units = "mins"))) |>
  mutate(across(c(deep_min, light_min, rem_min, awake_min, sleep_total_min,
                  sleep_onset_hr, sleep_offset_hr, time_in_bed_min),
                ~ if_else(no_sleep_data, NA_real_, as.numeric(.x)))) |>
  select(-no_sleep_data)

## ---- 3. Heart rate during sleep (from the minute-level sleep file) ---------
# In SLEEP_MINUTE.csv the 'date' label is one day AHEAD of the calendar date of
# the reading (checked against SLEEP.csv), so the real timestamp is (date - 1) + time.
# Each reading is then assigned to the night whose start/stop window contains it.
sleep_hr <- read_csv(find_file("^SLEEP_MINUTE_[0-9]+\\.csv$"),
                     col_types = cols(time = col_character()), show_col_types = FALSE) |>
  transmute(ts = make_ts(as_date(date) - 1, time), hr) |>
  filter(!is.na(hr)) |>
  inner_join(select(sleep, date, sleep_start, sleep_stop),
             join_by(between(ts, sleep_start, sleep_stop))) |>
  group_by(date) |>
  summarise(sleep_hr_mean = mean(hr),
            sleep_hr_p5   = quantile(hr, 0.05, names = FALSE),  # robust "lowest" overnight HR
            sleep_hr_sd   = sd(hr),
            sleep_hr_n    = n(),
            .groups = "drop")

## ---- 4. Heart rate across the day ------------------------------------------
sleep_windows <- select(sleep, sleep_date = date, sleep_start, sleep_stop)

hr_daily <- read_csv(find_file("^HEARTRATE_AUTO_[0-9]+\\.csv$"),
                     col_types = cols(time = col_character()), show_col_types = FALSE) |>
  transmute(ts = make_ts(date, time), heartRate) |>
  distinct(ts, .keep_all = TRUE) |>               # "24:00" duplicates the next day's 00:00
  filter(heartRate > 0) |>
  left_join(sleep_windows, join_by(between(ts, sleep_start, sleep_stop))) |>
  mutate(asleep = !is.na(sleep_date),
         day    = as_date(ts)) |>
  group_by(date = day) |>
  summarise(hr_day_mean   = mean(heartRate),
            hr_day_p5     = quantile(heartRate, 0.05, names = FALSE),  # resting-HR proxy
            hr_day_max    = max(heartRate),
            hr_day_sd     = sd(heartRate),
            hr_awake_mean = mean(heartRate[!asleep]),
            hr_n          = n(),                                       # readings that day (max about 1440)
            .groups = "drop") |>
  mutate(across(where(is.numeric), ~ if_else(is.nan(.x), NA_real_, .x)))

## ---- 5. Minute-level activity ----------------------------------------------
activity_detail <- read_csv(find_file("^ACTIVITY_MINUTE_[0-9]+\\.csv$"),
                            col_types = cols(time = col_character()), show_col_types = FALSE) |>
  mutate(day = as_date(make_ts(date, time))) |>
  group_by(date = day) |>
  summarise(active_records   = sum(steps > 0),
            brisk_records    = sum(steps >= BRISK_STEPS),
            peak_steps_record = max(steps),
            .groups = "drop")

## ---- 6. Workouts -----------------------------------------------------------
# Workout 'type' codes are not decoded here; all workouts are pooled per day.
# (startTime is stored as +0000, so a late-evening workout could land on the wrong day.)
workouts <- read_csv(find_file("^SPORT_[0-9]+\\.csv$"),
                     col_types = cols(startTime = col_character()), show_col_types = FALSE) |>
  mutate(date = as_date(substr(startTime, 1, 10))) |>
  group_by(date) |>
  summarise(n_workouts     = n(),
            workout_min    = sum(`sportTime(s)`, na.rm = TRUE) / 60,
            workout_kcal   = sum(`calories(kcal)`, na.rm = TRUE),
            workout_dist_m = sum(`distance(m)`, na.rm = TRUE),
            .groups = "drop")

## ---- 7. Body weight (sparse: only a few weigh-ins) -------------------------
body <- read_csv(find_file("^BODY_[0-9]+\\.csv$"),
                 col_types = cols(time = col_character()), show_col_types = FALSE) |>
  transmute(date = as_date(substr(time, 1, 10)), weight_kg = weight) |>
  group_by(date) |>
  summarise(weight_kg = mean(weight_kg), .groups = "drop")

## ---- 8. Biocharge scores (the targets) -------------------------------------
# Numbers adds a title line ("biocharge") above the headers, so skip it if present.
bio_file <- find_file("^biocharge.*\\.csv$")
skip_n   <- if (str_detect(tolower(read_lines(bio_file, n_max = 1)), "date")) 0 else 1

biocharge <- read_csv(bio_file, skip = skip_n, show_col_types = FALSE)
names(biocharge) <- c("date", "biocharge_waking", "biocharge_end")
biocharge <- biocharge |>
  mutate(date = as_date(date),
         biocharge_end_at_floor = biocharge_end == 5)   # 5 looks like the app's minimum (censored)

## ---- 9. Join everything: one row per day -----------------------------------
daily <- tibble(date = seq(min(activity$date), max(activity$date), by = "day")) |>
  left_join(biocharge,       by = "date") |>
  left_join(sleep,           by = "date") |>
  left_join(sleep_hr,        by = "date") |>
  left_join(activity,        by = "date") |>
  left_join(activity_detail, by = "date") |>
  left_join(hr_daily,        by = "date") |>
  left_join(workouts,        by = "date") |>
  left_join(body,            by = "date") |>
  arrange(date) |>
  fill(weight_kg, .direction = "down") |>                        # carry the last weigh-in forward
  mutate(across(c(n_workouts, workout_min, workout_kcal, workout_dist_m),
                ~ replace_na(.x, 0)))                            # no workout recorded = 0

## ---- 10. Quality checks and save -------------------------------------------
stopifnot(!anyDuplicated(daily$date))
message("Rows (days): ", nrow(daily), "  from ", min(daily$date), " to ", max(daily$date))
message("Missing values per column:")
print(colSums(is.na(daily))[colSums(is.na(daily)) > 0])

dir.create(here("data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("figures"), showWarnings = FALSE)
write_csv(daily, here("data", "processed", "daily.csv"))

daily |>
  select(date, biocharge_waking, biocharge_end) |>
  pivot_longer(-date, names_to = "score", values_to = "value") |>
  ggplot(aes(date, value, colour = score)) +
  geom_line(na.rm = TRUE) +
  geom_point(size = 0.8, na.rm = TRUE) +
  labs(title = "Biocharge over time", x = NULL, y = "Biocharge", colour = NULL) +
  theme_minimal()
ggsave(here("figures", "biocharge_over_time.png"), width = 9, height = 4.5, dpi = 150)


