# ── 0. Libraries ──────────────────────────────────────────────────────────────
library(tidyverse)   # data wrangling + ggplot2
library(zoo)         # yearmon class — clean monthly ts handling
library(tseries)     # later: adf.test, jarque.bera.test
library(forecast)    # Acf(), Pacf(), auto.arima() — nicer than base R acf()

# ── 1. Load ───────────────────────────────────────────────────────────────────
raw <- read_csv("PCEPI.csv") |>
  rename(date = observation_date, Y = PCEPI) |>
  mutate(date = as.Date(date))

# Quick sanity check
glimpse(raw)          # should be 806 rows, 2 cols
range(raw$date)       # 1959-01-01 to 2026-02-01

# ── 2. Build time series objects ──────────────────────────────────────────────
# ts() needs a start and frequency — monthly = 12
Y <- ts(raw$Y,
        start     = c(1959, 1),
        frequency = 12)

# ── 3. Transformation: log-differences × 100 ──────────────────────────────────
# Why log? Stabilises variance + gives % interpretation
# Why ×100? Converts log-ratio to approximate % change (monthly inflation)
# diff(log(Y)) drops the first obs → n = 805

y_logdiff <- diff(log(Y)) * 100
y_pctdiff <- (diff(Y) / stats::lag(Y, -1)) * 100

# ── 4. Attach dates to a tidy dataframe (useful for ggplot later) ─────────────
df <- tibble(
  date = raw$date,
  Y    = as.numeric(Y)
) |>
  mutate(
    y_logdiff = (log(Y) - log(lag(Y))) * 100,   # log-difference × 100
    y_pctdiff = (Y / lag(Y) - 1) * 100          # simple % change
  )

# ── 5. Quick plots ────────────────────────────────────────────────────────────
df |>
  select(date, y_logdiff, y_pctdiff) |>
  drop_na() |>
  pivot_longer(-date) |>
  ggplot(aes(x = date, y = value, colour = name)) +
  geom_line(alpha = 0.8) +
  labs(title = "Log-diff vs % diff",
       y = "% change", x = "") +
  theme_minimal()

par(mfrow = c(3, 1))
plot(Y,            main = "Level (Y)",              ylab = "Index")
plot(ts(df$y_logdiff[-1], start = c(1959,2), frequency = 12),
     main = "Log-diff × 100 (y)",     ylab = "% change")
plot(ts(df$y_pctdiff[-1], start = c(1959,2), frequency = 12),
     main = "Pct-diff (y_pct)",        ylab = "% change")
