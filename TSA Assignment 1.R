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

y <- diff(log(Y)) * 100

# ── 4. Attach dates to a tidy dataframe (useful for ggplot later) ─────────────
df <- tibble(
  date = raw$date,
  Y    = as.numeric(Y)
) |>
  mutate(
    y = c(NA, diff(log(Y)) * 100)   # NA for first row — no lag available
  )

# ── 5. Quick plots ────────────────────────────────────────────────────────────
par(mfrow = c(2, 1))

plot(Y,
     main = "PCE Price Index — Level (Y)",
     ylab = "Index (2017 = 100)", xlab = "")

plot(y,
     main = "Log-Differences × 100 — Monthly Inflation (y)",
     ylab = "% change", xlab = "")
abline(h = 0, col = "grey50", lty = 2)