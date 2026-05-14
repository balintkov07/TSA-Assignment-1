# 0. Libraries ----
library(tidyverse)
library(zoo)
library(tseries)
library(forecast)
library(lmtest)
library(ggplot2)

# 1. Load ----
raw <- read_csv("PCEPI.csv") |>
  rename(date = observation_date, Y = PCEPI) |>
  mutate(date = as.Date(date))

glimpse(raw)        # 806 rows, 2 cols
range(raw$date)     # 1959-01-01 to 2026-02-01

# 2. Time series object ----
Y <- ts(raw$Y,
        start     = c(1959, 1),
        frequency = 12)

# 3. Transformations ----

## 3a. Tidy dataframe with log-diff and pct-diff
df <- tibble(
  date = raw$date,
  Y    = as.numeric(Y)
) |>
  mutate(
    y_logdiff = (log(Y) - log(lag(Y))) * 100,  # log-difference × 100
    y_pctdiff = (Y / lag(Y) - 1) * 100         # simple % change
  )

## 3b. ts versions for time series functions
y_logdiff <- diff(log(Y)) * 100
y_pctdiff <- (diff(Y) / stats::lag(Y, -1)) * 100

# 4. De-trending ----

## 4a. Deterministic: OLS on time index
t <- 1:length(Y)
trend_model    <- lm(Y ~ t)
Y_det_detrend  <- ts(residuals(trend_model),
                     start     = c(1959, 1),
                     frequency = 12)

## 4b. Stochastic: first difference (already in y_logdiff)
# y_logdiff IS the stochastic de-trending of log(Y)

# 5. Plots ----

## 5a. Level + transformations
par(mfrow = c(3, 1))
plot(Y,
     main = "Level (Y)", ylab = "Index")
plot(ts(df$y_logdiff[-1], start = c(1959, 2), frequency = 12),
     main = "Log-diff × 100", ylab = "% change")
plot(ts(df$y_pctdiff[-1], start = c(1959, 2), frequency = 12),
     main = "Pct-diff", ylab = "% change")

## 5b. De-trending comparison
par(mfrow = c(3, 1))
plot(Y,
     main = "Level (Y) — raw", ylab = "Index")
plot(Y_det_detrend,
     main = "Deterministic de-trending (OLS residuals)", ylab = "Residual")
abline(h = 0, col = "grey50", lty = 2)
plot(y_logdiff,
     main = "Stochastic de-trending (log first-difference)", ylab = "% change")
abline(h = 0, col = "grey50", lty = 2)

avg_change <- mean(y_logdiff)
abline(h = avg_change, col = "red", lty = 2)

# 6. AR(1) on log-differences ----

## 6a. Estimate
ar1 <- arima(y_logdiff, order = c(1, 0, 0))
summary(ar1)

## 6b. Extract key output
coeftest(ar1)          # t-stats and p-values — needs library(lmtest)
ar1$coef               # phi_1 and intercept
ar1$sigma2             # residual variance

## 6c. Test 2 transformations
# Step 1 — test log(Y): does the log-level need differencing?
adf.test(log(Y))        # expect FAIL to reject H0 → unit root present → difference

# Step 2 — test y_logdiff: did one difference fix it?
adf.test(y_logdiff)     # expect REJECT H0 → now stationary → stop here

# 8.  ACF and PACF
## ── Compute ACF and PACF (no plot) ────────────────────────────────────────────
acf_vals  <- acf(Y,  lag.max = 20, plot = FALSE)
pacf_vals <- pacf(Y, lag.max = 20, plot = FALSE)

## ── Combine into a clean data frame ───────────────────────────────────────────
results <- data.frame(
  Lag  = 1:20,
  ACF  = round(acf_vals$acf[2:21],  4),   # lag 0 (=1) is excluded
  PACF = round(pacf_vals$acf[1:20], 4)
)

## ── Display as a formatted table ──────────────────────────────────────────────
## ── Compute ACF and PACF (no plot) ────────────────────────────────────────────
acf_vals  <- acf(Y,  lag.max = 20, plot = FALSE)
pacf_vals <- pacf(Y, lag.max = 20, plot = FALSE)

## ── Combine into a clean data frame ───────────────────────────────────────────
results <- data.frame(
  Lag  = 1:20,
  ACF  = round(acf_vals$acf[2:21],  4),   # lag 0 (=1) is excluded
  PACF = round(pacf_vals$acf[1:20], 4)
)

## ── Display as a formatted table ──────────────────────────────────────────────
print(results, row.names = FALSE)