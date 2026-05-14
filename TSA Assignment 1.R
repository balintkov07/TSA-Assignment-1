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

# 8. ACF and PACF ----

## 8a. Helper: compute + return clean data frame --------------------------------
acf_table <- function(series, lag.max = 20) {
  data.frame(
    Lag  = 1:lag.max,
    ACF  = round(acf(series,  lag.max = lag.max, plot = FALSE)$acf[2:(lag.max+1)], 4),
    PACF = round(pacf(series, lag.max = lag.max, plot = FALSE)$acf[1:lag.max],     4)
  )
}

results_Y <- acf_table(Y)
results_y <- acf_table(y_logdiff)

## 8b. Print tables -------------------------------------------------------------
cat("── ACF / PACF: Level series Y ──\n");  print(results_Y, row.names = FALSE)
cat("── ACF / PACF: Log-differences y ──\n"); print(results_y, row.names = FALSE)

## 8c. Plot: level vs differences side by side ----------------------------------
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

acf(Y,          lag.max = 20, main = "ACF — Level (Y)",        xlab = "Lag (months)")
pacf(Y,         lag.max = 20, main = "PACF — Level (Y)",       xlab = "Lag (months)")
acf(y_logdiff,  lag.max = 20, main = "ACF — Log-diff y (Q4)",  xlab = "Lag (months)")
pacf(y_logdiff, lag.max = 20, main = "PACF — Log-diff y (Q4)", xlab = "Lag (months)")

# 5. Model selection Question 5 ----


## 5a. Define candidates (motivated by Q4 ACF/PACF reading) --------------------
candidates <- list(
  "AR(1)"     = c(1, 0, 0),
  "AR(2)"     = c(2, 0, 0),
  "AR(3)"     = c(3, 0, 0),
  "MA(1)"     = c(0, 0, 1),
  "MA(2)"     = c(0, 0, 2),
  "MA(3)"     = c(0, 0, 3),
  "ARMA(1,1)" = c(1, 0, 1),
  "ARMA(1,2)" = c(1, 0, 2),
  "ARMA(2,1)" = c(2, 0, 1)
)

# Fits for y

## 5b. Fit all, collect AIC / BIC -----------------------------------------------
fits <- lapply(candidates, \(ord) arima(y_logdiff, order = ord))

ic_table <- data.frame(
  Model  = names(candidates),
  k      = sapply(fits, \(m) length(m$coef)),       # number of parameters
  LogLik = round(sapply(fits, \(m) as.numeric(logLik(m))), 2),
  AIC    = round(sapply(fits, AIC), 2),
  BIC    = round(sapply(fits, BIC), 2)
) |>
  arrange(AIC)

print(ic_table, row.names = FALSE)

## 5c. Best model by AIC --------------------------------------------------------
best_name  <- ic_table$Model[1]
best_model <- fits[[best_name]]

cat("\n── Best model (AIC):", best_name, "──\n")
print(best_model)

## 5d. Additional argument — BIC (Q5b) ------------------------------------------
# BIC penalises parameters more heavily than AIC
# If AIC and BIC agree → stronger case for the chosen model
# If they disagree → prefer the simpler model (parsimony principle)

cat("\n── Best model (BIC):", ic_table$Model[which.min(ic_table$BIC)], "──\n")
cat("── AIC and BIC agree?", ic_table$Model[1] == ic_table$Model[which.min(ic_table$BIC)], "──\n")

## 5e. Quick residual check on top 3 -------------------------------------------
top3 <- ic_table$Model[1:3]

par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
for (nm in top3) {
  acf(residuals(fits[[nm]]), lag.max = 20,
      main = paste("Residual ACF —", nm), xlab = "Lag (months)")
}

# 5b. Model selection — Level Y (comparison, force ML) ----
fits_Y <- lapply(candidates, \(ord) {
  tryCatch(
    arima(Y, order = ord, method = "ML"),
    error = \(e) NULL          # some models may still fail — skip them
  )
})

# Remove any NULLs (failed fits)
fits_Y <- Filter(Negate(is.null), fits_Y)

ic_table_Y <- data.frame(
  Model  = names(fits_Y),
  k      = sapply(fits_Y, \(m) length(m$coef)),
  LogLik = round(sapply(fits_Y, \(m) as.numeric(logLik(m))), 2),
  AIC    = round(sapply(fits_Y, AIC), 2),
  BIC    = round(sapply(fits_Y, BIC), 2)
) |>
  arrange(AIC)

print(ic_table_Y, row.names = FALSE)
cat("\n── Best by AIC:", ic_table_Y$Model[1], "──\n")
print(fits_Y[[ic_table_Y$Model[1]]])

# 6. Residual diagnostics ----

## 6a. Extract residuals from best model (ARMA(1,2)) ----------------------------
resid_best <- residuals(best_model)

## 6b. Visual: residuals over time ---------------------------------------------
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot(resid_best,
     main = "Residuals over time",
     ylab = "Residual", xlab = "")
abline(h = 0, col = "grey50", lty = 2)

## 6c. ACF and PACF of residuals ------------------------------------------------
acf(resid_best,
    lag.max = 20, main = "ACF — residuals",
    xlab = "Lag (months)")

pacf(resid_best,
     lag.max = 20, main = "PACF — residuals",
     xlab = "Lag (months)")

## 6d. QQ plot ------------------------------------------------------------------
qqnorm(resid_best, main = "QQ plot — residuals")
qqline(resid_best, col = "steelblue")

## 6e. Formal white noise: Ljung-Box test ---------------------------------------
# H0: residuals are white noise up to lag K
# Rule of thumb: K = min(10, n/5) for monthly data
lb_test <- Box.test(resid_best, lag = 12, type = "Ljung-Box")
cat("── Ljung-Box (lag 12) ──\n")
print(lb_test)

## 6f. Concrete improvement: check for ARCH effects ----------------------------
# ACF of squared residuals — if significant → variance clustering → GARCH needed

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

acf(resid_best^2,
    lag.max = 20, main = "ACF — squared residuals",
    xlab = "Lag (months)")

pacf(resid_best^2,
     lag.max = 20, main = "PACF — squared residuals",
     xlab = "Lag (months)")

# Formal ARCH test (requires FinTS package)
install.packages("FinTS")
library(FinTS)
arch_test <- ArchTest(resid_best, lags = 12)
cat("\n── ARCH LM test (lag 12) ──\n")
print(arch_test)
# FinTS masked forecast::Acf — restore it
detach("package:FinTS", unload = TRUE)

# 7. Forecasting and MSFE evaluation ----

## 7a. Setup -------------------------------------------------------------------
h_max   <- 4
y_vec   <- as.numeric(y_logdiff)    # strip ts class for clean indexing
y_times <- time(y_logdiff)

# Evaluation sample: actuals in Jan 2016 – Feb 2026
eval_idx <- which(
  y_times >= (2016 + 0/12) - 1e-6 &
    y_times <= (2026 + 1/12) + 1e-6
)

cat("Evaluation months:", length(eval_idx),
    "| From Jan 2016 to Feb 2026\n")

## 7b. Recursive MSFE function -------------------------------------------------
# For horizon h: fit on data up to (eval_idx - h), forecast h steps,
# compare to actual at eval_idx
recursive_msfe <- function(series, order, eval_idx, h_max,
                           min_train = 60) {
  sapply(1:h_max, function(h) {
    
    origins <- eval_idx - h          # forecast origin for each eval point
    valid   <- origins >= min_train  # need enough history to estimate
    
    sq_err <- mapply(function(t, t_act) {
      fit <- arima(series[1:t], order = order)
      fc  <- as.numeric(predict(fit, n.ahead = h)$pred)[h]
      (series[t_act] - fc)^2
    }, origins[valid], eval_idx[valid])
    
    mean(sq_err)
  })
}

## 7c. Run forecasts (this takes a few minutes — ~1000 model fits) -------------
cat("Running AR(1) forecasts...\n")
msfe_ar1    <- recursive_msfe(y_vec, c(1, 0, 0), eval_idx, h_max)

cat("Running ARMA(1,2) forecasts...\n")
msfe_arma12 <- recursive_msfe(y_vec, c(1, 0, 2), eval_idx, h_max)

## 7d. MSFE comparison table ---------------------------------------------------
msfe_table <- data.frame(
  Horizon    = 1:h_max,
  AR1        = round(msfe_ar1,    6),
  ARMA12     = round(msfe_arma12, 6),
  Ratio      = round(msfe_arma12 / msfe_ar1, 4)  # < 1: ARMA wins; > 1: AR wins
)

cat("\n── MSFE by horizon ──\n")
print(msfe_table, row.names = FALSE)

## 7e. Plot MSFE by horizon ----------------------------------------------------
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))

ylim <- range(c(msfe_ar1, msfe_arma12)) * c(0.95, 1.05)

plot(1:h_max, msfe_ar1,
     type = "b", pch = 16, col = "steelblue", lwd = 1.5,
     ylim = ylim, xaxt = "n",
     xlab = "Horizon (months ahead)", ylab = "MSFE",
     main = "MSFE by forecast horizon: AR(1) vs ARMA(1,2)")
axis(1, at = 1:h_max)
lines(1:h_max, msfe_arma12,
      type = "b", pch = 17, col = "coral", lwd = 1.5)
legend("topleft",
       legend = c("AR(1)", "ARMA(1,2)"),
       col    = c("steelblue", "coral"),
       pch    = c(16, 17), lty = 1, bty = "n")

## 7f. Shorter evaluation period (Q7c) -----------------------------------------
eval_idx_short <- which(
  y_times >= (2020 + 0/12) - 1e-6 &
    y_times <= (2026 + 1/12) + 1e-6
)

cat("\nShorter eval (Jan 2020 – Feb 2026):", length(eval_idx_short), "months\n")

msfe_ar1_short    <- recursive_msfe(y_vec, c(1, 0, 0), eval_idx_short, h_max)
msfe_arma12_short <- recursive_msfe(y_vec, c(1, 0, 2), eval_idx_short, h_max)

msfe_short <- data.frame(
  Horizon       = 1:h_max,
  AR1_short     = round(msfe_ar1_short,    6),
  ARMA12_short  = round(msfe_arma12_short, 6),
  Ratio_short   = round(msfe_arma12_short / msfe_ar1_short, 4)
)

cat("\n── MSFE — shorter evaluation period ──\n")
print(msfe_short, row.names = FALSE)
