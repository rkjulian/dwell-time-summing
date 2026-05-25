rm(list = ls())

# Exact 95% confidence intervals for CV via non-central t-distribution
# (Vangel 1996). Replaces the Monte Carlo simulation in cv_pdf_simulation.R.
#
# The statistic t = xbar*sqrt(n) / s has a non-central t distribution with
# df = n-1 and non-centrality parameter delta = mu*sqrt(n) / sigma = sqrt(n) / CV_true.
# Invert to find delta_L, delta_U such that the observed t falls at the
# alpha/2 and 1-alpha/2 quantiles, then CV = sqrt(n) / delta.

cv_ci <- function(cv_obs, n, alpha = 0.05) {
  df <- n - 1
  t_obs <- sqrt(n) / cv_obs

  # Find delta_upper: P(T < t_obs | delta_upper) = alpha/2
  # At delta_upper, t_obs is in the lower tail
  f_upper <- function(delta) pt(t_obs, df, ncp = delta) - alpha / 2
  delta_upper <- uniroot(f_upper, c(0.01, 1000))$root

  # Find delta_lower: P(T < t_obs | delta_lower) = 1 - alpha/2
  # At delta_lower, t_obs is in the upper tail
  f_lower <- function(delta) pt(t_obs, df, ncp = delta) - (1 - alpha / 2)
  delta_lower <- uniroot(f_lower, c(0.01, 1000))$root

  cv_lower <- sqrt(n) / delta_upper
  cv_upper <- sqrt(n) / delta_lower

  c(lower = cv_lower, upper = cv_upper)
}

# Data from Table 1, 50 ng/mL rows
n <- 3

# 96 ms dwell time: CV = 0.073
cv_96 <- 0.073
ci_96 <- cv_ci(cv_96, n)

# 6 ms dwell time: CV = 0.101
cv_06 <- 0.101
ci_06 <- cv_ci(cv_06, n)

cat("Exact 95% CI for CV (Vangel 1996, non-central t)\n")
cat("n =", n, "\n\n")
cat("96 ms (CV =", cv_96, "):",
    "[", round(ci_96["lower"], 3), ",", round(ci_96["upper"], 3), "]\n")
cat(" 6 ms (CV =", cv_06, "):",
    "[", round(ci_06["lower"], 3), ",", round(ci_06["upper"], 3), "]\n\n")
cat("Intervals overlap:",
    ci_06["lower"] < ci_96["upper"] && ci_96["lower"] < ci_06["upper"], "\n")
