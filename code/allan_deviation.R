# allan_deviation.R
#
# Overlapping Allan deviation computation and plotting for LC-MS rate
# time series. Used to characterize the stability regime of the
# detector chain and to demonstrate that the partitioning protocols
# tested operate within the white-noise branch of the Allan variance
# spectrum (radiometer equation regime).
#
# References:
#   Allan, D.W. (1966). Statistics of atomic frequency standards.
#     Proceedings of the IEEE, 54(2), 221-230.
#   Schieder, R. and Kramer, C. (2001). Optimization of heterodyne
#     observations using Allan variance measurements.
#     Astronomy and Astrophysics, 373, 746-756.
#
# Dependencies: dplyr, tibble, ggplot2.


# Core computation -------------------------------------------------------------

#' Compute overlapping Allan deviation for a uniformly sampled time series.
#'
#' For a sequence of measurements y_1, ..., y_N taken at uniform intervals
#' tau_0, returns the overlapping Allan deviation at a set of averaging
#' factors m. The averaging time at factor m is tau = m * tau_0.
#'
#' The overlapping estimator is:
#'   sigma_A^2(tau) = (1 / (2 * (N - 2m + 1))) *
#'                    sum_{i=1}^{N - 2m + 1} (ybar_{i+m} - ybar_i)^2
#' where ybar_i is the mean of y_{i}, ..., y_{i+m-1}.
#'
#' @param y Numeric vector of measurements.
#' @param tau_0 Sampling interval (seconds).
#' @param m_values Integer vector of averaging factors. Defaults to roughly
#'   log-spaced factors from 1 to floor(N / 5).
#' @return A tibble with columns m, tau, allan_deviation, n_pairs.

compute_allan_deviation <- function(y, tau_0, m_values = NULL) {
  n <- length(y)
  if (n < 5L) {
    stop("Need at least 5 samples to compute Allan deviation.")
  }

  if (is.null(m_values)) {
    max_m <- as.integer(floor(n / 5))
    m_values <- unique(as.integer(round(
      2 ^ seq(0, log2(max_m), length.out = 25L)
    )))
    m_values <- m_values[m_values >= 1L & m_values <= max_m]
  }

  m_values <- m_values[(2L * m_values + 1L) <= n]
  if (length(m_values) == 0L) {
    stop("No valid averaging factors for this time series length.")
  }

  # Cumulative sum allows O(1) computation of any window mean.
  cs <- c(0, cumsum(y))

  allan_var <- vapply(m_values, function(m) {
    y_bar <- (cs[(m + 1L):(n + 1L)] - cs[1L:(n - m + 1L)]) / m
    n_pairs <- n - 2L * m + 1L
    diffs <- y_bar[(m + 1L):(m + n_pairs)] - y_bar[1L:n_pairs]
    mean(diffs ^ 2) / 2
  }, numeric(1))

  n_pairs <- vapply(m_values, function(m) {
    as.integer(n - 2L * m + 1L)
  }, integer(1))

  tibble::tibble(
    m = m_values,
    tau = m_values * tau_0,
    allan_variance = allan_var,
    allan_deviation = sqrt(allan_var),
    n_pairs = n_pairs,
    n_samples = as.integer(n)
  )
}


# Confidence intervals ---------------------------------------------------------

#' Chi-squared confidence interval for Allan deviation.
#'
#' Uses conservative non-overlapping equivalent degrees of freedom:
#' edf = floor(N / (2m)) - 1 per replicate, pooled across n_reps
#' independent replicates. This underestimates the true edf for
#' overlapping estimators, giving wider (conservative) intervals.
#'
#' @param allan_variance Numeric vector of mean Allan variance estimates.
#' @param n_samples Integer, length of the original time series.
#' @param m Integer vector of averaging factors.
#' @param n_reps Number of independent replicates averaged.
#' @param alpha Significance level (default 0.05 for 95% CI).
#' @return Tibble with edf, ci_lower, ci_upper (on deviation scale).

compute_allan_ci <- function(allan_variance, n_samples, m,
                             n_reps = 1L, alpha = 0.05) {
  edf_single <- pmax(1L, as.integer(floor(n_samples / (2L * m)) - 1L))
  edf_total  <- n_reps * edf_single

  ci_lower_var <- edf_total * allan_variance / stats::qchisq(1 - alpha/2, edf_total)
  ci_upper_var <- edf_total * allan_variance / stats::qchisq(alpha/2, edf_total)

  tibble::tibble(
    edf = edf_total,
    ci_lower = sqrt(ci_lower_var),
    ci_upper = sqrt(ci_upper_var)
  )
}


# Multi-group wrapper ----------------------------------------------------------

#' Compute Allan deviation for grouped data.
#'
#' Applies compute_allan_deviation() to each group defined by group_cols.
#' The within-group time series is assumed uniformly sampled at tau_0.
#' Caller is responsible for trimming plateau edge effects before calling.
#'
#' @param data Data frame containing the time series and group columns.
#' @param value_col Name of the column with the rate measurements.
#' @param tau_0 Sampling interval (seconds).
#' @param group_cols Character vector of column names defining groups,
#'   e.g. c("concentration", "replicate").
#' @param m_values Optional shared m_values vector.
#' @return Tibble with group columns plus (m, tau, allan_deviation, n_pairs).

compute_allan_deviation_grouped <- function(
  data,
  value_col,
  tau_0,
  group_cols,
  m_values = NULL
) {
  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(function(group_data, group_keys) {
      compute_allan_deviation(
        y = group_data[[value_col]],
        tau_0 = tau_0,
        m_values = m_values
      )
    }) |>
    dplyr::ungroup()
}


# Plotting ---------------------------------------------------------------------

#' Plot Allan deviation curves on log-log axes with an optional white-noise
#' reference line of slope -1/2.
#'
#' @param allan_data Tibble with columns tau, allan_deviation, and optionally
#'   a grouping column for color.
#' @param color_col Optional column name to map to color and group.
#' @param white_noise_reference Logical; draw a slope -1/2 reference line.
#' @param reference_anchor Optional list(tau = ..., sigma = ...) to anchor
#'   the reference line. Defaults to the median of the data.
#' @return A ggplot object.

plot_allan_deviation <- function(
  allan_data,
  color_col = NULL,
  white_noise_reference = TRUE,
  reference_anchor = NULL
) {
  base_aes <- if (is.null(color_col)) {
    ggplot2::aes(x = tau, y = allan_deviation)
  } else {
    ggplot2::aes(
      x = tau,
      y = allan_deviation,
      color = .data[[color_col]],
      group = .data[[color_col]]
    )
  }

  p <- ggplot2::ggplot(allan_data, base_aes) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = expression(paste("Averaging time ", tau, " (s)")),
      y = expression(paste("Allan deviation ", sigma[A](tau)))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )

  if (white_noise_reference) {
    anchor <- reference_anchor %||% list(
      tau = stats::median(allan_data$tau),
      sigma = stats::median(allan_data$allan_deviation)
    )
    # On log10 axes, slope -1/2 line: log10(y) = -0.5 * log10(tau) + c
    # where c = log10(sigma_anchor) + 0.5 * log10(tau_anchor).
    log_intercept <- log10(anchor$sigma) + 0.5 * log10(anchor$tau)
    p <- p + ggplot2::geom_abline(
      slope = -0.5,
      intercept = log_intercept,
      linetype = "dashed",
      color = "grey40"
    )
  }

  p
}

# Null-coalescing operator (ggplot2/dplyr already provide one in recent
# versions; defined here for self-containment).
`%||%` <- function(a, b) if (is.null(a)) b else a


# Demo / sanity check ----------------------------------------------------------

#' Generate a simulated rate time series for testing.
#'
#' Returns a tibble with three concentration groups, each with Poisson
#' counting noise plus a proportional component matching the chemistry
#' paper's two-component model (no third in-band component). The Allan
#' deviation of these data should follow slope -1/2 on log-log axes for
#' all tau within the experimental window.
#'
#' @param n_samples Number of cycles per concentration.
#' @param tau_0 Cycle time in seconds.
#' @param cv_proportional Proportional noise CV (slow common-mode part).
#' @param seed Random seed.
#' @return Tibble with concentration_label, sample_index, rate.

simulate_demo_data <- function(
  n_samples = 1000L,
  tau_0 = 0.1,
  cv_proportional = 0.02,
  seed = 42L
) {
  set.seed(seed)
  concentrations <- c(low = 50, mid = 500, high = 1500)

  rows <- lapply(names(concentrations), function(label) {
    lambda <- concentrations[[label]]
    proportional_factor <- 1 + stats::rnorm(n_samples, sd = cv_proportional)
    rates <- stats::rpois(n_samples, lambda) * proportional_factor

    tibble::tibble(
      concentration_label = label,
      lambda = lambda,
      sample_index = seq_len(n_samples),
      rate = rates
    )
  })

  do.call(rbind, rows)
}


# Example usage ----------------------------------------------------------------

# To run the demo and verify the white-noise reference line lies on top
# of the simulated curves:
#
#   source("R/allan_deviation.R")
#   sim <- simulate_demo_data(n_samples = 2000L)
#   allan_results <- compute_allan_deviation_grouped(
#     data = sim,
#     value_col = "rate",
#     tau_0 = 0.1,
#     group_cols = "concentration_label"
#   )
#   p <- plot_allan_deviation(allan_results, color_col = "concentration_label")
#   print(p)
#
# For the actual paper figure, replace simulate_demo_data() with the
# function that loads the plateau rate time series from the trimmed
# flow-injection data (one rate per cycle per concentration), and
# pass tau_0 = 0.1 (100 ms cycle time).