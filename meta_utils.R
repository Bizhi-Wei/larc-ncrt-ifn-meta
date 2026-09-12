# Shared effect-size and random-effects meta-analysis helpers.
# Primary inference: REML heterogeneity + modified Hartung-Knapp confidence interval.

hedges_g <- function(x, group, good = "good", poor = "poor") {
  keep <- is.finite(x) & group %in% c(good, poor)
  x <- as.numeric(x[keep])
  group <- as.character(group[keep])
  x_good <- x[group == good]
  x_poor <- x[group == poor]
  n_good <- length(x_good)
  n_poor <- length(x_poor)
  if (n_good < 2L || n_poor < 2L) stop("Each group must contain at least two observations")

  pooled_sd <- sqrt(
    ((n_good - 1) * stats::var(x_good) + (n_poor - 1) * stats::var(x_poor)) /
      (n_good + n_poor - 2)
  )
  if (!is.finite(pooled_sd) || pooled_sd <= 0) stop("Pooled standard deviation is not positive")

  d <- (mean(x_good) - mean(x_poor)) / pooled_sd
  correction <- 1 - 3 / (4 * (n_good + n_poor) - 9)
  estimate <- correction * d
  se <- sqrt(
    (n_good + n_poor) / (n_good * n_poor) +
      estimate^2 / (2 * (n_good + n_poor - 2))
  )
  c(g = estimate, se = se, n_good = n_good, n_poor = n_poor)
}

pool_random_effects <- function(estimate, se, prediction = TRUE) {
  keep <- is.finite(estimate) & is.finite(se) & se > 0
  estimate <- as.numeric(estimate[keep])
  se <- as.numeric(se[keep])
  k <- length(estimate)
  # Prediction intervals use t with df = k - 2 and are undefined at k = 2.
  # Exploratory k = 2 subgroups may still request pooled estimates (prediction = FALSE).
  if (k < 2L) stop("At least two studies are required")
  if (prediction && k < 3L) stop("Prediction intervals require at least three studies")

  variance <- se^2
  restricted_nll <- function(tau2) {
    weight <- 1 / (variance + tau2)
    mu <- sum(weight * estimate) / sum(weight)
    0.5 * (
      sum(log(variance + tau2)) + log(sum(weight)) +
        sum(weight * (estimate - mu)^2)
    )
  }
  upper <- max(1, 100 * stats::var(estimate), 100 * max(variance))
  optimum <- stats::optimize(restricted_nll, interval = c(0, upper))
  tau2 <- if (restricted_nll(optimum$minimum) < restricted_nll(0)) optimum$minimum else 0

  weight <- 1 / (variance + tau2)
  mu <- sum(weight * estimate) / sum(weight)
  hk_q <- sum(weight * (estimate - mu)^2) / (k - 1)
  hk_scale <- max(1, hk_q)
  se_mu <- sqrt(hk_scale / sum(weight))
  df <- k - 1
  critical <- stats::qt(0.975, df = df)
  ci_lo <- mu - critical * se_mu
  ci_hi <- mu + critical * se_mu
  p_value <- 2 * stats::pt(-abs(mu / se_mu), df = df)

  fixed_weight <- 1 / variance
  fixed_mu <- sum(fixed_weight * estimate) / sum(fixed_weight)
  Q <- sum(fixed_weight * (estimate - fixed_mu)^2)
  I2 <- max(0, (Q - df) / max(Q, .Machine$double.eps)) * 100
  Q_p <- stats::pchisq(Q, df = df, lower.tail = FALSE)

  pi_lo <- pi_hi <- NA_real_
  if (prediction) {
    pi_critical <- stats::qt(0.975, df = k - 2)
    pi_half <- pi_critical * sqrt(tau2 + se_mu^2)
    pi_lo <- mu - pi_half
    pi_hi <- mu + pi_half
  }

  list(
    g = mu, se = se_mu, lo = ci_lo, hi = ci_hi, p = p_value,
    I2 = I2, Q = Q, Qp = Q_p, tau2 = tau2,
    pi_lo = pi_lo, pi_hi = pi_hi, k = k,
    method = "REML + modified Hartung-Knapp"
  )
}

pool_dl <- function(estimate, se) {
  keep <- is.finite(estimate) & is.finite(se) & se > 0
  estimate <- as.numeric(estimate[keep])
  se <- as.numeric(se[keep])
  k <- length(estimate)
  if (k < 3L) stop("At least three studies are required")

  weight <- 1 / se^2
  fixed_mu <- sum(weight * estimate) / sum(weight)
  Q <- sum(weight * (estimate - fixed_mu)^2)
  df <- k - 1
  denominator <- sum(weight) - sum(weight^2) / sum(weight)
  tau2 <- max(0, (Q - df) / denominator)
  random_weight <- 1 / (se^2 + tau2)
  mu <- sum(random_weight * estimate) / sum(random_weight)
  se_mu <- sqrt(1 / sum(random_weight))
  list(
    g = mu, se = se_mu, lo = mu - 1.96 * se_mu, hi = mu + 1.96 * se_mu,
    p = 2 * stats::pnorm(-abs(mu / se_mu)),
    I2 = max(0, (Q - df) / max(Q, .Machine$double.eps)) * 100,
    Q = Q, Qp = stats::pchisq(Q, df = df, lower.tail = FALSE), tau2 = tau2,
    k = k, method = "DerSimonian-Laird"
  )
}
