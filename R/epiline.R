# EpiLine (github.com/BDI-pathogens/EpiLine), fitted to the cases only, for
# the comparison of nowcasts (R/comparison.R). Its model: the number of
# people with onset on each day grows at a rate which follows a random walk
# (in steps of two days); the delay from onset to report is Johnson SU, its
# four parameters each following a random walk over onset days; the cases
# reported each day are negative binomial about the onsets convolved with
# the delay; and a line list of (onset, report) pairs informs the delay,
# allowing for its truncation by the days of report seen. Its Stan model is
# taken from the package's repository at a fixed commit, and fitted with
# cmdstanr.
#
# The reports here are snapshots of the counts by onset day, the first (12
# July) holding all the cases reported until then, and later ones 1 to 6 days
# apart. The days of report fitted start the day after the first report, so
# the cases in it are left out (the model allows for the days of report
# starting after the onsets: it only loses what they tell); each later
# report's new cases of each onset day (from the counts made monotone,
# comparison_monotone) are line-list pairs on the day of the report, and its
# new cases in all are spread evenly over the days since the report before.

epiline_commit <- "4f5b938a106df2d5a2104d2981c4058fc36a3615"

#' The settings of EpiLine: those of its function symptom_report.fit, but the
#' longest delay from onset to report (t_symptom_pre) 60 days, as for
#' epinowcast
epiline_settings <- function() {
  list(t_symptom_pre = 60L, t_symptom_post = 5L, chains = 4L, iter_warmup = 1000L,
       iter_sampling = 1000L,
       priors = list(prior_gamma_min = -20, prior_delta_min = -20, prior_xi_min = -20,
                     prior_lambda_min = 0.1, prior_gamma_max = 20, prior_delta_max = 20,
                     prior_xi_max = 20, prior_lambda_max = 20, prior_xi_gp_sd_max = 0.6,
                     prior_lambda_gp_sd_max = 0.4, prior_gamma_gp_sd_max = 0.4,
                     prior_delta_gp_sd_max = 0.2, prior_r_0_min = -0.1, prior_r_0_max = 0.2,
                     prior_r_gp_sd_max = 0.055, prior_phi_od_max = 1,
                     prior_log_symptoms0_min = log(0.5), prior_log_symptoms0_max = log(50),
                     hyper_gp_period_r = 2L, hyper_gp_period_dist = 2L))
}

#' The data of a fit
#'
#' @param snapshots as build_bvd_snapshots
#' @param date_now the fit date
#' @param s as epiline_settings
#' @return a list: 'stan' (the data for the Stan model), 'start' (the first
#'   day of report fitted, day 1), 'last' (the last report fitted) and
#'   'delays' (the delays of the line list, one per case)
epiline_data <- function(snapshots, date_now, s = epiline_settings()) {
  date_now <- as.Date(date_now)
  counts <- comparison_monotone(snapshots, "cases_reported", date_now)
  reports <- sort(unique(counts$date_reported))
  start <- reports[1] + 1
  last <- max(reports)
  new <- counts %>%
    dplyr::group_by(.data$date_onset) %>%
    dplyr::mutate(n = .data$count - dplyr::lag(.data$count, default = 0)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$date_reported >= start, .data$n > 0)
  days <- seq(start, last, by = "day")
  # each report's new cases spread over the days since the report before
  totals <- new %>% dplyr::group_by(.data$date_reported) %>% dplyr::summarise(n = sum(.data$n))
  reported <- integer(length(days))
  for(k in seq_len(nrow(totals))) {
    d <- totals$date_reported[k]
    before <- max(reports[reports < d])
    span <- match(seq(before + 1, d, by = "day"), days)
    n <- totals$n[k]
    reported[span] <- reported[span] + (n %/% length(span)) +
      as.integer(seq_along(span) > length(span) - n %% length(span))
  }
  rep_day <- as.integer(new$date_reported - start) + 1L
  onset_day <- as.integer(new$date_onset - start) + 1L
  onset_day <- pmin(pmax(onset_day, rep_day - s$t_symptom_pre), rep_day + s$t_symptom_post)
  ll <- dplyr::tibble(report = rep_day, symptom = onset_day, N = new$n) %>%
    dplyr::group_by(.data$report, .data$symptom) %>%
    dplyr::summarise(N = sum(.data$N), .groups = "drop")
  stan <- c(list(t_rep = length(days), t_symptom_pre = s$t_symptom_pre,
                 t_symptom_post = s$t_symptom_post, reported = reported, n_ll = nrow(ll),
                 ll_report = ll$report, ll_symptoms = ll$symptom, ll_N = ll$N,
                 t_static_dist = 0L), s$priors)
  list(stan = stan, start = start, last = last, delays = rep(ll$report - ll$symptom, ll$N))
}

#' Starting values as symptom_report.fit makes them: the Johnson SU
#' parameters matching the moments of the delays of the line list (but gamma
#' -1 if the skewness is beyond reach), the growth rate 0.02 and the random
#' walks flat
#'
#' @param data as epiline_data
#' @param chains the number of chains
epiline_inits <- function(data, chains) {
  x <- data$delays
  m <- mean(x); v <- stats::var(x); sk <- mean((x - m)^3) / v^1.5
  p <- data$stan
  f_var <- function(lambda, gamma, delta) {
    omega <- exp(1 / delta^2)
    0.5 * lambda^2 * (omega - 1) * (omega * cosh(2 * gamma / delta) + 1)
  }
  f_skew <- function(gamma, delta) {
    omega <- exp(1 / delta^2)
    -0.25 * sqrt(omega) * (omega - 1)^2 * (omega * (omega + 2) * sinh(3 * gamma / delta) +
                                             3 * sinh(gamma / delta)) / f_var(1, gamma, delta)^1.5
  }
  t_max <- p$t_rep + p$t_symptom_pre + p$t_symptom_post
  eps <- 1e-6
  lapply(seq_len(chains), function(k) {
    delta <- stats::runif(1, 1.5, 2)
    g <- function(gamma) f_skew(gamma, delta) - sk
    gamma <- if(g(p$prior_gamma_min) * g(p$prior_gamma_max) > 0) {
      if(sk > 0) p$prior_gamma_min else p$prior_gamma_max
    } else stats::uniroot(g, c(p$prior_gamma_min, p$prior_gamma_max))$root
    gamma <- min(max(gamma, p$prior_gamma_min + eps), p$prior_gamma_max - eps)
    # if no scale then matches the variance (the skewness is beyond reach, so
    # gamma is at a bound), a moderate right skew instead
    h <- function(l, gamma) f_var(l, gamma, delta) - v
    if(h(p$prior_lambda_min, gamma) * h(p$prior_lambda_max, gamma) > 0)
      gamma <- -1
    lambda <- stats::uniroot(h, c(p$prior_lambda_min, p$prior_lambda_max), gamma = gamma)$root
    omega <- exp(1 / delta^2)
    xi <- m + lambda * sqrt(omega) * sinh(gamma / delta)
    xi <- min(max(xi, p$prior_xi_min + eps), p$prior_xi_max - eps)
    n_dist <- ceiling(t_max / p$hyper_gp_period_dist)
    list(r_gp = rep(0, ceiling(t_max / p$hyper_gp_period_r)), r_0 = 0.02,
         gamma0 = gamma, delta0 = delta, lambda0 = lambda, xi0 = xi,
         xi_gp = rep(0, n_dist), lambda_gp = rep(0, n_dist),
         gamma_gp = rep(0, n_dist), delta_gp = rep(0, n_dist))
  })
}

#' EpiLine's Stan model, compiled
epiline_model <- function() {
  path <- file.path(tempdir(), "symptom_report_model.stan")
  if(!file.exists(path))
    utils::download.file(sprintf("https://raw.githubusercontent.com/BDI-pathogens/EpiLine/%s/inst/stan/symptom_report_model.stan",
                                 epiline_commit), path, quiet = TRUE)
  cmdstanr::cmdstan_model(path)
}

#' Fits EpiLine to the cases reported up to a fit date, and projects the
#' reports scored: for each draw, the new cases of each onset day reported
#' after the last report fitted and up to a report are negative binomial
#' (with the model's overdispersion) about the onsets of that day times the
#' probability of a delay in those days
#'
#' @param snapshots as build_bvd_snapshots
#' @param date_now the fit date
#' @param path where the projections are saved
#' @param settings as comparison_settings
#' @param s as epiline_settings
#' @param seed of the sampler and the projections
#' @return path; the file holds 'projections' (a list of one, for cases) and
#'   'diagnostics'
fit_epiline <- function(snapshots, date_now, path, settings = comparison_settings(),
                        s = epiline_settings(), seed = 1) {
  date_now <- as.Date(date_now)
  data <- epiline_data(snapshots, date_now, s)
  withr::local_seed(seed)
  start_time <- Sys.time()
  fit <- epiline_model()$sample(data = data$stan, init = epiline_inits(data, s$chains),
                                chains = s$chains, parallel_chains = s$chains,
                                iter_warmup = s$iter_warmup, iter_sampling = s$iter_sampling,
                                seed = seed, refresh = 0, show_messages = FALSE)
  minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  f <- function(name) unclass(posterior::as_draws_matrix(fit$draws(name)))
  symptoms <- f("symptoms"); xi <- f("xi"); lambda <- f("lambda"); gamma <- f("gamma")
  delta <- f("delta"); phi <- as.vector(f("phi_od"))
  p <- data$stan
  scored <- seq(date_now - settings$max_age, date_now, by = "day")
  # the column of an onset day: day 1 of report is column 1 + t_symptom_pre
  col <- as.integer(scored - data$start) + 1L + p$t_symptom_pre
  stopifnot(all(col >= 1), all(col <= ncol(symptoms)))
  counts <- comparison_monotone(snapshots, "cases_reported", date_now)
  now <- counts %>% dplyr::filter(.data$date_onset %in% scored) %>%
    dplyr::group_by(.data$date_onset) %>%
    dplyr::summarise(now = .data$count[which.max(.data$date_reported)])
  now <- dplyr::left_join(dplyr::tibble(date_onset = scored), now, by = "date_onset") %>%
    dplyr::mutate(now = dplyr::coalesce(.data$now, 0L))
  reports <- comparison_reports(snapshots, date_now, settings)
  djsu <- function(x, k) {
    z <- (x - xi[, k]) / lambda[, k]
    delta[, k] / lambda[, k] / sqrt(2 * pi) / sqrt(1 + z^2) * exp(-0.5 * (gamma[, k] + delta[, k] * asinh(z))^2)
  }
  draws <- lapply(reports$date_report, function(r) {
    m <- vapply(seq_along(scored), function(j) {
      tau <- as.integer(seq(data$last + 1, as.Date(r), by = "day") - scored[j])
      tau <- tau[tau >= -p$t_symptom_post & tau <= p$t_symptom_pre]
      prob <- if(length(tau) == 0) 0 else rowSums(vapply(tau, djsu, numeric(nrow(xi)), k = col[j]))
      mu <- symptoms[, col[j]] * prob
      now$now[j] + stats::rnbinom(length(mu), size = 1 / phi, mu = mu)
    }, numeric(nrow(symptoms)))
    comparison_thin(m, settings$n_draws)
  })
  names(draws) <- as.character(reports$date_report)
  summ <- fit$summary(c("symptoms", "xi", "lambda", "gamma", "delta"), "rhat")
  keep <- summ$variable %in% sprintf("%s[%d]", rep(c("symptoms", "xi", "lambda", "gamma", "delta"),
                                                   each = length(col)), col)
  sampler <- fit$diagnostic_summary(quiet = TRUE)
  out <- list(projections = list(list(model = "epiline", series = "cases", date_now = date_now,
                                      onsets = scored, draws = draws)),
              diagnostics = dplyr::tibble(model = "epiline", series = "cases", date_now = date_now,
                                          max_rhat = max(summ$rhat[keep], na.rm = TRUE),
                                          divergent = sum(sampler$num_divergent),
                                          at_max_treedepth = sum(sampler$num_max_treedepth),
                                          minutes = minutes))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, path)
  path
}
