# The data-driven model of cases and deaths (stan/nowcast_joint.stan): the
# delay from onset to a record is learned from the reports and may change
# smoothly over onset time, each report has an effect of its own, the new
# count of an onset day between reports is additions less removals, and the
# expected deaths of an onset day are its expected cases times a case
# fatality ratio which changes slowly. Each chain is fitted on its own (so
# that chains can run as separate jobs), its draws are saved, and the
# projections to any later day are made from the saved draws.

#' The new counts of each onset day between consecutive reports
#'
#' The first report of an onset day gives its count for the delays from 0
#' to that report (for the onset days the correction relabelled, the first
#' report of the corrected series); each later report the change since the
#' report before.
#'
#' @param snapshots as prepare_bvd_reports()$snapshots
#' @param column the name of the count
#' @return a tibble of t, lo, hi (the first and last delay of the interval),
#'   time_reported, n (the new count, which may be negative) and count
nowcast_triangle <- function(snapshots, column) {
  snapshots %>%
    dplyr::arrange(.data$time_onset, .data$time_reported) %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::mutate(
      count = .data[[column]],
      n = .data$count - dplyr::lag(.data$count, default = 0L),
      lo = dplyr::if_else(dplyr::row_number() == 1, 0L,
                          dplyr::lag(.data$time_reported) - .data$time_onset + 1L),
      hi = .data$time_reported - .data$time_onset) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(t = .data$time_onset, lo = .data$lo, hi = .data$hi,
                     time_reported = .data$time_reported, n = .data$n,
                     count = .data$count)
}

#' Knots of the basis in delay: quantiles of the delays at which records were added
#'
#' Each rise in a count between consecutive reports is taken as records added
#' at the delay of the later report, weighted by the size of the rise. First
#' counts of onset days made more than first_max days after onset are left
#' out, since the delays of their records are not known.
#'
#' @param tri as nowcast_triangle
#' @param probs the quantiles
#' @param first_max the longest delay of a first count which is used
#' @return the delays (in days) of the interior knots
nowcast_delay_knots <- function(tri, probs = c(0.25, 0.5, 0.75), first_max = 7) {
  added <- dplyr::filter(tri, .data$n > 0, !(.data$lo == 0 & .data$hi > first_max))
  o <- order(added$hi)
  cumulative <- cumsum(added$n[o]) / sum(added$n)
  knots <- vapply(probs, function(p) added$hi[o][which(cumulative >= p)[1]], numeric(1))
  knots <- unique(knots)
  knots[knots > 0 & knots < max(tri$hi)]
}

#' The log probabilities of the reference delay distribution
#'
#' A gamma distribution, discretised to whole days as the reporting delay of
#' the mechanistic model is (the probability of (d - 1/2, d + 1/2], from 0),
#' and renormalised over the delays 0 to D.
#'
#' @param D the longest delay
#' @param mean,sd of the gamma distribution, in days
#' @return a vector of log probabilities for the delays 0 to D
nowcast_reference_delay <- function(D, mean, sd) {
  log(nowcast_gamma_days(D, mean, sd))
}

#' A gamma distribution discretised to whole days 0 to D: the probability it
#' assigns to (d - 1/2, d + 1/2] (to [0, 1/2] for d = 0), renormalised
#' @inheritParams nowcast_reference_delay
nowcast_gamma_days <- function(D, mean, sd) {
  shape <- mean^2 / sd^2
  rate <- mean / sd^2
  p <- diff(stats::pgamma(c(0, seq_len(D + 1) - 0.5), shape, rate))
  p / sum(p)
}

#' The reference delay from onset to a death being recorded, as the
#' mechanistic model has it at the means of its priors: a fatal case counted
#' on its day of onset is recorded after the time from onset to death (gamma,
#' the probability of (z - 1, z] for z = 1, ..., its 99.9th percentile
#' rounded up, renormalised) plus the delay from death to recording; one
#' counted on its day of death (the fraction f, at its prior mean) after the
#' delay from death to recording alone. The delay from death to recording is
#' prompt (gamma, at the prior means of its centres; nowcast_gamma_days) for
#' a share 1 - b of deaths; the share b in the backlog are recorded from
#' hold_days after the death, with probability 1 - exp(-rate) each day (b and
#' rate at the centres of their priors). Truncated to D days and
#' renormalised.
#'
#' @param D the longest delay
#' @return the log probabilities of delays 0 to D
nowcast_reference_death <- function(D) {
  priors <- mechanistic_delay_priors()
  m <- priors$onset_to_death$mean
  s <- priors$onset_to_death$sd
  z_max <- ceiling(stats::qgamma(0.999, m^2 / s^2, m / s^2))
  p_die <- diff(stats::pgamma(0:z_max, m^2 / s^2, m / s^2))
  p_die <- p_die / sum(p_die)
  p_prompt <- nowcast_gamma_days(D, priors$recording$mean_mu, priors$recording$sd_mu)
  hold <- mechanistic_settings()$recording_backlog$hold_days
  rate <- priors$backlog$rate
  p_backlog <- ifelse(0:D >= hold, (1 - exp(-rate)) * exp(-rate * pmax(0, 0:D - hold)), 0)
  p_rec <- (1 - priors$backlog$share) * p_prompt + priors$backlog$share * p_backlog / sum(p_backlog)
  p_onset <- vapply(0:D, function(d) {
    z <- seq_len(min(d, z_max))
    sum(p_die[z] * p_rec[d - z + 1])
  }, numeric(1))
  fraction <- mechanistic_settings()$dated_by_death$fraction
  f <- fraction$alpha / (fraction$alpha + fraction$beta)
  p <- (1 - f) * p_onset / sum(p_onset) + f * p_rec
  log(p / sum(p))
}

#' Priors of each series (Table tab:nowcast-priors of the paper)
#'
#' The reference delays are those at the prior centres of the mechanistic
#' model (mechanistic_delay_priors): for cases, a gamma distribution with the
#' prior means of the centres of the mean and sd of the delay to a case being
#' reported; for deaths, the delay from onset to a death being recorded,
#' with fatal cases counted on their day of onset or of death
#' (nowcast_reference_death).
#' Cases have a tail beyond d_c = 28 days common to all onset days, with the
#' probability q of a record beyond it Beta(2, 38); deaths have no such tail.
#'
#' @param column "cases_reported" or "deaths_reported"
#' @return a named list
nowcast_priors <- function(column) {
  is_deaths <- column == "deaths_reported"
  list(sigma_lambda_scale = 0.3,
       # the delay of deaths changed more over the outbreak (the backlog)
       tau_t_scale = if(is_deaths) 0.5 else 0.2,
       sigma_u_scale = 0.5, rho_scale = 0.25,
       beta_g_scale = 3, log_phi_mu = log(20), log_phi_sd = 1.5,
       d_c = 28L, q_alpha = 2, q_beta = 38)
}

#' The data of one series
#'
#' The basis in delay is cubic in log(1 + d), with interior knots at quantiles
#' of the delays at which records were added (nowcast_delay_knots); the basis
#' in onset time is cubic, with interior knots every onset_knot_days days.
#'
#' @inheritParams nowcast_triangle
#' @param onset_knot_days the days of onset time between knots of h
#' @return a named list
nowcast_series_data <- function(snapshots, column, onset_knot_days = 14) {
  priors <- nowcast_priors(column)
  tri <- nowcast_triangle(snapshots, column)
  n_days <- max(tri$t)
  stopifnot(identical(sort(unique(tri$t)), seq_len(n_days)))
  delay_max <- max(tri$hi)
  reports <- sort(unique(tri$time_reported))
  delay_knots <- nowcast_delay_knots(tri)
  Bd <- splines::bs(log1p(0:delay_max), knots = log1p(delay_knots), degree = 3,
                    intercept = TRUE, Boundary.knots = c(0, log1p(delay_max)))
  onset_knots <- seq(1 + onset_knot_days, n_days - 1, by = onset_knot_days)
  Bt <- splines::bs(seq_len(n_days), knots = onset_knots, degree = 3, intercept = TRUE,
                    Boundary.knots = c(1, n_days))
  last <- tri %>% dplyr::group_by(.data$t) %>%
    dplyr::filter(.data$time_reported == max(.data$time_reported)) %>%
    dplyr::ungroup() %>% dplyr::arrange(.data$t)
  c(list(
    T = n_days, D = delay_max, Kt = ncol(Bt), Kd = ncol(Bd),
    Bt = unclass(Bt), Bd = unclass(Bd),
    R = length(reports), N = nrow(tri), t = tri$t, lo = tri$lo, hi = tri$hi,
    r = match(tri$time_reported, reports), report_day = as.integer(reports), n = tri$n,
    # the removals summed over, enough for removals of up to about half of
    # the additions
    w_max = pmax(0L, -tri$n) + 6L + as.integer(ceiling(0.5 * abs(tri$n))),
    count_last = last$count, delay_last = last$hi,
    log_p_ref = if(column == "deaths_reported") nowcast_reference_death(delay_max) else
      nowcast_reference_delay(delay_max, mechanistic_delay_priors()$reporting$mean_mu,
                              mechanistic_delay_priors()$reporting$sd_mu),
    delay_knots = delay_knots, onset_knots = onset_knots),
    priors)
}

#' A variant of the settings, for tuning them by the scores of the nowcasts
#' (set by the environment variable BVD_DD_VARIANT)
#'
#' @param name 'base' or "" (knots of the basis in onset time every 14 days
#'   and the half-normal prior of the sd of the case fatality ratio's walk with
#'   scale 0.05), 'knots7' and 'knots21' (knots every 7 or 21 days), 'cfr02'
#'   and 'cfr15' (the scale 0.02 or 0.15)
#' @return a list of onset_knot_days and sigma_cfr_scale
data_driven_variant <- function(name) {
  base <- list(onset_knot_days = 14, sigma_cfr_scale = 0.05)
  variants <- list(base = list(), knots7 = list(onset_knot_days = 7),
                   knots21 = list(onset_knot_days = 21),
                   cfr02 = list(sigma_cfr_scale = 0.02), cfr15 = list(sigma_cfr_scale = 0.15))
  if(name == "")
    name <- "base"
  if(!name %in% names(variants))
    stop(paste("Unknown variant:", name))
  utils::modifyList(base, variants[[name]])
}

#' The data for stan/nowcast_joint.stan
#'
#' Those of each series, with the suffixes _c and _d, and the onset days,
#' longest delay, reports and basis in onset time shared. Deaths have no
#' common tail (d_c set to the longest delay). The prior of the case fatality
#' ratio: its logit on the first onset day is N(0, 1.5^2), and it follows a
#' random walk over onset days whose sd has a half-normal prior of scale 0.05.
#'
#' @param analysis as prepare_bvd_reports
#' @param day_target the day (counted as the onset days are) to which counts
#'   are projected, or 0 for none
#' @param grainsize of reduce_sum
#' @return a named list, with elements 'cases' and 'deaths' holding the data
#'   of each series and 'stan' the data passed to Stan
nowcast_stan_data <- function(analysis, day_target = 0L, grainsize = 64L) {
  variant <- data_driven_variant(Sys.getenv("BVD_DD_VARIANT", ""))
  cases <- nowcast_series_data(analysis$snapshots, "cases_reported", variant$onset_knot_days)
  deaths <- nowcast_series_data(analysis$snapshots, "deaths_reported", variant$onset_knot_days)
  deaths$d_c <- deaths$D
  stopifnot(cases$T == deaths$T, cases$D == deaths$D, cases$R == deaths$R,
            isTRUE(all.equal(cases$Bt, deaths$Bt)))
  shared <- c("T", "D", "Kt", "Bt", "R")
  own <- c("Kd", "Bd", "N", "t", "lo", "hi", "r", "n", "w_max", "count_last",
           "delay_last", "log_p_ref", "d_c", "q_alpha", "q_beta", "beta_g_scale",
           "log_phi_mu", "log_phi_sd", "tau_t_scale", "sigma_u_scale", "rho_scale")
  f_own <- function(x, suffix) stats::setNames(x[own], paste0(own, suffix))
  stan <- c(cases[shared], list(day_target = as.integer(day_target)),
            f_own(cases, "_c"), list(sigma_lambda_scale_c = cases$sigma_lambda_scale),
            f_own(deaths, "_d"),
            list(sigma_cfr_scale = variant$sigma_cfr_scale,
                 grainsize = as.integer(grainsize)))
  list(stan = stan, cases = cases, deaths = deaths)
}

#' Initial values
#'
#' Stan's default (uniform on (-2, 2) unconstrained) can start the walk of
#' the log expected counts with a sd of several per day, which summed over
#' the onset days overflows; so the chains start near the level of the first
#' counts, with the scales at half their prior scales, small z, and the case
#' fatality ratio at that of all the deaths to all the cases.
#'
#' @param data as nowcast_stan_data
#' @param chains the number of chains
#' @param seed for the jitter
#' @return a list of lists, one for each chain
nowcast_inits <- function(data, chains, seed = 1) {
  f_series <- function(x, seed) withr::with_seed(seed, lapply(seq_len(chains), function(k) {
    f_z <- function(n) stats::rnorm(n, 0, 0.1)
    first <- x$n[x$t <= 7 & x$lo == 0]
    list(log_lambda_1 = log(mean(c(first, 1)) + 1), z_lambda = f_z(x$T - 1),
         sigma_lambda = x$sigma_lambda_scale / 2, beta_g_free = f_z(x$Kd - 1),
         z_h = matrix(f_z(x$Kt * (x$Kd - 1)), x$Kt), tau_t = x$tau_t_scale / 2,
         z_u = f_z(x$R), sigma_u = x$sigma_u_scale / 2, log_phi = x$log_phi_mu, rho = 0.05,
         q = x$q_alpha / (x$q_alpha + x$q_beta))
  }))
  cases <- f_series(data$cases, seed)
  deaths <- f_series(data$deaths, seed + 1)
  own <- c("beta_g_free", "z_h", "tau_t", "z_u", "sigma_u", "log_phi", "rho")
  ratio <- sum(data$deaths$count_last) / sum(data$cases$count_last)
  withr::with_seed(seed, lapply(seq_len(chains), function(k) {
    out <- c(cases[[k]][c("log_lambda_1", "z_lambda", "sigma_lambda")],
             stats::setNames(cases[[k]][own], paste0(own, "_c")),
             stats::setNames(deaths[[k]][own], paste0(own, "_d")),
             list(logit_cfr_1 = stats::qlogis(ratio), z_cfr = stats::rnorm(data$stan$T - 1, 0, 0.1),
                  sigma_cfr = data$stan$sigma_cfr_scale / 2))
    if(data$stan$d_c_c < data$stan$D)
      out$q_c <- array(cases[[k]]$q)
    out
  }))
}

#' The model, compiled with threads for reduce_sum
#'
#' @param model_file the Stan file
#' @return a CmdStanModel
nowcast_model <- function(model_file) {
  cmdstanr::cmdstan_model(model_file, cpp_options = list(stan_threads = TRUE))
}

#' Fits one chain and saves its draws
#'
#' The chains of a fit share their seed and differ in chain_id (so in their
#' random numbers) and in their initial values (those of nowcast_inits for
#' that chain).
#'
#' @param data as nowcast_stan_data
#' @param chain the chain (1 to chains)
#' @param path the file to which the draws are saved
#' @param iter_warmup,iter_sampling,adapt_delta,seed as cmdstanr's sample
#' @param chains the number of chains of the fit
#' @param threads the threads of reduce_sum
#' @param model_file the Stan file
#' @return path
fit_nowcast_chain <- function(data, chain, path, iter_warmup, iter_sampling, chains = 4,
                              threads = 4, adapt_delta = 0.95, seed = 1,
                              model_file = "stan/nowcast_joint.stan") {
  model <- nowcast_model(model_file)
  time_start <- Sys.time()
  fit <- model$sample(data = data$stan, chains = 1, chain_ids = chain,
                      threads_per_chain = threads, iter_warmup = iter_warmup,
                      iter_sampling = iter_sampling, adapt_delta = adapt_delta,
                      init = nowcast_inits(data, chains, seed)[chain], seed = seed,
                      refresh = 100, show_messages = FALSE)
  if(any(fit$return_codes() != 0))
    stop("The chain failed: ", paste(fit$return_codes(), collapse = " "))
  # the parameters with at least one element (q_d has none: deaths have no
  # common tail)
  present <- unique(sub("\\[.*", "", posterior::variables(fit$draws())))
  parameters <- intersect(names(model$variables()$parameters), present)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(chain = chain, minutes = as.numeric(Sys.time() - time_start, units = "mins"),
               draws = fit$draws(c("lp__", parameters)),
               sampler = fit$sampler_diagnostics()), path)
  path
}

#' The chains of a fit, joined
#'
#' @param paths the files of fit_nowcast_chain, one for each chain
#' @return a named list: 'draws' (a draws_array of the parameters, chains in
#'   order), 'sampler' (the sampler diagnostics) and 'minutes' (of each chain)
nowcast_join_chains <- function(paths) {
  chains <- lapply(paths, readRDS)
  chains <- chains[order(vapply(chains, `[[`, numeric(1), "chain"))]
  list(draws = posterior::bind_draws(lapply(chains, `[[`, "draws"), along = "chain"),
       sampler = posterior::bind_draws(lapply(chains, `[[`, "sampler"), along = "chain"),
       minutes = vapply(chains, `[[`, numeric(1), "minutes"))
}

#' Whether a fit has converged
#'
#' The R-hat and bulk and tail effective sample sizes (Vehtari et al. 2021)
#' of every parameter, and the divergent transitions. A fit is taken to have
#' converged when every R-hat is below 1.01, every effective sample size
#' above 400 and there are no divergent transitions.
#'
#' @param fit as nowcast_join_chains
#' @param date_now the fit date
#' @return a one-row tibble
nowcast_convergence <- function(fit, date_now) {
  s <- posterior::summarise_draws(fit$draws, "rhat", "ess_bulk", "ess_tail")
  sampler <- posterior::as_draws_df(fit$sampler)
  worst <- s$variable[which.max(s$rhat)]
  dplyr::tibble(
    date_now = as.Date(date_now), chains = posterior::nchains(fit$draws),
    draws_per_chain = posterior::niterations(fit$draws),
    minutes_max = max(fit$minutes),
    divergences = sum(sampler$divergent__), max_treedepth = sum(sampler$treedepth__ >= 10),
    rhat_max = max(s$rhat, na.rm = TRUE), worst = worst,
    ess_bulk_min = min(s$ess_bulk, na.rm = TRUE), ess_tail_min = min(s$ess_tail, na.rm = TRUE),
    converged = .data$rhat_max < 1.01 & .data$ess_bulk_min > 400 & .data$ess_tail_min > 400 &
      .data$divergences == 0)
}

#' Summaries of the parameters of a fit
#'
#' @inheritParams nowcast_convergence
#' @return a tibble of the scalar parameters, by series, with their quantiles
nowcast_parameters <- function(fit, date_now) {
  scalars <- c("sigma_lambda", "sigma_cfr", paste0(rep(c("tau_t", "sigma_u", "log_phi", "rho"), 2),
                                                    rep(c("_c", "_d"), each = 4)), "q_c[1]")
  draws <- posterior::subset_draws(fit$draws, intersect(scalars, posterior::variables(fit$draws)))
  draws <- posterior::mutate_variables(draws, phi_c = exp(log_phi_c), phi_d = exp(log_phi_d))
  posterior::summarise_draws(draws, "mean", "median", "sd", ~quantile(.x, c(0.05, 0.95)), "rhat",
                             "ess_bulk", "ess_tail") %>%
    dplyr::filter(!grepl("^log_phi", .data$variable)) %>%
    dplyr::mutate(column = dplyr::case_when(grepl("_d(\\[|$)", .data$variable) ~ "deaths_reported",
                                            .data$variable %in% c("sigma_cfr") ~ "deaths_reported",
                                            TRUE ~ "cases_reported"),
                  variable = sub("_(c|d)(\\[|$)", "\\2", .data$variable),
                  date_now = as.Date(date_now)) %>%
    dplyr::rename(q5 = "5%", q95 = "95%")
}

#' Projects the counts of a fit to a later day, from its saved draws
#'
#' Standalone generated quantities: for each draw, the reports still to come
#' are simulated up to day_target, each with the effect of a future report.
#'
#' @param fit as nowcast_join_chains
#' @param data as nowcast_stan_data (of the fit)
#' @param day_target the day (counted as the onset days are) projected to
#' @param seed of the simulation
#' @param model_file the Stan file
#' @return a named list of draws matrices (draws by onset day): projected_c,
#'   projected_d, mean_delay_c and mean_delay_d
nowcast_project <- function(fit, data, day_target, seed = 1,
                            model_file = "stan/nowcast_joint.stan") {
  model <- nowcast_model(model_file)
  stan <- data$stan
  stan$day_target <- as.integer(day_target)
  gq <- model$generate_quantities(fitted_params = posterior::subset_draws(
    fit$draws, setdiff(posterior::variables(fit$draws), "lp__")),
    data = stan, seed = seed, threads_per_chain = 1)
  f <- function(name) unclass(posterior::as_draws_matrix(gq$draws(name)))
  list(projected_c = f("projected_c"), projected_d = f("projected_d"),
       mean_delay_c = f("mean_delay_c"), mean_delay_d = f("mean_delay_d"))
}

#' The day of a date, counted as the onset days are
#'
#' @param analysis as prepare_bvd_reports
#' @param date the date
#' @return an integer
nowcast_day <- function(analysis, date) {
  as.integer(as.Date(date) - analysis$dates$date[1]) + 1L
}

#' The counts of a fit's onset days shown by a later report
#'
#' @param snapshots as build_bvd_snapshots
#' @param analysis the analysis of the fit (as prepare_bvd_reports)
#' @param date the later report
#' @return a tibble of time_onset, cases and deaths: the last count of each
#'   onset day of the fit shown by that report
nowcast_shown <- function(snapshots, analysis, date) {
  prepare_bvd_reports(snapshots, date)$snapshots %>%
    dplyr::filter(.data$time_onset <= nrow(analysis$dates)) %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::filter(.data$time_reported == max(.data$time_reported)) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(time_onset = .data$time_onset, cases = .data$cases_reported,
                     deaths = .data$deaths_reported)
}

#' Projections of a fit to each later report, against what the report showed
#'
#' @param fit as nowcast_join_chains
#' @param data,analysis of the fit
#' @param snapshots as build_bvd_snapshots
#' @param date_last the last report
#' @return a named list: 'totals' (by later report and series, the quantiles
#'   of the projected total of the fit's onset days and the total shown) and
#'   'by_day' (by later report, series and onset day, the quantiles of the
#'   projected count, the count shown, and the count of the fit)
nowcast_vantage <- function(fit, data, analysis, snapshots, date_last) {
  dates <- sort(unique(snapshots$date_reported))
  dates <- dates[dates > analysis$date_now & dates <= as.Date(date_last)]
  out <- lapply(seq_along(dates), function(i) {
    proj <- nowcast_project(fit, data, nowcast_day(analysis, dates[i]), seed = i)
    shown <- nowcast_shown(snapshots, analysis, dates[i])
    by_series <- lapply(c(cases = "c", deaths = "d"), function(s) {
      m <- proj[[paste0("projected_", s)]]
      series <- if(s == "c") "cases" else "deaths"
      counts <- dplyr::tibble(time_onset = seq_len(ncol(m))) %>%
        dplyr::left_join(dplyr::select(shown, "time_onset", shown = dplyr::all_of(series)),
                         by = "time_onset")
      q <- function(x, p) apply(x, 2, stats::quantile, p, names = FALSE)
      tq <- stats::quantile(rowSums(m), c(0.05, 0.5, 0.95), names = FALSE)
      count_last <- data[[series]]$count_last
      list(totals = dplyr::tibble(series = series, date_vantage = dates[i],
                                  p05 = tq[1], p50 = tq[2], p95 = tq[3],
                                  shown = sum(counts$shown, na.rm = TRUE),
                                  fitted = sum(count_last)),
           by_day = dplyr::tibble(series = series, date_vantage = dates[i],
                                  time_onset = counts$time_onset,
                                  date = analysis$dates$date[counts$time_onset],
                                  now = count_last,
                                  p05 = q(m, 0.05), p50 = q(m, 0.5), p95 = q(m, 0.95),
                                  shown = counts$shown))
    })
    list(totals = dplyr::bind_rows(lapply(by_series, `[[`, "totals")),
         by_day = dplyr::bind_rows(lapply(by_series, `[[`, "by_day")))
  })
  list(date_now = analysis$date_now,
       totals = dplyr::bind_rows(lapply(out, `[[`, "totals")) %>%
         dplyr::mutate(date_now = analysis$date_now),
       by_day = dplyr::bind_rows(lapply(out, `[[`, "by_day")) %>%
         dplyr::mutate(date_now = analysis$date_now))
}

#' The eventual counts of a fit, until every delay has elapsed
#'
#' The projection to a day a year after the last report fitted, beyond the
#' longest delay of any onset day. The case fatality ratio of an onset week
#' is the ratio of its projected deaths to its projected cases, in each draw.
#' Weeks whose last day is fewer than min_age days before the fit date are
#' marked as recent: most of their deaths are still to come.
#'
#' @param fit,data,analysis as nowcast_vantage
#' @param min_age the days from the end of a week to the fit date below which
#'   it is recent
#' @return a named list: 'by_day' (by series and onset day, the quantiles of
#'   the eventual count and of the mean delay), 'totals' and 'cfr' (by week)
nowcast_eventual <- function(fit, data, analysis, min_age = 21) {
  proj <- nowcast_project(fit, data, nowcast_day(analysis, analysis$date_now + 365))
  q <- function(x, p) apply(x, 2, stats::quantile, p, names = FALSE)
  dates <- analysis$dates$date[seq_len(data$stan$T)]
  by_day <- dplyr::bind_rows(lapply(c("c", "d"), function(s) {
    m <- proj[[paste0("projected_", s)]]
    d <- proj[[paste0("mean_delay_", s)]]
    series <- if(s == "c") "cases" else "deaths"
    count_last <- data[[series]]$count_last
    dplyr::tibble(series = series, time_onset = seq_len(ncol(m)), date = dates,
                  now = count_last,
                  p05 = q(m, 0.05), p50 = q(m, 0.5), p95 = q(m, 0.95),
                  mean_delay = q(d, 0.5), mean_delay_p05 = q(d, 0.05),
                  mean_delay_p95 = q(d, 0.95))
  })) %>% dplyr::mutate(date_now = analysis$date_now)
  totals <- dplyr::bind_rows(lapply(c(cases = "c", deaths = "d"), function(s) {
    tq <- stats::quantile(rowSums(proj[[paste0("projected_", s)]]), c(0.05, 0.5, 0.95), names = FALSE)
    dplyr::tibble(p05 = tq[1], p50 = tq[2], p95 = tq[3])
  }), .id = "series") %>% dplyr::mutate(date_now = analysis$date_now)
  week <- as.Date(cut(dates, "week"))
  f_week <- function(m) {
    w <- t(rowsum(t(m), week))
    colnames(w) <- as.character(sort(unique(week)))
    w
  }
  ratio <- f_week(proj$projected_d) / pmax(f_week(proj$projected_c), 1)
  weeks <- as.Date(colnames(ratio))
  cfr <- dplyr::tibble(date_now = analysis$date_now, week = weeks,
                       p05 = q(ratio, 0.05), p50 = q(ratio, 0.5), p95 = q(ratio, 0.95),
                       recent = as.numeric(analysis$date_now - (weeks + 6)) < min_age)
  list(by_day = by_day, totals = totals, cfr = cfr)
}

#' Scores the projections of a fit to the last report against what it showed
#'
#' @param vantage as nowcast_vantage
#' @param date_last the last report
#' @return a tibble by series: the projected total and its interval, the
#'   total shown, whether it is inside the 90% interval, and the coverage of
#'   the 90% intervals of the onset days (all, and the last 28 days)
nowcast_score <- function(vantage, date_last) {
  by_day <- dplyr::filter(vantage$by_day, .data$date_vantage == as.Date(date_last)) %>%
    dplyr::mutate(inside = .data$shown >= .data$p05 & .data$shown <= .data$p95,
                  recent = .data$time_onset > max(.data$time_onset) - 28)
  dplyr::filter(vantage$totals, .data$date_vantage == as.Date(date_last)) %>%
    dplyr::left_join(by_day %>% dplyr::group_by(.data$series) %>%
                       dplyr::summarise(coverage_all = mean(.data$inside, na.rm = TRUE),
                                        coverage_recent = mean(.data$inside[.data$recent], na.rm = TRUE)),
                     by = "series") %>%
    dplyr::mutate(inside = .data$shown >= .data$p05 & .data$shown <= .data$p95)
}

#' The number of columns of panels in a plot: two for more than three panels
#'
#' @param panel the panel of each row
#' @return 1 or 2
nowcast_columns <- function(panel) {
  if(dplyr::n_distinct(panel) > 3) 2 else 1
}

#' Saves a plot with a panel for each of several values
#'
#' @param path the file
#' @param g the plot
#' @param panel the panel of each row
nowcast_save <- function(path, g, panel) {
  ncol <- nowcast_columns(panel)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, g, width = c(9, 15)[ncol],
                  height = 2.6 * ceiling(dplyr::n_distinct(panel) / ncol) + 1)
}

f_fit_label <- function(date) {
  dates <- sort(unique(as.Date(date)))
  factor(format(as.Date(date), "Fit to %d %b"), levels = format(dates, "Fit to %d %b"))
}

#' Plots the projections of each fit to the last report against what it showed
#'
#' @param vantage a list of results of nowcast_vantage
#' @param date_last the last report
#' @param out the folder
#' @return the paths of the files written
plot_nowcast_backtests <- function(vantage, date_last, out) {
  by_day <- dplyr::bind_rows(lapply(vantage, `[[`, "by_day")) %>%
    dplyr::filter(.data$date_vantage == as.Date(date_last)) %>%
    dplyr::mutate(fit = f_fit_label(.data$date_now))
  vapply(c("cases", "deaths"), function(s) {
    df <- dplyr::filter(by_day, .data$series == s)
    g <- ggplot2::ggplot(df, ggplot2::aes(.data$date)) +
      ggplot2::geom_col(ggplot2::aes(y = .data$now, fill = "Shown in the last report fitted"),
                        width = 1) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$p05, ymax = .data$p95,
                                        fill = "Projected to the last report (90% interval)"),
                           alpha = 0.6) +
      ggplot2::geom_line(ggplot2::aes(y = .data$p50, colour = "Projected (median)")) +
      ggplot2::geom_line(ggplot2::aes(y = .data$shown, colour = "Shown in the last report")) +
      ggplot2::facet_wrap(~fit, ncol = nowcast_columns(df$fit), scales = "free_y") +
      ggplot2::scale_fill_manual(NULL, values = c("#a6c8f0", "grey80")) +
      ggplot2::scale_colour_manual(NULL, values = c("#1f4e9c", "black")) +
      ggplot2::labs(x = "Onset date", y = s) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
    path <- file.path(out, paste0("backtest_", s, ".png"))
    nowcast_save(path, g, df$fit)
    path
  }, character(1))
}

#' Plots, for each fit, the projected total of its onset days at each later
#' report against the total that report showed
#'
#' @inheritParams plot_nowcast_backtests
#' @return the paths of the files written
plot_nowcast_vantage <- function(vantage, out) {
  totals <- dplyr::bind_rows(lapply(vantage, `[[`, "totals")) %>%
    dplyr::mutate(fit = f_fit_label(.data$date_now))
  vapply(c("cases", "deaths"), function(s) {
    df <- dplyr::filter(totals, .data$series == s)
    g <- ggplot2::ggplot(df, ggplot2::aes(.data$date_vantage)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$fitted, linetype = "Shown in the last report fitted"),
                          colour = "grey50") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$p05, ymax = .data$p95,
                                        fill = "Projected (90% interval)"), alpha = 0.6) +
      ggplot2::geom_line(ggplot2::aes(y = .data$p50, colour = "Projected (median)")) +
      ggplot2::geom_point(ggplot2::aes(y = .data$shown, shape = "Shown in the report"), size = 1.4) +
      ggplot2::facet_wrap(~fit, ncol = nowcast_columns(df$fit), scales = "free") +
      ggplot2::scale_fill_manual(NULL, values = "#a6c8f0") +
      ggplot2::scale_colour_manual(NULL, values = "#1f4e9c") +
      ggplot2::scale_shape_manual(NULL, values = 16) +
      ggplot2::scale_linetype_manual(NULL, values = "dashed") +
      ggplot2::labs(x = "Date of the later report",
                    y = paste("Total", s, "of the onset days fitted")) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
    path <- file.path(out, paste0("vantage_", s, ".png"))
    nowcast_save(path, g, df$fit)
    path
  }, character(1))
}

#' Plots, for one fit, the projected count of each onset day at a few later
#' reports against what each showed
#'
#' @param vantage as nowcast_vantage, for one fit
#' @param weeks the later reports shown: the first report at least this many
#'   weeks after the fit date, for each value
#' @param days_shown the last days of onset shown
#' @param out the folder
#' @return the paths of the files written
plot_nowcast_vantage_days <- function(vantage, out, weeks = c(1, 2, 4, 6), days_shown = 56) {
  available <- sort(unique(vantage$by_day$date_vantage))
  chosen <- unique(vapply(weeks, function(w) {
    later <- available[available >= vantage$date_now + 7 * w]
    if(length(later) == 0) NA_real_ else as.numeric(min(later))
  }, numeric(1)))
  chosen <- as.Date(chosen[!is.na(chosen)], origin = "1970-01-01")
  df <- dplyr::filter(vantage$by_day, .data$date_vantage %in% chosen,
                      .data$date > vantage$date_now - days_shown) %>%
    dplyr::mutate(panel = factor(format(.data$date_vantage, "Report of %d %b"),
                                 levels = format(chosen, "Report of %d %b")))
  vapply(c("cases", "deaths"), function(s) {
    d <- dplyr::filter(df, .data$series == s)
    g <- ggplot2::ggplot(d, ggplot2::aes(.data$date)) +
      ggplot2::geom_col(ggplot2::aes(y = .data$now, fill = "Shown in the last report fitted"), width = 1) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$p05, ymax = .data$p95,
                                        fill = "Projected to the report (90% interval)"), alpha = 0.6) +
      ggplot2::geom_line(ggplot2::aes(y = .data$p50, colour = "Projected (median)")) +
      ggplot2::geom_line(ggplot2::aes(y = .data$shown, colour = "Shown in the report")) +
      ggplot2::facet_wrap(~panel, ncol = 2) +
      ggplot2::scale_fill_manual(NULL, values = c("#a6c8f0", "grey80")) +
      ggplot2::scale_colour_manual(NULL, values = c("#1f4e9c", "black")) +
      ggplot2::labs(x = "Onset date", y = s,
                    title = format(vantage$date_now, "Fit to %d %B")) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
    path <- file.path(out, paste0("vantage_days_", format(vantage$date_now, "%m%d"), "_", s, ".png"))
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(path, g, width = 12, height = 2.6 * ceiling(length(chosen) / 2) + 1.5)
    path
  }, character(1))
}

#' Plots the case fatality ratio by onset week of each fit, from its eventual
#' counts
#'
#' @param eventual a list of results of nowcast_eventual
#' @param out the folder
#' @return the path of the file written
plot_nowcast_cfr <- function(eventual, out) {
  df <- dplyr::bind_rows(lapply(eventual, `[[`, "cfr")) %>%
    dplyr::mutate(fit = f_fit_label(.data$date_now),
                  weeks = ifelse(.data$recent, "Recent (under 3 weeks old)", "Older"))
  g <- ggplot2::ggplot(df, ggplot2::aes(.data$week)) +
    ggplot2::geom_pointrange(ggplot2::aes(y = .data$p50, ymin = .data$p05, ymax = .data$p95,
                                          colour = .data$weeks), size = 0.25) +
    ggplot2::facet_wrap(~fit, ncol = nowcast_columns(df$fit)) +
    ggplot2::scale_colour_manual(NULL, values = c("Older" = "#1f4e9c",
                                                  "Recent (under 3 weeks old)" = "#e6550d")) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Week of onset", y = "Case fatality ratio (median, 90% interval)") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  path <- file.path(out, "cfr.png")
  nowcast_save(path, g, df$fit)
  path
}

#' Writes a table of what the fits infer, in meaningful units: the mean
#' delay from onset to a record for onset days through the outbreak, and the
#' other parameters re-expressed as quantities which can be read directly
#'
#' Each row is a quantity for a series (cases, deaths, or both when the
#' quantity joins them), and each column a fit date, with the posterior median
#' and 90% interval. The re-expressed parameters are monotone functions of a
#' single parameter, so their quantiles are those of the parameter carried
#' through. The mean delay is that of the onset day on the 15th of each month,
#' omitted for onset days within 14 days of the fit date, whose delay is
#' mostly carried forward from earlier days.
#'
#' @param pars the rows of nowcast_parameters for each fit
#' @param by_day the rows of by_day of nowcast_eventual for each fit
#' @param out the folder
#' @return the paths of the files written (csv and markdown)
write_nowcast_parameters <- function(pars, by_day, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  f_num <- function(x) sub("\\.$", "", formatC(signif(x, 2), digits = 2, format = "fg", flag = "#"))
  f_cell <- function(m, lo, hi) paste0(f_num(m), " (", f_num(pmin(lo, hi)), "-", f_num(pmax(lo, hi)), ")")
  delays <- by_day %>%
    dplyr::mutate(date = as.Date(.data$date), date_now = as.Date(.data$date_now)) %>%
    dplyr::filter(format(.data$date, "%d") == "15", .data$date <= .data$date_now - 14) %>%
    dplyr::transmute(group = "Delay", order = as.numeric(.data$date),
                     quantity = paste("Mean delay from onset to record (days), onset", format(.data$date, "%d %b")),
                     series = .data$series, date_now = .data$date_now,
                     value = f_cell(.data$mean_delay, .data$mean_delay_p05, .data$mean_delay_p95))
  inv_logit <- function(x) 1 / (1 + exp(-x))
  re <- list(
    sigma_lambda = list("Factor by which expected cases can change over a week (95% of weeks)",
                        function(s) exp(1.96 * s * sqrt(7))),
    sigma_cfr = list("Case fatality ratio of 0.5 after four weeks, upper end of 95% range",
                     function(s) inv_logit(1.96 * s * sqrt(28))),
    `q[1]` = list("Cases recorded more than 28 days after onset (%)", function(q) 100 * q),
    sigma_u = list("Factor by which a report can add more (or fewer) records than expected (95% of reports)",
                   function(s) exp(1.96 * s)),
    phi = list("Variance over mean of the records added, when 50 are expected",
               function(phi) 1 + 50 / phi),
    rho = list("Records removed for each 100 added", function(rho) 100 * rho))
  shared <- c("sigma_lambda", "sigma_cfr")
  others <- pars %>%
    dplyr::filter(.data$variable %in% names(re)) %>%
    dplyr::mutate(series = ifelse(.data$variable %in% shared, "cases and deaths",
                                  sub("_reported$", "", .data$column)),
                  date_now = as.Date(.data$date_now)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(group = "Other", order = match(.data$variable, names(re)),
                  quantity = re[[.data$variable]][[1]],
                  value = {f <- re[[.data$variable]][[2]]; f_cell(f(.data$median), f(.data$q5), f(.data$q95))}) %>%
    dplyr::ungroup()
  wide <- dplyr::bind_rows(delays, dplyr::select(others, dplyr::all_of(names(delays)))) %>%
    dplyr::mutate(series = factor(.data$series, c("cases and deaths", "cases", "deaths")),
                  group = factor(.data$group, c("Delay", "Other")), fit = format(.data$date_now, "%d %b")) %>%
    dplyr::arrange(.data$group, .data$series, .data$order, .data$date_now) %>%
    dplyr::select("quantity", "series", "fit", "value") %>%
    tidyr::pivot_wider(names_from = "fit", values_from = "value")
  paths <- file.path(out, c("parameters_table.csv", "parameters_table.md"))
  utils::write.csv(wide, paths[1], row.names = FALSE)
  md <- as.matrix(dplyr::mutate(wide, dplyr::across(dplyr::everything(), ~ifelse(is.na(.x), "", as.character(.x)))))
  writeLines(c("Posterior median (90% interval), for each fit date.", "",
               paste("|", paste(colnames(md), collapse = " | "), "|"),
               paste0("|", paste(rep("---", ncol(md)), collapse = "|"), "|"),
               apply(md, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))), paths[2])
  paths
}
