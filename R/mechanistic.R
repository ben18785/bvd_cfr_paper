# The mechanistic model of cases and deaths (the package incidenceinflation,
# version 2.0): the true cases of each onset day arise by renewal, are
# reported after a delay which changes at changepoints, and a share of them
# die (the case fatality ratio). A death is recorded after the time from
# onset to death and a delay which is estimated; a fraction of fatal cases are
# counted on their day of death, not of onset (the exact form); a death
# brings about the report of its case if it has not already been reported
# (the share of such fatal cases held at 0.95); and deaths recorded after a
# long delay (the backlog) have a delay of their own. The fits sample all of
# this by MCMC, one chain at a time, carried on across runs since a chain can
# need more than the six hours a GitHub job allows.

#' The reports up to a date in the form the mechanistic model takes
#'
#' The reports of prepare_bvd_reports, with what the sampler also needs: the
#' weekday of the first onset day, the cases of the onset days before it, and
#' the onset days over which the delay to a case being reported cannot change
#' (those first seen long after onset).
#'
#' @inheritParams prepare_bvd_reports
#' @param serial_max the days before the first onset day whose cases are given
#' @param days_observed onset days within this many days before the first
#'   report are taken as seen by it from their onset
#' @param delay_min_length the shortest span of onset days with one reporting
#'   delay
#' @return a named list as prepare_bvd_reports, with weekday_first,
#'   cases_before, dates_before, date_observed and reporting_first_length
mechanistic_analysis <- function(snapshots, date_now, date_first = "2026-05-11",
                                 serial_max = 40, days_observed = 7,
                                 delay_min_length = 14) {
  analysis <- prepare_bvd_reports(snapshots, date_now, date_first)
  date_start <- as.Date(date_first)
  f_time <- function(date) as.integer(as.Date(date) - date_start) + 1L
  analysis$weekday_first <- as.integer(format(date_start, "%u"))
  known <- snapshots %>%
    dplyr::filter(.data$date_reported <= analysis$date_now) %>%
    dplyr::filter(.data$date_reported == max(.data$date_reported),
                  .data$date_onset < date_start,
                  .data$date_onset >= date_start - serial_max)
  before <- dplyr::tibble(date = seq(date_start - serial_max, date_start - 1, by = 1)) %>%
    dplyr::left_join(dplyr::select(known, date = "date_onset", "cases_reported"), by = "date") %>%
    dplyr::mutate(cases = dplyr::coalesce(as.numeric(.data$cases_reported), 0))
  analysis$cases_before <- before$cases
  analysis$dates_before <- before$date
  analysis$date_observed <- min(snapshots$date_reported) - days_observed
  analysis$reporting_first_length <- max(
    delay_min_length, f_time(analysis$date_observed) - 1L + delay_min_length)
  analysis
}

#' The priors of the delays, shared with the data-driven model
#'
#' The centres of the gamma delay to a case being reported: its mean has a
#' gamma prior with mean 16 and sd 6.5 days, between the mean delay from
#' onset to notification of the outbreaks in the DRC of 1995 to 2014 (9.5
#' days; Rosello et al. 2015) and that from onset to a situation report in
#' North Kivu and Ituri in September 2018 (23 days; Tariq et al. 2019); its
#' sd has a gamma prior with mean 13 and sd 6.5 days, the sd of those delays
#' being about 0.8 of their mean. The centres of the delay from death to a
#' death being recorded have gamma priors with mean 4 and sd 4 days; the time
#' from onset to death is gamma with mean 9.3 and sd 5.5 days, from past
#' outbreaks. The backlog's share and rate are the centres of the package's
#' priors on them (logit share ~ N(logit 0.2, 1), log rate ~ N(log 1/60, 1)),
#' which the fits use.
#'
#' @return a named list: 'reporting' and 'recording' (mean_mu, mean_sigma,
#'   sd_mu and sd_sigma: the mean and sd of the gamma priors on the centres
#'   of the mean and sd of each delay), 'onset_to_death' (mean, sd) and
#'   'backlog' (share and rate: the centres of their priors)
mechanistic_delay_priors <- function() {
  list(reporting = list(mean_mu = 16, mean_sigma = 6.5, sd_mu = 13, sd_sigma = 6.5),
       recording = list(mean_mu = 4, mean_sigma = 4, sd_mu = 4, sd_sigma = 4),
       onset_to_death = list(mean = 9.3, sd = 5.5),
       backlog = list(share = 0.2, rate = 1 / 60))
}

#' The settings of the model of deaths
#'
#' The time from onset to death is gamma with mean 9.3 and sd 5.5 days (from
#' past outbreaks); the fraction of fatal cases counted on their day of death
#' has a uniform prior; the share of fatal cases reported when their death is
#' recorded is held at 0.95; the delay to a death being recorded is the same
#' for all onset days up to 31 July; and deaths are held in the backlog for at
#' least 21 days.
#'
#' @param ridge whether the sampler also makes the joint move of the fraction
#'   counted on the day of death and the delay to a death being recorded
#' @return a list, as the 'settings' of the deaths of
#'   incidenceinflation::mcmc_single
mechanistic_settings <- function(ridge = FALSE) {
  list(
    onset_to_death = mechanistic_delay_priors()$onset_to_death,
    dated_by_death = c(list(fraction = list(alpha = 1, beta = 1, initial = 0.25), exact = TRUE),
                       if(ridge) list(ridge = TRUE)),
    initial_delay = list(mean = 3, sd = 2.5),
    delay_prior = mechanistic_delay_priors()$recording,
    triggered = list(initial = 0.95, n_steps = 0L),
    fatal_delay = FALSE,
    recording_fixed_until = "2026-07-31",
    recording_backlog = list(hold_days = 21L))
}

#' The starting values of a chain: near one another, scaled by the chain
#' @param chain the chain (1, 2, ...)
mechanistic_starts <- function(chain) {
  scale <- c(1, 0.9, 1.1, 0.95)[(chain - 1) %% 4 + 1]
  list(inflation = 2 * scale, inflation_days = 4, cfr = 0.37 * scale,
       delay = list(mean = 15 * scale, sd = 13 * scale),
       reporting = list(mean = 7.6 * scale, sd = 5.4 * scale),
       overdispersion = 30 * scale)
}

#' Samples, or carries on, one chain of the mechanistic model
#'
#' If the state file exists (the progress of an earlier run of the chain), the
#' chain is carried on from it. The chain stops when it has its iterations or
#' after time_limit minutes, whichever is first, and saves the draws of the
#' second half of the iterations done (thinned to n_draws) to path and what is
#' needed to carry it on to the state file.
#'
#' @param analysis as mechanistic_analysis
#' @param chain the chain (1, 2, ...)
#' @param path the file of the chain's draws
#' @param state the file of the chain's progress
#' @param iterations the iterations the chain is to have in all
#' @param time_limit minutes after which the chain stops, or NULL
#' @param ridge as mechanistic_settings
#' @param n_draws draws kept
#' @return path
fit_mechanistic_chain <- function(analysis, chain, path, state, iterations,
                                  time_limit = NULL, ridge = FALSE, n_draws = 300) {
  set.seed(1000 * chain + as.integer(analysis$date_now) %% 1000)
  resume <- if(file.exists(state)) readRDS(state)
  snapshots <- analysis$snapshots
  largest <- snapshots %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::summarise(cases_true = max(.data$cases_reported))
  days_before_end <- max(snapshots$time_onset) - largest$time_onset
  start <- mechanistic_starts(chain)
  initial_cases_true <- largest %>%
    dplyr::mutate(cases_true = ceiling(
      .data$cases_true * (1 + start$inflation * exp(-days_before_end / start$inflation_days)) +
        5 * (days_before_end < 5)))
  # the largest fall in a count between reports, which erroneous records
  # have to explain
  f_largest_fall <- function(column) {
    falls <- snapshots %>%
      dplyr::group_by(.data$time_onset) %>%
      dplyr::arrange(.data$time_reported, .by_group = TRUE) %>%
      dplyr::summarise(fall = max(c(0, -diff(.data[[column]]))))
    max(falls$fall)
  }
  max_erroneous <- max(3, f_largest_fall("cases_reported") + 1,
                       f_largest_fall("deaths_reported") + 1)
  delay_changepoints <- list(min_length = 14, n_partition_moves = 50)
  priors <- list(
    Rt = list(shape = 1, rate = 1),
    reporting = mechanistic_delay_priors()$reporting,
    max_cases = ceiling(4 * max(initial_cases_true$cases_true) + 50),
    overdispersion = list(mean = 20, sd = 20),
    revisions = list(eta_mu = 0.005, eta_sigma = 0.005, rho_alpha = 1, rho_beta = 4),
    max_erroneous = max_erroneous,
    reporting_effects = list(tail_weight_alpha = 1, tail_weight_beta = 19,
                             tail_mean_mu = 40, tail_mean_sigma = 30),
    Rt_changepoints = list(),
    reporting_changepoints = c(delay_changepoints,
                               list(first_length = analysis$reporting_first_length)))
  deaths_settings <- utils::modifyList(
    list(initial_cfr = start$cfr, initial_delay = start$delay,
         delay_changepoints = delay_changepoints),
    mechanistic_settings(ridge))
  # no change in the delay to a death being recorded among the onset days up
  # to recording_fixed_until
  deaths_settings$delay_changepoints$first_length <- min(
    as.integer(as.Date(deaths_settings$recording_fixed_until) - analysis$dates$date[1]) + 1L,
    nrow(analysis$dates))
  deaths_settings$recording_fixed_until <- NULL

  done <- if(is.null(resume)) 0 else resume$iterations
  minutes_before <- if(is.null(resume)) 0 else resume$minutes
  time_start <- Sys.time()
  draws <- resume$draws
  state_sampler <- resume$state
  if(done < iterations) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    progress_file <- file.path(dirname(path), sub("\\.rds$", "_progress.txt", basename(path)))
    fit <- incidenceinflation::mcmc_single(
      niterations = iterations - done,
      snapshot_with_Rt_index_df = snapshots %>%
        dplyr::select("time_onset", "time_reported", "cases_reported") %>%
        dplyr::mutate(Rt_index = 1),
      priors = priors,
      serial_parameters = list(mean = 15.3, sd = 9.3),
      initial_cases_true = initial_cases_true,
      initial_reporting_parameters = start$reporting,
      initial_Rt = dplyr::tibble(Rt_index = 1, Rt = 1),
      reporting_metropolis_parameters = list(mean_step = 0.3, sd_step = 0.25),
      print_to_screen = FALSE,
      initial_overdispersion = start$overdispersion,
      is_negative_binomial = TRUE,
      overdispersion_metropolis_sd = 8,
      initial_revision_parameters = list(eta = 0.01, rho = 0.2, weekday = TRUE,
                                         weekday_first = analysis$weekday_first),
      revision_metropolis_parameters = list(eta_step = 0.3, rho_step = 0.3, pi_step = 0.3),
      cases_before = analysis$cases_before,
      n_parameter_steps = 5,
      progress_file = progress_file, progress_every = 10,
      deaths = list(snapshot = snapshots %>%
                      dplyr::select("time_onset", "time_reported", "deaths_reported"),
                    settings = deaths_settings),
      return_state = TRUE, resume_from = state_sampler, time_limit = time_limit)
    names_draws <- setdiff(names(fit), c("state", "iterations_done"))
    draws <- lapply(stats::setNames(names_draws, names_draws), function(name)
      dplyr::bind_rows(draws[[name]], fit[[name]]))
    state_sampler <- fit$state
    done <- fit$state$iterations
  }
  minutes <- minutes_before + as.numeric(Sys.time() - time_start, units = "mins")
  dir.create(dirname(state), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(draws = draws, state = state_sampler, iterations = done, minutes = minutes),
          state, compress = "xz")
  kept <- unique(round(seq(floor(done / 2) + 1, done, length.out = n_draws)))
  out <- list(chain = chain, date_now = analysis$date_now, analysis = analysis,
              ridge = ridge, niterations = done, niterations_target = iterations,
              minutes = minutes, max_erroneous = max_erroneous)
  for(name in intersect(names(draws), c(
    "cases", "Rt", "reporting_days", "overdispersion", "deaths", "deaths_parameters",
    "Rt_changepoints", "reporting_changepoints", "revisions", "reporting_effects")))
    out[[name]] <- draws[[name]] %>%
      dplyr::filter(.data$iteration %in% kept) %>%
      dplyr::select(-dplyr::any_of("chain")) %>%
      dplyr::mutate(chain = chain)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, path)
  path
}

#' The saved chains of a fit
#' @param paths the files of fit_mechanistic_chain
#' @return a list of the chains which exist
mechanistic_join_chains <- function(paths) {
  lapply(paths[file.exists(paths)], readRDS)
}

#' The draws of each chain of a fit as a draws array (iterations by chains by
#' variables): the scalar parameters, and the total cases and deaths of the
#' fit's onset days
#'
#' Only the draws which all chains have are kept (chains which have not
#' finished have fewer), the last ones of each.
#'
#' @param chains as mechanistic_join_chains
mechanistic_draws <- function(chains) {
  parts <- c("overdispersion", "revisions", "Rt_changepoints", "reporting_changepoints",
             "reporting_effects", "deaths_parameters")
  one <- function(x) {
    tables <- lapply(parts, function(p) {
      if(is.null(x[[p]])) return(NULL)
      df <- dplyr::select(x[[p]], -dplyr::any_of(c("chain", "variant")))
      df <- df[, vapply(df, is.numeric, logical(1)), drop = FALSE]
      names(df)[names(df) != "iteration"] <- paste0(p, ".", names(df)[names(df) != "iteration"])
      df
    })
    tables <- c(tables, list(
      x$cases %>% dplyr::group_by(.data$iteration) %>%
        dplyr::summarise(total.cases = sum(.data$cases_true)),
      x$deaths %>% dplyr::group_by(.data$iteration) %>%
        dplyr::summarise(total.deaths = sum(.data$deaths_true))))
    Reduce(function(a, b) dplyr::inner_join(a, b, by = "iteration"),
           Filter(Negate(is.null), tables)) %>%
      dplyr::arrange(.data$iteration)
  }
  each <- lapply(chains, one)
  n <- min(vapply(each, nrow, integer(1)))
  variables <- Reduce(intersect, lapply(each, names))
  variables <- setdiff(variables, "iteration")
  # the backlog's rate of each report day is not a parameter of interest
  variables <- variables[!grepl("backlog_rate_[0-9]", variables)]
  a <- array(NA_real_, c(n, length(each), length(variables)),
             dimnames = list(NULL, NULL, variables))
  for(k in seq_along(each))
    a[, k, ] <- as.matrix(utils::tail(each[[k]], n)[, variables])
  posterior::as_draws_array(a)
}

#' Checks of the convergence of a fit
#'
#' R-hat and the effective sample sizes of the scalar parameters and the
#' totals, from the draws of the second half of each chain.
#'
#' @param chains as mechanistic_join_chains
#' @param date_now the fit date
#' @return a one-row tibble
mechanistic_convergence <- function(chains, date_now) {
  draws <- mechanistic_draws(chains)
  s <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail") %>%
    dplyr::filter(is.finite(.data$rhat))
  # NA where nothing can be said (too few draws)
  f_total <- function(v, column) c(s[[column]][s$variable == v], NA)[1]
  f_ext <- function(f, x) if(any(is.finite(x))) f(x[is.finite(x)]) else NA_real_
  f_chains <- function(name) vapply(chains, function(x) as.numeric(x[[name]]), numeric(1))
  done <- f_chains("niterations")
  target <- f_chains("niterations_target")
  minutes <- f_chains("minutes")
  dplyr::tibble(
    date_now = as.Date(date_now), chains = length(done),
    iterations = paste(done, collapse = "/"), iterations_target = max(target),
    finished = all(done >= target), minutes_max = max(minutes),
    draws_per_chain = posterior::niterations(draws),
    rhat_max = f_ext(max, s$rhat), worst = c(s$variable[which.max(s$rhat)], NA)[1],
    rhat_cases = f_total("total.cases", "rhat"), rhat_deaths = f_total("total.deaths", "rhat"),
    ess_bulk_min = f_ext(min, s$ess_bulk), ess_tail_min = f_ext(min, s$ess_tail),
    converged = isTRUE(.data$rhat_max < 1.01 & .data$ess_bulk_min > 400 & .data$ess_tail_min > 400))
}

#' Summaries of the scalar parameters of a fit
#' @inheritParams mechanistic_convergence
#' @return a tibble of each parameter's part, name, quantiles and R-hat
mechanistic_parameters <- function(chains, date_now) {
  posterior::summarise_draws(mechanistic_draws(chains), "median",
                             ~quantile(.x, c(0.05, 0.95)), "rhat", "ess_bulk") %>%
    dplyr::rename(q5 = "5%", q95 = "95%") %>%
    tidyr::separate_wider_delim("variable", ".", names = c("part", "parameter")) %>%
    dplyr::mutate(date_now = as.Date(date_now))
}

#' The eventual cases and deaths of each onset day of a fit, and the delays
#' and case fatality ratio of each onset day
#'
#' @inheritParams mechanistic_convergence
#' @param snapshots as build_bvd_snapshots
#' @param date_last the last report, whose counts are set against the fit's
#' @return a named list: 'by_day' (by series and onset day, the quantiles of
#'   the eventual count, the count in the last report fitted and that in the
#'   last report), 'totals' (the same for the totals of the fit's onset days)
#'   and 'days' (by onset day, the quantiles of the case fatality ratio and of
#'   the mean delays to a case being reported and to a death being recorded)
mechanistic_eventual <- function(chains, snapshots, date_last) {
  analysis <- chains[[1]]$analysis
  date0 <- analysis$dates$date[1]
  f_last <- function(reports) reports %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::filter(.data$time_reported == max(.data$time_reported)) %>%
    dplyr::ungroup()
  now <- f_last(analysis$snapshots) %>%
    dplyr::select("time_onset", cases = "cases_reported", deaths = "deaths_reported") %>%
    tidyr::pivot_longer(-"time_onset", names_to = "series", values_to = "now")
  shown <- snapshots %>%
    dplyr::filter(.data$date_reported == as.Date(date_last)) %>%
    dplyr::transmute(time_onset = as.integer(.data$date_onset - date0) + 1L,
                     cases = .data$cases_reported, deaths = .data$deaths_reported) %>%
    tidyr::pivot_longer(-"time_onset", names_to = "series", values_to = "shown")
  draws <- dplyr::bind_rows(lapply(chains, function(x)
    dplyr::inner_join(
      dplyr::select(x$cases, "chain", "iteration", "time_onset", cases = "cases_true"),
      dplyr::select(x$deaths, "chain", "iteration", "time_onset", deaths = "deaths_true"),
      by = c("chain", "iteration", "time_onset")))) %>%
    tidyr::pivot_longer(c("cases", "deaths"), names_to = "series")
  f_q <- function(df) dplyr::summarise(df, p05 = stats::quantile(.data$value, 0.05),
                                       p50 = stats::median(.data$value),
                                       p95 = stats::quantile(.data$value, 0.95), .groups = "drop")
  by_day <- draws %>% dplyr::group_by(.data$series, .data$time_onset) %>% f_q() %>%
    dplyr::left_join(now, by = c("series", "time_onset")) %>%
    dplyr::left_join(shown, by = c("series", "time_onset")) %>%
    dplyr::mutate(date = date0 + .data$time_onset - 1, date_now = analysis$date_now)
  totals <- draws %>%
    dplyr::group_by(.data$series, .data$chain, .data$iteration) %>%
    dplyr::summarise(value = sum(.data$value), .groups = "drop") %>%
    dplyr::group_by(.data$series) %>% f_q() %>%
    dplyr::left_join(by_day %>% dplyr::group_by(.data$series) %>%
                       dplyr::summarise(now = sum(.data$now, na.rm = TRUE),
                                        shown = sum(.data$shown, na.rm = TRUE)),
                     by = "series") %>%
    dplyr::mutate(date_now = analysis$date_now)
  f_day <- function(df, name) df %>% dplyr::group_by(.data$time_onset) %>%
    dplyr::summarise(p05 = stats::quantile(.data$value, 0.05), p50 = stats::median(.data$value),
                     p95 = stats::quantile(.data$value, 0.95)) %>%
    dplyr::mutate(quantity = name)
  all_of <- function(part, column) dplyr::bind_rows(lapply(chains, function(x)
    dplyr::select(x[[part]], "time_onset", value = dplyr::all_of(column))))
  days <- dplyr::bind_rows(
    f_day(all_of("deaths", "cfr"), "Case fatality ratio"),
    f_day(all_of("reporting_days", "mean"), "Mean delay to a case being reported (days)"),
    f_day(all_of("deaths", "delay_mean"), "Mean delay to a death being recorded (days)"),
    f_day(all_of("Rt", "Rt"), "Rt")) %>%
    dplyr::mutate(date = date0 + .data$time_onset - 1, date_now = analysis$date_now)
  list(by_day = by_day, totals = totals, days = days)
}

f_fit_mechanistic <- function(date) {
  dates <- sort(unique(as.Date(date)))
  factor(format(as.Date(date), "Fit to %d %b"), levels = format(dates, "Fit to %d %b"))
}

#' Plots each fit's eventual cases and deaths of each onset day against those
#' shown in the last report
#'
#' @param eventual a list of results of mechanistic_eventual
#' @param out the folder
#' @param date_from the first onset day drawn
#' @return the paths of the files written
plot_mechanistic_backtests <- function(eventual, out, date_from = as.Date("2026-06-15")) {
  by_day <- dplyr::bind_rows(lapply(eventual, `[[`, "by_day")) %>%
    dplyr::filter(.data$date >= date_from) %>%
    dplyr::mutate(fit = f_fit_mechanistic(.data$date_now))
  vapply(c("cases", "deaths"), function(s) {
    df <- dplyr::filter(by_day, .data$series == s)
    g <- ggplot2::ggplot(df, ggplot2::aes(.data$date)) +
      ggplot2::geom_col(ggplot2::aes(y = .data$now, fill = "Shown in the last report fitted"),
                        width = 1) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$p05, ymax = .data$p95,
                                        fill = "Eventual (90% interval)"), alpha = 0.6) +
      ggplot2::geom_line(ggplot2::aes(y = .data$p50, colour = "Eventual (median)")) +
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

#' Plots each fit's eventual total of its onset days against the total shown
#' in the last report
#' @inheritParams plot_mechanistic_backtests
plot_mechanistic_totals <- function(eventual, out) {
  df <- dplyr::bind_rows(lapply(eventual, `[[`, "totals"))
  g <- ggplot2::ggplot(df, ggplot2::aes(.data$date_now)) +
    ggplot2::geom_pointrange(ggplot2::aes(y = .data$p50, ymin = .data$p05, ymax = .data$p95,
                                          colour = "Eventual (median, 90% interval)")) +
    ggplot2::geom_line(ggplot2::aes(y = .data$shown, colour = "Shown in the last report")) +
    ggplot2::facet_wrap(~series, scales = "free_y") +
    ggplot2::scale_colour_manual(NULL, values = c("#1f4e9c", "black")) +
    ggplot2::labs(x = "Fit date", y = "Total of the onset days up to the fit date") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  path <- file.path(out, "totals.png")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, g, width = 10, height = 4.5)
  path
}

#' Plots, for each fit, the case fatality ratio, the mean delays and Rt by
#' onset day
#' @inheritParams plot_mechanistic_backtests
plot_mechanistic_days <- function(eventual, out) {
  df <- dplyr::bind_rows(lapply(eventual, `[[`, "days")) %>%
    dplyr::mutate(fit = f_fit_mechanistic(.data$date_now))
  g <- ggplot2::ggplot(df, ggplot2::aes(.data$date, group = .data$fit)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$p05, ymax = .data$p95, fill = .data$fit),
                         alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = .data$p50, colour = .data$fit)) +
    ggplot2::facet_wrap(~quantity, ncol = 2, scales = "free_y") +
    ggplot2::labs(x = "Onset date", y = "Median (90% interval)", colour = NULL, fill = NULL) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  path <- file.path(out, "by_onset_day.png")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, g, width = 11, height = 7.5)
  path
}

#' Writes a table of what the fits infer, in meaningful units: the case
#' fatality ratio and the mean delays for onset days on the 15th of each
#' month (omitted within 14 days of the fit date), and the other parameters
#' which can be read directly
#'
#' @param pars the rows of mechanistic_parameters for each fit
#' @param days the rows of 'days' of mechanistic_eventual for each fit
#' @param out the folder
#' @return the paths of the files written (csv and markdown)
write_mechanistic_parameters <- function(pars, days, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  f_num <- function(x) sub("\\.$", "", formatC(signif(x, 2), digits = 2, format = "fg", flag = "#"))
  f_cell <- function(m, lo, hi) paste0(f_num(m), " (", f_num(lo), "-", f_num(hi), ")")
  by_day <- days %>%
    dplyr::filter(.data$quantity != "Rt", format(.data$date, "%d") == "15",
                  .data$date <= .data$date_now - 14) %>%
    dplyr::transmute(order = match(.data$quantity, unique(days$quantity)) * 1e6 + as.numeric(.data$date),
                     quantity = paste0(.data$quantity, ", onset ", format(.data$date, "%d %b")),
                     date_now = .data$date_now, value = f_cell(.data$p50, .data$p05, .data$p95))
  re <- list(
    deaths_parameters.dated_fraction = list("Fatal cases counted on their day of death (%)", function(x) 100 * x),
    deaths_parameters.wrong_dead = list("Erroneous records shown as dead (%)", function(x) 100 * x),
    deaths_parameters.n_delay_changes = list("Changes in the delay to a death being recorded", identity),
    deaths_parameters.n_cfr_changes = list("Changes in the case fatality ratio", identity),
    reporting_changepoints.n_changes = list("Changes in the delay to a case being reported", identity),
    Rt_changepoints.n_changes = list("Changes in Rt", identity),
    overdispersion.overdispersion = list("Variance over mean of the true cases, when 50 are expected",
                                         function(k) 1 + 50 / k),
    revisions.eta = list("Erroneous records added per true case (%)", function(x) 100 * x),
    revisions.rho = list("Chance an erroneous record is withdrawn at a report (%)", function(x) 100 * x))
  others <- pars %>%
    dplyr::mutate(key = paste(.data$part, .data$parameter, sep = ".")) %>%
    dplyr::filter(.data$key %in% names(re)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(order = 1e9 + match(.data$key, names(re)),
                  quantity = re[[.data$key]][[1]],
                  value = {
                    f <- re[[.data$key]][[2]]
                    v <- c(f(.data$median), sort(c(f(.data$q5), f(.data$q95))))
                    # the numbers of changes are whole numbers
                    cell <- if(grepl("n_(delay_|cfr_)?changes$", .data$key))
                      sprintf("%d (%d-%d)", round(v[1]), floor(v[2]), ceiling(v[3]))
                    else f_cell(v[1], v[2], v[3])
                    paste0(cell, " [", if(is.na(.data$rhat)) "-" else sprintf("%.2f", .data$rhat), "]")}) %>%
    dplyr::ungroup()
  wide <- dplyr::bind_rows(by_day, dplyr::select(others, dplyr::all_of(names(by_day)))) %>%
    dplyr::arrange(.data$order, .data$date_now) %>%
    dplyr::mutate(fit = format(.data$date_now, "%d %b")) %>%
    dplyr::select("quantity", "fit", "value") %>%
    tidyr::pivot_wider(names_from = "fit", values_from = "value")
  paths <- file.path(out, c("parameters_table.csv", "parameters_table.md"))
  utils::write.csv(wide, paths[1], row.names = FALSE)
  md <- as.matrix(dplyr::mutate(wide, dplyr::across(dplyr::everything(), ~ifelse(is.na(.x), "", .x))))
  writeLines(c("Posterior median (90% interval), for each fit date; [R-hat] for the parameters.", "",
               paste("|", paste(colnames(md), collapse = " | "), "|"),
               paste0("|", paste(rep("---", ncol(md)), collapse = "|"), "|"),
               apply(md, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))), paths[2])
  paths
}
